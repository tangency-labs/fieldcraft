//! Desktop native panel viewports, drag-out preview bridge, exit-safe redock, and
//! post-Step-7 cross-viewport synchronization (`#1088`).

use std::collections::{BTreeMap, BTreeSet};
use std::hash::{Hash, Hasher};

use egui::{Vec2, ViewportBuilder, ViewportClass, ViewportCommand, ViewportId, pos2};
use solidcraft_panel_layout::{
    Availability, ContentMeasureRequest, Destination, HostEffect, HostId, InstanceId, LayoutAction, MAX_LOGICAL, Mutation, PanelResult, Point, Rule,
    Section, SectionLocation, SideHeightPolicy, Size, Zone, error, next_section_id_with_reservations, presented_tab, tab_bounds,
};
use solidcraft_ui_kit::panel_runtime::{
    ContentServices, DesktopNativeAdapter, EmptyPanelHost, HostContentMetrics, NativeHostAdapter, NativeMonitor, NativeWindowAction,
    PanelLayoutStore, PanelRuntime, RuntimeActionSink, WindowCard, WorkspacePanelHost, native_window_chrome, paint_native_drop, paint_native_window,
    rule_key,
};
use solidcraft_ui_text::tr;

use super::content::{self, FieldCraftCommands, FieldCraftFrameContent, FieldCraftLocal, FieldCraftModels};
use super::store::FieldCraftLayoutStore;
use super::{MIN_PANEL_SIZE, PANEL_INSPECTOR, to_egui_rect};
use crate::ModelApp;

/// Maximum monitor count accepted from the desktop host.
const MAX_MONITORS: usize = 64;
/// Universal GPU 2D surface cap (physical pixels) for secondary native viewports.
const MAX_NATIVE_SURFACE_PX: f32 = 8192.0;
/// Floor on DPI scale used when converting physical caps into logical points.
const MIN_DPI_SCALE: f32 = 0.5;
/// Minimum usable native window dimensions (points).
const MIN_NATIVE_WIDTH: f32 = 360.0;
const MIN_NATIVE_HEIGHT: f32 = 440.0;
/// Default OS window frame allowance used when seeding saved `Zone::Window` origins.
const DEFAULT_WINDOW_DECORATION: Vec2 = Vec2::new(16.0, 48.0);
/// Minimum observed logical-point delta before updating native window geometry in the layout.
const GEOMETRY_EPSILON_PT: f32 = 0.5;
/// Root frames to defer `Prefs` persistence after an OS-`×` `NativeClose` so
/// child-before-root OS shutdown cannot wipe saved multi-monitor placements.
const NATIVE_CLOSE_PERSIST_DELAY_FRAMES: u8 = 2;

/// Native window adapter switching cleanly between in-window capability refusal and
/// desktop OS viewports.
pub(crate) enum NativeAdapter {
    InApp(solidcraft_ui_kit::panel_runtime::InAppNativeAdapter),
    Desktop(DesktopNativeAdapter),
}

impl NativeAdapter {
    pub(crate) fn desktop(&self) -> Option<&DesktopNativeAdapter> {
        match self {
            Self::Desktop(d) => Some(d),
            Self::InApp(_) => None,
        }
    }
}

impl NativeHostAdapter for NativeAdapter {
    fn capabilities(&self) -> solidcraft_panel_layout::HostCapabilities {
        match self {
            Self::InApp(a) => a.capabilities(),
            Self::Desktop(a) => a.capabilities(),
        }
    }

    fn prepare(&self, effects: &[HostEffect]) -> PanelResult<()> {
        match self {
            Self::InApp(a) => a.prepare(effects),
            Self::Desktop(a) => a.prepare(effects),
        }
    }

    fn commit(&mut self, effects: &[HostEffect]) {
        match self {
            Self::InApp(a) => a.commit(effects),
            Self::Desktop(a) => a.commit(effects),
        }
    }
}

/// Filter, bound, and deterministically sort physical-pixel monitor records
/// (primary first, then physical origin `(x, y)`, size `(w, h)`, and name).
pub fn canonicalize_monitors(monitors: impl IntoIterator<Item = NativeMonitor>) -> Vec<NativeMonitor> {
    let mut out: Vec<NativeMonitor> = monitors
        .into_iter()
        .filter(|m| m.bounds.is_finite() && m.bounds.is_positive() && (0.1..=16.0).contains(&m.scale))
        .take(MAX_MONITORS)
        .collect();
    out.sort_by(|a, b| {
        b.primary
            .cmp(&a.primary)
            .then_with(|| a.bounds.min.x.total_cmp(&b.bounds.min.x))
            .then_with(|| a.bounds.min.y.total_cmp(&b.bounds.min.y))
            .then_with(|| a.bounds.width().total_cmp(&b.bounds.width()))
            .then_with(|| a.bounds.height().total_cmp(&b.bounds.height()))
            .then_with(|| a.name.as_deref().unwrap_or("").cmp(b.name.as_deref().unwrap_or("")))
    });
    out
}

fn clone_monitors(monitors: &[NativeMonitor]) -> Vec<NativeMonitor> {
    monitors.iter().map(|m| NativeMonitor { name: m.name.clone(), bounds: m.bounds, scale: m.scale, primary: m.primary }).collect()
}

/// Fingerprint monitor topology for transition logging.
fn monitors_fingerprint(monitors: &[NativeMonitor]) -> u64 {
    let mut h = std::hash::DefaultHasher::new();
    monitors.len().hash(&mut h);
    for m in monitors {
        m.name.hash(&mut h);
        m.primary.hash(&mut h);
        m.bounds.min.x.to_bits().hash(&mut h);
        m.bounds.min.y.to_bits().hash(&mut h);
        m.bounds.width().to_bits().hash(&mut h);
        m.bounds.height().to_bits().hash(&mut h);
        m.scale.to_bits().hash(&mut h);
    }
    h.finish()
}

/// Whether a secondary viewport is minimized or hidden with zero/invalid area.
pub(crate) fn is_minimized_viewport(inner_size: Vec2, visible: Option<bool>) -> bool {
    visible == Some(false) || !inner_size.is_finite() || inner_size.x <= 0.0 || inner_size.y <= 0.0
}

/// Clamp native window preferred size and min/max bounds against physical GPU texture
/// limits and monitor physical bounds divided by DPI scale, guaranteeing `min <= max`.
pub(crate) fn clamp_window_size(preferred: Size, (min, max): (Size, Size), monitor_phys: Option<Vec2>, scale: f32) -> (Size, (Size, Size)) {
    let safe_scale = if scale.is_finite() { scale.max(MIN_DPI_SCALE) } else { 1.0 };
    let cap_w =
        (monitor_phys.filter(|m| m.is_finite() && m.x > 0.0).map_or(MAX_NATIVE_SURFACE_PX, |m| m.x.min(MAX_NATIVE_SURFACE_PX)) / safe_scale).max(1.0);
    let cap_h =
        (monitor_phys.filter(|m| m.is_finite() && m.y > 0.0).map_or(MAX_NATIVE_SURFACE_PX, |m| m.y.min(MAX_NATIVE_SURFACE_PX)) / safe_scale).max(1.0);
    let floor_w = MIN_PANEL_SIZE.min(cap_w);
    let floor_h = MIN_PANEL_SIZE.min(cap_h);

    let max_w = if max.width.is_finite() { max.width.clamp(floor_w, cap_w) } else { cap_w };
    let max_h = if max.height.is_finite() { max.height.clamp(floor_h, cap_h) } else { cap_h };
    let min_w = if min.width.is_finite() { min.width.clamp(1.0, max_w) } else { floor_w.min(max_w) };
    let min_h = if min.height.is_finite() { min.height.clamp(1.0, max_h) } else { floor_h.min(max_h) };

    let clamped_min = Size { width: min_w, height: min_h };
    let clamped_max = Size { width: max_w.max(min_w), height: max_h.max(min_h) };
    let pref_w = if preferred.width.is_finite() { preferred.width.clamp(clamped_min.width, clamped_max.width) } else { clamped_min.width };
    let pref_h = if preferred.height.is_finite() { preferred.height.clamp(clamped_min.height, clamped_max.height) } else { clamped_min.height };
    (Size { width: pref_w, height: pref_h }, (clamped_min, clamped_max))
}

fn preferred_window_size(size: Size, (min, max): (Size, Size)) -> Size {
    Size {
        width: size.width.max(MIN_NATIVE_WIDTH).clamp(min.width, max.width.max(min.width)),
        height: size.height.max(MIN_NATIVE_HEIGHT).clamp(min.height, max.height.max(min.height)),
    }
}

fn portrait(size: Size) -> Size {
    Size { height: size.height.max(size.width * 4.0 / 3.0), ..size }
}

fn content_height(measure: Option<&dyn HostContentMetrics>, instance: &InstanceId, width: f32, chrome: f32) -> Option<f32> {
    let request = ContentMeasureRequest { instance: instance.clone(), width, chrome_height: chrome };
    let measure = measure?.measure(&request).ok()?;
    (measure.policy == SideHeightPolicy::Content && measure.needed_height.is_finite() && measure.needed_height < MAX_LOGICAL)
        .then_some(measure.needed_height)
}

fn drop_window_size(zone: Zone, card: Size, header: f32, chrome: f32, side_width: f32, column_height: f32) -> Size {
    if zone == Zone::Bottom {
        Size { width: side_width, height: column_height }
    } else {
        Size { width: card.width, height: card.height - header.clamp(0.0, card.height) + chrome }
    }
}

#[derive(Default)]
pub(crate) struct WindowFrameResult {
    pub(crate) action: Option<NativeWindowAction>,
    pub(crate) embedded_fallback: bool,
    pub(crate) minimized_or_hidden: bool,
    pub(crate) focused: bool,
    pub(crate) observed_geometry: Option<(Option<Point>, Size)>,
}

/// Abstracts egui secondary-viewport execution for headless testing.
pub(crate) trait ViewportDriver {
    fn root_close_requested(&self, ctx: &egui::Context) -> bool;
    fn root_focused(&self, ctx: &egui::Context) -> Option<bool>;
    fn escape_pressed(&self, ctx: &egui::Context) -> bool;
    fn modal_input_blocked(&self, ctx: &egui::Context) -> bool;
    fn native_ready(&mut self, ctx: &egui::Context, desktop: &DesktopNativeAdapter, host: &HostId, builder: &ViewportBuilder) -> PanelResult<bool>;
    fn drop_ready(&mut self, ctx: &egui::Context, desktop: &DesktopNativeAdapter, builder: &ViewportBuilder) -> PanelResult<bool>;
    #[allow(clippy::too_many_arguments)]
    fn show_native_window(
        &mut self,
        ctx: &egui::Context,
        desktop: &DesktopNativeAdapter,
        viewport: ViewportId,
        builder: ViewportBuilder,
        host: &HostId,
        section: &Section,
        active: &InstanceId,
        runtime: &mut PanelRuntime,
        services: &mut ContentServices<'_>,
        app: &mut ModelApp,
        enabled: bool,
    ) -> PanelResult<WindowFrameResult>;
    fn show_drop_preview(
        &mut self,
        ctx: &egui::Context,
        viewport: ViewportId,
        builder: ViewportBuilder,
        card: Option<&WindowCard>,
        chrome: f32,
    ) -> bool;
    fn send_viewport_cmd_to(&mut self, ctx: &egui::Context, viewport: ViewportId, cmd: ViewportCommand);
    fn request_repaint(&mut self, ctx: &egui::Context);
    fn request_repaint_of(&mut self, ctx: &egui::Context, viewport: ViewportId);
}

pub(crate) struct EguiViewportDriver;

impl ViewportDriver for EguiViewportDriver {
    fn root_close_requested(&self, ctx: &egui::Context) -> bool {
        ctx.input(|i| i.viewport().close_requested())
    }

    fn root_focused(&self, ctx: &egui::Context) -> Option<bool> {
        ctx.input(|i| i.viewport().focused)
    }

    fn escape_pressed(&self, ctx: &egui::Context) -> bool {
        ctx.input(|i| i.key_pressed(egui::Key::Escape))
    }

    fn modal_input_blocked(&self, ctx: &egui::Context) -> bool {
        solidcraft_ui_kit::Container::modal_input_blocked(ctx)
    }

    fn native_ready(&mut self, ctx: &egui::Context, desktop: &DesktopNativeAdapter, host: &HostId, builder: &ViewportBuilder) -> PanelResult<bool> {
        desktop.ready(ctx, host, builder)
    }

    fn drop_ready(&mut self, ctx: &egui::Context, desktop: &DesktopNativeAdapter, builder: &ViewportBuilder) -> PanelResult<bool> {
        desktop.drop_ready(ctx, builder)
    }

    fn show_native_window(
        &mut self,
        ctx: &egui::Context,
        desktop: &DesktopNativeAdapter,
        viewport: ViewportId,
        builder: ViewportBuilder,
        host: &HostId,
        section: &Section,
        active: &InstanceId,
        runtime: &mut PanelRuntime,
        services: &mut ContentServices<'_>,
        app: &mut ModelApp,
        enabled: bool,
    ) -> PanelResult<WindowFrameResult> {
        ctx.show_viewport_immediate(viewport, builder, |ui, class| {
            if class == ViewportClass::EmbeddedWindow {
                return Ok(WindowFrameResult { embedded_fallback: true, ..Default::default() });
            }
            let close_requested = ui.input(|i| i.viewport().close_requested());
            if close_requested {
                return Ok(WindowFrameResult { action: Some(NativeWindowAction::Redock), ..Default::default() });
            }
            let (inner_size, visible, focused, observed_geometry) = ui.input(|i| {
                let vp = i.viewport();
                let size = vp.inner_rect.map_or_else(|| i.content_rect().size(), |r| r.size());
                let normal = vp.visible() != Some(false) && vp.maximized != Some(true) && vp.fullscreen != Some(true);
                let geom = normal
                    .then(|| {
                        let s = Size { width: size.x, height: size.y };
                        s.validate().ok()?;
                        let origin = vp
                            .outer_rect
                            .filter(|r| r.is_finite() && r.is_positive() && r.size().max_elem() <= MAX_LOGICAL)
                            .map(|r| Point { x: r.min.x, y: r.min.y })
                            .filter(|p| p.validate().is_ok());
                        Some((origin, s))
                    })
                    .flatten();
                (size, vp.visible(), vp.focused == Some(true), geom)
            });
            if is_minimized_viewport(inner_size, visible) {
                runtime.mount_native(host, &Default::default(), &mut EmptyPanelHost, services)?;
                return Ok(WindowFrameResult { minimized_or_hidden: true, focused, ..Default::default() });
            }
            desktop.observe(ui.ctx(), host, active);
            if enabled && focused && !ui.ctx().egui_wants_keyboard_input() {
                app.shortcuts(ui.ctx());
            }
            let mut content = FieldCraftFrameContent { app };
            let action = solidcraft_ui_registry::isolated(|| paint_native_window(runtime, services, &mut content, host, section, ui, enabled))?;
            Ok(WindowFrameResult { action, embedded_fallback: false, minimized_or_hidden: false, focused, observed_geometry })
        })
    }

    fn show_drop_preview(
        &mut self,
        ctx: &egui::Context,
        viewport: ViewportId,
        builder: ViewportBuilder,
        card: Option<&WindowCard>,
        chrome: f32,
    ) -> bool {
        ctx.show_viewport_immediate(viewport, builder, |ui, class| {
            paint_native_drop(ui, card, chrome);
            class != ViewportClass::EmbeddedWindow && ui.input(|i| i.viewport().native_pixels_per_point.is_some())
        })
    }

    fn send_viewport_cmd_to(&mut self, ctx: &egui::Context, viewport: ViewportId, cmd: ViewportCommand) {
        ctx.send_viewport_cmd_to(viewport, cmd);
    }

    fn request_repaint(&mut self, ctx: &egui::Context) {
        ctx.request_repaint();
    }

    fn request_repaint_of(&mut self, ctx: &egui::Context, viewport: ViewportId) {
        ctx.request_repaint_of(viewport);
    }
}

/// Coordinates desktop native OS panel windows, drag-out drop preview, exit-safe redock,
/// capability fallback, and cross-viewport repaint/focus synchronization.
#[derive(Default)]
pub(crate) struct NativeViewportCoordinator {
    monitors: Vec<NativeMonitor>,
    monitors_fingerprint: u64,
    global_coords_available: Option<bool>,
    embedded_downgraded: bool,
    exiting: bool,
    was_exit_in_flight: bool,
    pending_close_persist: Option<u8>,
    titles: BTreeMap<HostId, String>,
    active_viewports: BTreeSet<ViewportId>,
    seeded_hosts: BTreeSet<HostId>,
    any_child_focused: bool,
    last_epoch: Option<u64>,
}

impl NativeViewportCoordinator {
    pub(crate) fn update_monitors(
        &mut self,
        desktop: Option<&DesktopNativeAdapter>,
        monitors: impl FnOnce() -> Vec<NativeMonitor>,
        global_coords_available: bool,
    ) -> bool {
        let canonical = canonicalize_monitors(monitors());
        let fp = monitors_fingerprint(&canonical);
        let topology_changed = fp != self.monitors_fingerprint || self.global_coords_available != Some(global_coords_available);
        if topology_changed {
            log::info!("desktop monitor topology updated: {} monitor(s), global_coords={global_coords_available}", canonical.len());
            fieldcraft_engine::crash::activity("panel monitor topology updated");
            self.monitors_fingerprint = fp;
            self.global_coords_available = Some(global_coords_available);
        }
        self.monitors = canonical;
        if let Some(desktop) = desktop {
            let snapshot = clone_monitors(&self.monitors);
            desktop.update_monitors(|| snapshot);
        }
        !self.embedded_downgraded && global_coords_available && !self.monitors.is_empty()
    }

    pub(crate) fn mark_exiting(&mut self) {
        if !self.exiting {
            log::info!("desktop panel coordinator marked exiting");
            fieldcraft_engine::crash::activity("panel coordinator exiting");
        }
        self.exiting = true;
        self.pending_close_persist = None;
    }

    pub(crate) fn is_embedded_downgraded(&self) -> bool {
        self.embedded_downgraded
    }

    pub(crate) fn has_pending_close_persist(&self) -> bool {
        self.pending_close_persist.is_some()
    }

    pub(crate) fn clear_active_state(&mut self) {
        self.titles.clear();
        self.active_viewports.clear();
        self.seeded_hosts.clear();
        self.any_child_focused = false;
    }

    pub(crate) fn push_to_desktop(&self, desktop: &DesktopNativeAdapter) {
        let snapshot = clone_monitors(&self.monitors);
        desktop.update_monitors(|| snapshot);
    }

    pub(crate) fn exit_in_flight(&self, app: &ModelApp, root_close_requested: bool) -> bool {
        self.exiting || app.quit_requested || root_close_requested
    }

    fn monitor_bounds_and_scale(&self, ctx: &egui::Context, section: &Section) -> (Option<Vec2>, f32) {
        let target = match &section.location {
            SectionLocation::Window { origin, .. } => self.monitors.iter().find(|m| {
                let scale = (ctx.zoom_factor() * m.scale).max(MIN_DPI_SCALE);
                let center = pos2(origin.x + section.size.width * 0.5, origin.y + section.size.height * 0.5);
                let phys_center = (center.to_vec2() * scale).to_pos2();
                m.bounds.contains(phys_center)
            }),
            _ => None,
        }
        .or_else(|| self.monitors.iter().find(|m| m.primary))
        .or_else(|| self.monitors.first());
        match target {
            Some(m) => (Some(m.bounds.size()), ctx.zoom_factor() * m.scale),
            None => (None, ctx.pixels_per_point()),
        }
    }

    /// Advance 2-frame deferred `NativeClose` persistence or cancel if shutdown began.
    pub(crate) fn tick_close_persist(&mut self, store: &mut FieldCraftLayoutStore, app: &ModelApp, root_close_requested: bool) {
        let in_flight = self.exit_in_flight(app, root_close_requested);
        if in_flight {
            if self.pending_close_persist.take().is_some() {
                log::info!("cancelled pending NativeClose prefs write because exit is in flight");
                fieldcraft_engine::crash::activity("panel native close persist cancelled on exit");
            }
            self.was_exit_in_flight = true;
            return;
        }
        if self.was_exit_in_flight && !in_flight {
            self.was_exit_in_flight = false;
        }
        if let Some(remaining) = self.pending_close_persist {
            if remaining <= 1 {
                self.pending_close_persist = None;
                store.mark_dirty();
            } else {
                self.pending_close_persist = Some(remaining - 1);
            }
        }
    }

    /// Render active native OS windows and drag-out drop preview window.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn show_desktop_windows<D: ViewportDriver>(
        &mut self,
        driver: &mut D,
        app: &mut ModelApp,
        ctx: &egui::Context,
        runtime: &mut PanelRuntime,
        host: &mut WorkspacePanelHost,
        native: &mut NativeAdapter,
        store: &mut FieldCraftLayoutStore,
        local: &mut FieldCraftLocal,
    ) -> PanelResult<()> {
        let root_close_requested = driver.root_close_requested(ctx);
        self.tick_close_persist(store, app, root_close_requested);
        let Some(desktop) = native.desktop().cloned() else {
            return Ok(());
        };
        self.show_drop_window(driver, ctx, runtime, host, &desktop);

        let windows: Vec<(HostId, Section)> = runtime
            .layout()
            .layout()
            .sections
            .iter()
            .filter_map(|section| {
                let SectionLocation::Window { host, .. } = &section.location else { return None };
                Some((host.clone(), section.clone()))
            })
            .collect();

        if windows.is_empty() {
            self.active_viewports.clear();
            self.titles.clear();
            self.seeded_hosts.clear();
            self.any_child_focused = false;
            return Ok(());
        }

        self.push_to_desktop(&desktop);
        self.seeded_hosts.retain(|h| windows.iter().any(|(id, _)| id == h));

        let enabled = !driver.modal_input_blocked(ctx);
        let mut next_active = BTreeSet::new();
        let mut any_focused = false;
        let mut pending_geometry_updates: Vec<(HostId, Option<Point>, Size)> = Vec::new();
        let mut failure = None;

        for (host_id, section) in windows {
            let step = (|| -> PanelResult<()> {
                let layout = runtime.layout();
                let available: BTreeMap<InstanceId, Availability> = section.tabs.iter().map(|id| (id.clone(), layout.availability(id))).collect();
                let active = (!app.start_shown()).then(|| presented_tab(&section, &available)).flatten();
                let Some(active) = active else {
                    desktop.suspend(&host_id);
                    self.titles.remove(&host_id);
                    let models = FieldCraftModels;
                    let mut commands = FieldCraftCommands;
                    let mut services = ContentServices { models: &models, commands: &mut commands, local };
                    runtime.mount_native(&host_id, &Default::default(), &mut EmptyPanelHost, &mut services)?;
                    return Ok(());
                };

                let bounds = tab_bounds(layout.layout(), layout.registry(), &section)?;
                if !desktop.owns(&host_id, &section.id) {
                    return Err(error(Rule::Lifecycle, "native window lost its layout registration"));
                }
                let title =
                    if active.0 == PANEL_INSPECTOR { content::inspector_window_title(app) } else { runtime.instance_presentation(active)?.title };
                let (monitor_phys, scale) = self.monitor_bounds_and_scale(ctx, &section);
                let raw_pref = preferred_window_size(section.size, bounds);
                let (preferred, clamped_bounds) = clamp_window_size(raw_pref, bounds, monitor_phys, scale);
                if !self.seeded_hosts.contains(&host_id)
                    && let SectionLocation::Window { origin, .. } = &section.location
                    && origin.validate().is_ok()
                {
                    let outer =
                        egui::Rect::from_min_size(pos2(origin.x, origin.y), Vec2::new(preferred.width, preferred.height) + DEFAULT_WINDOW_DECORATION);
                    desktop.place(ctx, active, outer);
                    self.seeded_hosts.insert(host_id.clone());
                }
                let Some(builder) = desktop.builder(ctx, &host_id, active, title.clone(), preferred, clamped_bounds) else {
                    return Ok(());
                };
                let Some(viewport) = desktop.viewport_id(&host_id) else {
                    return Ok(());
                };

                match driver.native_ready(ctx, &desktop, &host_id, &builder) {
                    Ok(false) => return Ok(()),
                    Err(reason) => {
                        log::warn!("desktop panel window failed to open for {}: {reason}", host_id.0);
                        fieldcraft_engine::crash::activity("panel native open failed");
                        self.seeded_hosts.remove(&host_id);
                        self.recover_unopened_window(runtime, native, store, local, &host_id, &section)?;
                        app.set_status(tr(rule_key(reason.rule)).to_owned(), true);
                        return Err(reason);
                    }
                    Ok(true) => {}
                }

                if self.titles.get(&host_id) != Some(&title) {
                    self.titles.insert(host_id.clone(), title.clone());
                    driver.send_viewport_cmd_to(ctx, viewport, ViewportCommand::Title(title));
                }

                let models = FieldCraftModels;
                let mut commands = FieldCraftCommands;
                let mut services = ContentServices { models: &models, commands: &mut commands, local };
                let result =
                    driver.show_native_window(ctx, &desktop, viewport, builder, &host_id, &section, active, runtime, &mut services, app, enabled)?;

                if result.embedded_fallback {
                    log::warn!("EmbeddedWindow detected for native panel {}; downgrading session to Floating", host_id.0);
                    fieldcraft_engine::crash::activity("panel native embedded downgrade");
                    self.embedded_downgraded = true;
                    self.clear_active_state();
                    return Ok(());
                }

                if !result.minimized_or_hidden {
                    next_active.insert(viewport);
                    if result.focused {
                        any_focused = true;
                    }
                    if result.action.is_none()
                        && !self.exit_in_flight(app, root_close_requested)
                        && let Some((obs_origin, obs_size)) = result.observed_geometry
                    {
                        let clamped_size = Size {
                            width: obs_size.width.clamp(bounds.0.width, bounds.1.width),
                            height: obs_size.height.clamp(bounds.0.height, bounds.1.height),
                        };
                        let size_changed = (clamped_size.width - section.size.width).abs() > GEOMETRY_EPSILON_PT
                            || (clamped_size.height - section.size.height).abs() > GEOMETRY_EPSILON_PT;
                        let origin_changed = matches!(
                            (&section.location, obs_origin),
                            (SectionLocation::Window { origin, .. }, Some(o))
                                if (o.x - origin.x).abs() > GEOMETRY_EPSILON_PT || (o.y - origin.y).abs() > GEOMETRY_EPSILON_PT
                        );
                        if size_changed || origin_changed {
                            pending_geometry_updates.push((host_id.clone(), obs_origin, clamped_size));
                        }
                    }
                }

                let Some(action) = result.action else {
                    return Ok(());
                };

                let in_flight = self.exit_in_flight(app, root_close_requested);
                if matches!(action, NativeWindowAction::Redock) && in_flight {
                    log::info!("suppressed NativeClose redock for {} while exit is in flight", host_id.0);
                    return Ok(());
                }

                let close = matches!(action, NativeWindowAction::Redock);
                let pre_close_layout = close.then(|| runtime.layout().layout().clone());
                let registry = runtime.layout().registry().clone();
                let (mutation, persist_layout) = match action {
                    NativeWindowAction::Activate(instance) => (Mutation::Activate { instance }, true),
                    NativeWindowAction::Expand(section) => (Mutation::Collapse { section, collapsed: false }, true),
                    NativeWindowAction::Content(instance, action) => (Mutation::HeaderAction { instance, action }, false),
                    NativeWindowAction::Redock => (Mutation::NativeClose { host: host_id.clone() }, false),
                };

                let models = FieldCraftModels;
                let mut commands = FieldCraftCommands;
                let mut services = ContentServices { models: &models, commands: &mut commands, local };
                let action_req = LayoutAction::Apply { expected_revision: None, mutations: vec![mutation] };
                let dispatch_res = if persist_layout {
                    let mut sink = RuntimeActionSink { runtime, native, services: &mut services };
                    store.execute_action(action_req, &mut sink, &registry)
                } else {
                    runtime.dispatch(action_req, native, &mut services).inspect(|r| {
                        if close {
                            if let Some(prev) = &pre_close_layout {
                                store.preserve_before_deferred_close(prev, &registry);
                            }
                            let _ = store.replace_live(&r.layout, &registry, false);
                        }
                    })
                };

                match dispatch_res {
                    Ok(_) => {
                        if close {
                            log::info!("native panel host {} redocked by user", host_id.0);
                            fieldcraft_engine::crash::activity("panel native redock");
                            self.titles.remove(&host_id);
                            self.seeded_hosts.remove(&host_id);
                            next_active.remove(&viewport);
                            self.pending_close_persist = Some(NATIVE_CLOSE_PERSIST_DELAY_FRAMES);
                        }
                        driver.request_repaint(ctx);
                    }
                    Err(reason) => {
                        if close {
                            if reason.rule == Rule::NoSpace && self.recover_unopened_window(runtime, native, store, local, &host_id, &section).is_ok()
                            {
                                self.titles.remove(&host_id);
                                self.seeded_hosts.remove(&host_id);
                                next_active.remove(&viewport);
                                self.pending_close_persist = Some(NATIVE_CLOSE_PERSIST_DELAY_FRAMES);
                            } else {
                                driver.send_viewport_cmd_to(ctx, viewport, ViewportCommand::CancelClose);
                            }
                        }
                        return Err(reason);
                    }
                }
                Ok(())
            })();
            if self.embedded_downgraded {
                return Ok(());
            }
            if let Err(reason) = step {
                failure.get_or_insert(reason);
            }
        }

        if !pending_geometry_updates.is_empty() {
            let mut updated = runtime.layout().layout().clone();
            let mut updated_tabs = Vec::new();
            for (host_id, new_origin, new_size) in &pending_geometry_updates {
                if let Some(sec) =
                    updated.sections.iter_mut().find(|s| matches!(&s.location, SectionLocation::Window { host, .. } if host == host_id))
                {
                    sec.size = *new_size;
                    if let Some(o) = new_origin
                        && let SectionLocation::Window { origin, .. } = &mut sec.location
                    {
                        *origin = *o;
                    }
                    updated_tabs.push((sec.tabs.clone(), *new_origin, *new_size));
                }
            }
            for (tabs, new_origin, new_size) in updated_tabs {
                for entry in &mut updated.instances {
                    if tabs.contains(&entry.id) {
                        entry.memory.size = new_size;
                        if let Some(o) = new_origin
                            && let SectionLocation::Window { origin, .. } = &mut entry.memory.location
                        {
                            *origin = o;
                        }
                    }
                }
            }
            let registry = runtime.layout().registry().clone();
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local };
            if let Ok(reply) = runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout: updated }, native, &mut services) {
                let _ = store.replace_live(&reply.layout, &registry, true);
            }
        }

        self.active_viewports = next_active;
        self.any_child_focused = any_focused;
        failure.map_or(Ok(()), Err)
    }

    fn show_drop_window<D: ViewportDriver>(
        &mut self,
        driver: &mut D,
        ctx: &egui::Context,
        runtime: &mut PanelRuntime,
        host: &mut WorkspacePanelHost,
        desktop: &DesktopNativeAdapter,
    ) {
        if let Some((instance, outer)) = host.take_window_landing() {
            desktop.place(ctx, &instance, to_egui_rect(outer));
            let opened = runtime.layout().layout().containing(&instance).ok().and_then(|section| match &section.location {
                SectionLocation::Window { host, .. } => Some(host.clone()),
                _ => None,
            });
            if let Some(host_id) = opened
                && let Some(viewport) = desktop.adopt(&host_id, &instance)
            {
                self.seeded_hosts.insert(host_id.clone());
                log::info!("drag-out panel {} adopted into native host {}", instance.0, host_id.0);
                fieldcraft_engine::crash::activity("panel native drag-out adopted");
                driver.send_viewport_cmd_to(ctx, viewport, ViewportCommand::Focus);
            }
        }

        let Some((instance, corner)) = host.window_landing_preview() else {
            desktop.end_drop();
            return;
        };

        // Tear down preview on Escape or when focus leaves every FieldCraft viewport.
        let root_lost_focus = driver.root_focused(ctx) == Some(false) && !self.any_child_focused;
        if driver.escape_pressed(ctx) || root_lost_focus {
            host.cancel_gesture();
            desktop.end_drop();
            return;
        }

        let side = host.leaving_side_card().filter(|card| card.instance == instance).map(|card| (card.outer.size, card.header.size.height));
        let in_app = host.frame().and_then(|frame| {
            let card = frame.layout.sections.iter().find(|card| card.instance == instance)?;
            let measure = frame.layout.geometry_policy.measure.clone();
            Some((card.zone, card.outer.size, card.header.size.height, frame.layout.capacity.size.height, measure))
        });
        let mut bootstrap = false;
        let mut chrome = 0.0;
        let layout = runtime.layout();
        let landing = layout.layout().containing(&instance).ok().zip(layout.layout().instance(&instance).ok()).and_then(|(section, panel)| {
            let spec = layout.registry().get(&panel.kind).ok()?;
            let bounds = spec.size_bounds(Zone::Window);
            chrome = native_window_chrome(ctx, !spec.header_actions.is_empty());
            let size = match in_app {
                None => preferred_window_size(section.size, bounds),
                Some((zone, size, header, column, measure)) => {
                    let width = spec.default_size.width;
                    let fitted = (zone == Zone::Bottom).then(|| content_height(measure.as_deref(), &instance, width, chrome)).flatten();
                    bootstrap = zone == Zone::Bottom && fitted.is_none();
                    let size = match side {
                        Some((size, header)) => drop_window_size(Zone::Left, size, header, 0.0, width, column),
                        None => drop_window_size(zone, size, header, chrome, width, fitted.map_or(column, |needed| needed.min(column))),
                    };
                    if zone == Zone::Bottom { portrait(size) } else { size }
                }
            };
            let (monitor_phys, scale) = self.monitor_bounds_and_scale(ctx, section);
            let (size, clamped_bounds) = clamp_window_size(size, bounds, monitor_phys, scale);
            desktop.landing(ctx, &instance, pos2(corner.x, corner.y), size, clamped_bounds)
        });
        if bootstrap {
            driver.request_repaint(ctx);
        }
        let Some(outer) = landing else {
            return;
        };
        let title = runtime.instance_presentation(&instance).map(|p| p.title).unwrap_or_default();
        let (viewport, builder) = desktop.drop_window(ctx, &instance, title, outer);
        if !driver.drop_ready(ctx, desktop, &builder).unwrap_or(false) {
            return;
        }
        let shown = driver.show_drop_preview(ctx, viewport, builder, host.window_card(), if side.is_some() { 0.0 } else { chrome });
        if shown {
            host.confirm_window_shown(corner, super::logical_rect(outer));
        } else {
            driver.request_repaint(ctx);
        }
    }

    pub(crate) fn recover_unopened_window(
        &self,
        runtime: &mut PanelRuntime,
        native: &mut NativeAdapter,
        store: &mut FieldCraftLayoutStore,
        local: &mut FieldCraftLocal,
        host_id: &HostId,
        section: &Section,
    ) -> PanelResult<()> {
        let close = LayoutAction::Apply { expected_revision: None, mutations: vec![Mutation::NativeClose { host: host_id.clone() }] };
        let models = FieldCraftModels;
        let mut commands = FieldCraftCommands;
        let mut services = ContentServices { models: &models, commands: &mut commands, local };
        match runtime.dispatch(close, native, &mut services) {
            Ok(reply) => {
                let _ = store.replace_live(&reply.layout, runtime.layout().registry(), false);
                return Ok(());
            }
            Err(reason) if reason.rule != Rule::NoSpace => {}
            Err(_) => {}
        }

        let layout = runtime.layout().layout().clone();
        let reservations = store.preparation_reservations(&layout, runtime.layout().registry()).unwrap_or_default();
        let id = next_section_id_with_reservations(&layout, &reservations)?;
        let insert_index = layout.sections.len().saturating_sub(1);
        let mut mutations = Vec::new();
        for (index, instance) in section.tabs.iter().enumerate() {
            let destination = if index == 0 {
                Destination::NewSection {
                    section: id.clone(),
                    index: insert_index,
                    location: SectionLocation::Floating { origin: Point { x: 0.0, y: 0.0 } },
                    size: section.size,
                }
            } else {
                Destination::JoinTab { section: id.clone(), index }
            };
            mutations.push(Mutation::Move { instance: instance.clone(), destination });
        }
        mutations.push(Mutation::Activate { instance: section.active.clone() });
        if section.collapsed {
            mutations.push(Mutation::Collapse { section: id, collapsed: true });
        }
        let models = FieldCraftModels;
        let mut commands = FieldCraftCommands;
        let mut services = ContentServices { models: &models, commands: &mut commands, local };
        let reply = runtime.dispatch(LayoutAction::Apply { expected_revision: None, mutations }, native, &mut services)?;
        let _ = store.replace_live(&reply.layout, runtime.layout().registry(), false);
        Ok(())
    }

    /// Compute post-Step-7 mutation epoch and fan out repaints + root focus handoff across viewports.
    pub(crate) fn sync_post_step7<D: ViewportDriver>(
        &mut self,
        driver: &mut D,
        app: &ModelApp,
        ctx: &egui::Context,
        had_dialog: bool,
        had_refused: bool,
        had_drained_commands: bool,
    ) {
        let epoch = mutation_epoch(app, had_drained_commands);
        let changed = self.last_epoch.is_some_and(|prev| prev != epoch);
        self.last_epoch = Some(epoch);
        if changed && !self.active_viewports.is_empty() {
            driver.request_repaint_of(ctx, ViewportId::ROOT);
            for vp in &self.active_viewports {
                driver.request_repaint_of(ctx, *vp);
            }
        }
        let newly_opened_transient = (!had_dialog && app.dialog.is_some()) || (!had_refused && app.refused.is_some());
        if newly_opened_transient && !self.active_viewports.is_empty() {
            driver.send_viewport_cmd_to(ctx, ViewportId::ROOT, ViewportCommand::Focus);
        }
    }

    #[cfg(test)]
    pub(crate) fn active_viewports(&self) -> &BTreeSet<ViewportId> {
        &self.active_viewports
    }

    #[cfg(test)]
    pub(crate) fn pending_close_persist(&self) -> Option<u8> {
        self.pending_close_persist
    }
}

fn mutation_epoch(app: &ModelApp, had_drained_commands: bool) -> u64 {
    let mut h = std::hash::DefaultHasher::new();
    app.session.revision().hash(&mut h);
    app.revision.hash(&mut h);
    app.pick.generation.hash(&mut h);
    app.pick.hover_generation.hash(&mut h);
    app.probe.as_ref().map(|p| (p.node, p.step.as_str(), p.mode, p.frame)).hash(&mut h);
    app.path_editor.as_ref().map(|e| (e.step.as_str(), e.field.id(), e.points.len(), e.picking, e.selected, e.samples, e.frame, e.mode)).hash(&mut h);
    app.dialog.as_ref().map(|d| (d.command, d.opening, d.text.len(), d.error.is_some())).hash(&mut h);
    app.refused.as_ref().map(|u| u.command).hash(&mut h);
    app.ui.dark.hash(&mut h);
    app.ui.style.hash(&mut h);
    app.ui.show_nodes.hash(&mut h);
    app.ui.show_geometry.hash(&mut h);
    app.ui.show_mesh.hash(&mut h);
    app.ui.show_loads.hash(&mut h);
    app.ui.show_interactions.hash(&mut h);
    app.session.units().label().hash(&mut h);
    app.prefs.language.hash(&mut h);
    app.selected.as_ref().map(|(k, n)| (*k as u8, n.as_str())).hash(&mut h);
    app.selected_seed.is_some().hash(&mut h);
    app.selected_mesh.hash(&mut h);
    app.reply.as_ref().map(|r| r.command.as_str()).hash(&mut h);
    app.session.results_view.as_ref().map(|v| (v.step.as_str(), v.mode, v.deformed)).hash(&mut h);
    had_drained_commands.hash(&mut h);
    h.finish()
}

#[cfg(test)]
pub(crate) fn pop_out_panel(host: &mut crate::panels::PanelHost, panel_id: &str, origin: Point) -> PanelResult<()> {
    let instance = InstanceId::new(panel_id)?;
    let layout = host.runtime.layout().layout().clone();
    let current = layout.containing(&instance)?.clone();
    let spec = host.runtime.layout().registry().get(&layout.instance(&instance)?.kind)?;
    let reservations = host.store.preparation_reservations(&layout, host.runtime.layout().registry()).unwrap_or_default();
    let section_id = next_section_id_with_reservations(&layout, &reservations)?;
    let mut num = 0usize;
    let host_id = loop {
        let candidate = HostId::new(format!("host.window.{num}"))?;
        if !layout.sections.iter().any(|s| matches!(&s.location, SectionLocation::Window { host, .. } if *host == candidate)) {
            break candidate;
        }
        num += 1;
    };
    let dock_zone = if current.location.zone().docked() {
        current.location.zone()
    } else if spec.default_zone.docked() {
        spec.default_zone
    } else {
        Zone::Left
    };
    let return_index = layout.sections.iter().position(|s| s.id == current.id).unwrap_or(0);
    let insert_index = layout.sections.len().saturating_sub(usize::from(current.tabs.len() == 1));
    let destination = Destination::NewSection {
        section: section_id,
        index: insert_index,
        location: SectionLocation::Window {
            host: host_id,
            origin,
            return_to: solidcraft_panel_layout::DockReturn {
                zone: dock_zone,
                index: return_index,
                size: current.size,
                collapsed: current.collapsed,
                section: Some(current.id),
                column: None,
                column_index: 0,
                column_width: None,
                weight: 1.0,
                preferred_extent: None,
                extent_priority: 0,
                join: None,
            },
        },
        size: current.size,
    };
    host.apply_mutations(vec![Mutation::Move { instance, destination }], true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::files::Prefs;
    use crate::panels::{PANEL_BROWSER, PANEL_INSPECTOR, PANEL_RESULTS, PanelHost, panel_kind};
    use egui::{Rect, vec2};
    use fieldcraft_engine::Session;
    use serde_json::json;
    use solidcraft_panel_layout::{HostCapabilities, Mutation, PanelError};

    struct MockDriver {
        root_close_requested: bool,
        root_focused: Option<bool>,
        escape: bool,
        modal_blocked: bool,
        ready_result: Result<bool, PanelError>,
        drop_ready_result: Result<bool, PanelError>,
        next_window_result: WindowFrameResult,
        drop_shown_return: bool,
        sent_commands: Vec<(ViewportId, ViewportCommand)>,
        repaints: Vec<ViewportId>,
        root_repaints: usize,
        windows_shown: Vec<(HostId, InstanceId)>,
        last_builders: BTreeMap<HostId, ViewportBuilder>,
        drop_previews_shown: usize,
    }

    impl Default for MockDriver {
        fn default() -> Self {
            Self {
                root_close_requested: false,
                root_focused: Some(true),
                escape: false,
                modal_blocked: false,
                ready_result: Ok(true),
                drop_ready_result: Ok(true),
                next_window_result: WindowFrameResult::default(),
                drop_shown_return: true,
                sent_commands: Vec::new(),
                repaints: Vec::new(),
                root_repaints: 0,
                windows_shown: Vec::new(),
                last_builders: BTreeMap::new(),
                drop_previews_shown: 0,
            }
        }
    }

    impl ViewportDriver for MockDriver {
        fn root_close_requested(&self, _ctx: &egui::Context) -> bool {
            self.root_close_requested
        }
        fn root_focused(&self, _ctx: &egui::Context) -> Option<bool> {
            self.root_focused
        }
        fn escape_pressed(&self, _ctx: &egui::Context) -> bool {
            self.escape
        }
        fn modal_input_blocked(&self, _ctx: &egui::Context) -> bool {
            self.modal_blocked
        }
        fn native_ready(
            &mut self,
            _ctx: &egui::Context,
            _desktop: &DesktopNativeAdapter,
            _host: &HostId,
            _builder: &ViewportBuilder,
        ) -> PanelResult<bool> {
            self.ready_result.clone()
        }
        fn drop_ready(&mut self, _ctx: &egui::Context, _desktop: &DesktopNativeAdapter, _builder: &ViewportBuilder) -> PanelResult<bool> {
            self.drop_ready_result.clone()
        }
        fn show_native_window(
            &mut self,
            _ctx: &egui::Context,
            _desktop: &DesktopNativeAdapter,
            _viewport: ViewportId,
            builder: ViewportBuilder,
            host: &HostId,
            _section: &Section,
            active: &InstanceId,
            _runtime: &mut PanelRuntime,
            _services: &mut ContentServices<'_>,
            _app: &mut ModelApp,
            _enabled: bool,
        ) -> PanelResult<WindowFrameResult> {
            self.windows_shown.push((host.clone(), active.clone()));
            self.last_builders.insert(host.clone(), builder);
            Ok(std::mem::take(&mut self.next_window_result))
        }
        fn show_drop_preview(
            &mut self,
            _ctx: &egui::Context,
            _viewport: ViewportId,
            _builder: ViewportBuilder,
            _card: Option<&WindowCard>,
            _chrome: f32,
        ) -> bool {
            self.drop_previews_shown += 1;
            self.drop_shown_return
        }
        fn send_viewport_cmd_to(&mut self, _ctx: &egui::Context, viewport: ViewportId, cmd: ViewportCommand) {
            self.sent_commands.push((viewport, cmd));
        }
        fn request_repaint(&mut self, _ctx: &egui::Context) {
            self.root_repaints += 1;
        }
        fn request_repaint_of(&mut self, _ctx: &egui::Context, viewport: ViewportId) {
            self.repaints.push(viewport);
        }
    }

    fn monitor(name: Option<&str>, x: f32, y: f32, w: f32, h: f32, scale: f32, primary: bool) -> NativeMonitor {
        NativeMonitor { name: name.map(str::to_string), bounds: Rect::from_min_size(pos2(x, y), vec2(w, h)), scale, primary }
    }

    fn desktop_app() -> ModelApp {
        let mut app = ModelApp::new(Session::default());
        app.panels = Some(PanelHost::from_prefs(&Prefs::default(), true).unwrap());
        let monitors =
            vec![monitor(Some("Primary"), 0.0, 0.0, 1920.0, 1080.0, 1.0, true), monitor(Some("Secondary"), 1920.0, 0.0, 2560.0, 1440.0, 1.0, false)];
        app.refresh_panel_monitors(|| clone_monitors(&monitors), true);
        app
    }

    fn pop_out(host: &mut PanelHost, panel_id: &str, origin: Point) {
        pop_out_panel(host, panel_id, origin).unwrap();
    }

    #[test]
    fn canonicalize_monitors_sorts_deterministically_and_filters_invalid() {
        let raw = vec![
            monitor(Some("BadNaN"), f32::NAN, 0.0, 100.0, 100.0, 1.0, false),
            monitor(Some("BadScale"), 0.0, 0.0, 100.0, 100.0, 0.0, false),
            monitor(Some("Right"), 1920.0, 0.0, 2560.0, 1440.0, 1.0, false),
            monitor(Some("Left"), -1920.0, 0.0, 1920.0, 1080.0, 1.0, false),
            monitor(Some("Main"), 0.0, 0.0, 1920.0, 1080.0, 1.0, true),
        ];
        let out = canonicalize_monitors(raw);
        let names: Vec<_> = out.iter().map(|m| m.name.as_deref().unwrap_or("")).collect();
        assert_eq!(names, vec!["Main", "Left", "Right"]);
    }

    #[test]
    fn clamp_window_size_handles_empty_monitors_and_wgpu_limits() {
        // Empty monitor fallback never panics and guarantees min <= max.
        let (pref, (min, max)) = clamp_window_size(
            Size { width: 9000.0, height: 9000.0 },
            (Size { width: 200.0, height: 200.0 }, Size { width: 10000.0, height: 10000.0 }),
            None,
            2.0,
        );
        assert!(min.width <= max.width && min.height <= max.height);
        assert!(pref.width <= max.width && pref.height <= max.height);
        assert!((max.width - 4096.0).abs() < 1e-3);

        // Tiny monitor smaller than MIN_PANEL_SIZE keeps min/pref/max within the monitor bounds.
        let (pref, (min, max)) = clamp_window_size(
            Size { width: 500.0, height: 500.0 },
            (Size { width: 300.0, height: 300.0 }, Size { width: 600.0, height: 600.0 }),
            Some(vec2(100.0, 100.0)),
            1.0,
        );
        assert!(min.width <= max.width && min.height <= max.height);
        assert!(pref.width <= max.width && pref.height <= max.height);
        assert!(max.width <= 100.0 && max.height <= 100.0);
        assert!(min.width <= 100.0 && min.height <= 100.0);
        assert!(pref.width <= 100.0 && pref.height <= 100.0);
        assert!(is_minimized_viewport(Vec2::ZERO, Some(true)));
        assert!(is_minimized_viewport(vec2(300.0, 300.0), Some(false)));
        assert!(!is_minimized_viewport(vec2(300.0, 300.0), Some(true)));
    }

    #[test]
    fn popout_redock_preserves_live_state_and_defers_prefs_two_frames() {
        let mut app = desktop_app();
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();
        app.open_path_plot();
        app.select_definition(Some((crate::browser::Def::Material, "Steel".into())));

        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = true;
        host.shared.borrow_mut().inspector_available = true;
        host.reconcile_runtime().unwrap();
        pop_out(&mut host, PANEL_RESULTS, Point { x: 2000.0, y: 120.0 });
        host.ack_dirty();

        let ctx = egui::Context::default();
        // Frame 1: user clicks OS close (Redock).
        let mut driver = MockDriver {
            next_window_result: WindowFrameResult { action: Some(NativeWindowAction::Redock), ..Default::default() },
            ..Default::default()
        };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();

        // Visually redocked immediately to Zone::Right, preserving plot + path_editor + selection.
        let results_id = InstanceId::new(PANEL_RESULTS).unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&results_id).unwrap().location.zone(), Zone::Right);
        assert!(app.plot().is_some() && app.path_editor.is_some());
        assert_eq!(app.selected.as_ref().map(|(_, n)| n.as_str()), Some("Steel"));
        // Prefs persistence deferred by 2 root frames.
        assert!(!host.is_dirty());
        assert_eq!(host.coordinator.pending_close_persist(), Some(2));

        // Frame 2 + Frame 3 complete the deferral and mark store dirty.
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(!host.is_dirty());
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.is_dirty());
    }

    #[test]
    fn conditional_panel_clear_closes_os_window_without_redock_and_keeps_multitab_open() {
        let mut app = desktop_app();
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();
        app.select_definition(Some((crate::browser::Def::Material, "Steel".into())));

        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();

        // Pop out results alone first -> clearing results closes its OS window without redock.
        pop_out(&mut host, PANEL_RESULTS, Point { x: 2000.0, y: 120.0 });
        let ctx = egui::Context::default();
        let mut driver = MockDriver::default();
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert_eq!(driver.windows_shown.len(), 1);

        driver.windows_shown.clear();
        app.run("results.clear", json!({})).unwrap();
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(driver.windows_shown.is_empty(), "empty conditional host must not draw OS viewport");
        let results_id = InstanceId::new(PANEL_RESULTS).unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&results_id).unwrap().location.zone(), Zone::Window);

        // Join inspector as a tab inside the same Window section: clearing results keeps the window open for inspector.
        let win_sec = host.runtime.layout().layout().containing(&results_id).unwrap().id.clone();
        host.apply_mutations(
            vec![Mutation::Move {
                instance: InstanceId::new(PANEL_INSPECTOR).unwrap(),
                destination: Destination::JoinTab { section: win_sec, index: 1 },
            }],
            true,
        )
        .unwrap();
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert_eq!(driver.windows_shown.len(), 1);
        assert_eq!(driver.windows_shown[0].1.0, PANEL_INSPECTOR);
    }

    #[test]
    fn embedded_window_downgrades_to_floating_and_latches_session() {
        let mut app = desktop_app();
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();

        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();
        pop_out(&mut host, PANEL_BROWSER, Point { x: 2000.0, y: 100.0 });
        pop_out(&mut host, PANEL_RESULTS, Point { x: 2400.0, y: 120.0 });

        let ctx = egui::Context::default();
        let mut driver = MockDriver { next_window_result: WindowFrameResult { embedded_fallback: true, ..Default::default() }, ..Default::default() };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();

        assert!(!host.native_windows);
        assert!(host.native.desktop().is_none());
        assert!(!host.runtime.layout().capabilities().native_windows);
        assert!(host.coordinator.is_embedded_downgraded());
        let browser_id = InstanceId::new(PANEL_BROWSER).unwrap();
        let results_id = InstanceId::new(PANEL_RESULTS).unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&browser_id).unwrap().location.zone(), Zone::Floating);
        assert_eq!(host.runtime.layout().layout().containing(&results_id).unwrap().location.zone(), Zone::Floating);
        // Subsequent monitor updates cannot re-enable native_windows in this session.
        assert!(!host.sync_desktop_monitors(&Prefs::default(), || vec![monitor(Some("Main"), 0.0, 0.0, 1920.0, 1080.0, 1.0, true)], true));
        assert!(!host.native_windows);
    }

    #[test]
    fn exit_in_flight_cancels_pending_close_persist_and_cancel_prompt_resets() {
        let mut app = desktop_app();
        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();
        pop_out(&mut host, PANEL_BROWSER, Point { x: 2000.0, y: 100.0 });
        host.ack_dirty();

        let ctx = egui::Context::default();
        // OS child close happens first, starting 2-frame deferral.
        let mut driver = MockDriver {
            next_window_result: WindowFrameResult { action: Some(NativeWindowAction::Redock), ..Default::default() },
            ..Default::default()
        };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert_eq!(host.coordinator.pending_close_persist(), Some(2));

        // Child-before-root shutdown arrives before countdown expires -> cancels prefs write
        // and preserves the pre-close Zone::Window placement in to_prefs_value()!
        app.quit_requested = true;
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert_eq!(host.coordinator.pending_close_persist(), None);
        assert!(!host.is_dirty());
        let saved_on_exit = host.to_prefs_value().unwrap();
        let decoded_on_exit = solidcraft_panel_layout::decode(
            &serde_json::to_vec(&saved_on_exit).unwrap(),
            &crate::panels::build_registry().unwrap(),
            &crate::panels::seeded_default_layout(&crate::panels::build_registry().unwrap(), &Prefs::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(decoded_on_exit.envelope.layout.containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap().location.zone(), Zone::Window);

        // Cancelling Prompt::Unsaved restores normal OS-× redock behavior.
        app.prompt = Some(crate::files::Prompt::Unsaved(crate::files::Replace::New));
        app.answer_unsaved(crate::files::Answer::Cancel);
        assert!(!app.quit_requested);
        pop_out(&mut host, PANEL_BROWSER, Point { x: 2000.0, y: 100.0 });
        driver.next_window_result = WindowFrameResult { action: Some(NativeWindowAction::Redock), ..Default::default() };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap().location.zone(), Zone::Left);
    }

    #[test]
    fn inspector_window_title_updates_and_epoch_syncs_viewports_and_focus() {
        let mut app = desktop_app();
        let _lang = crate::i18n::scope(&app.prefs.language);
        app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210e9, "poisson": 0.3}})).unwrap();
        app.select_definition(Some((crate::browser::Def::Material, "Steel".into())));

        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();
        host.shared.borrow_mut().inspector_available = true;
        host.shared.borrow_mut().inspector_heading = crate::context::inspector_heading(&app);
        host.reconcile_runtime().unwrap();
        pop_out(&mut host, PANEL_INSPECTOR, Point { x: 2000.0, y: 120.0 });

        let ctx = egui::Context::default();
        let mut driver = MockDriver::default();
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(driver.sent_commands.iter().any(|(_, c)| matches!(c, ViewportCommand::Title(t) if t == "Inspector — Material: Steel")));

        // Establish initial epoch, then verify unchanged epoch does not spin repaints.
        host.finish_cross_viewport_sync_with(&mut driver, &app, &ctx, false, false, false);
        driver.repaints.clear();
        host.finish_cross_viewport_sync_with(&mut driver, &app, &ctx, false, false, false);
        assert!(driver.repaints.is_empty());

        // Mutating model & opening dialog triggers cross-viewport repaint and ROOT focus handoff.
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        app.open_form("model.material.create", json!({}));
        host.finish_cross_viewport_sync_with(&mut driver, &app, &ctx, false, false, false);
        assert!(!driver.repaints.is_empty());
        assert!(driver.sent_commands.iter().any(|(vp, c)| *vp == ViewportId::ROOT && matches!(c, ViewportCommand::Focus)));
    }

    #[test]
    fn desktop_native_false_refuses_popout_with_capability_and_preserves_disconnected_coords() {
        let mut host = PanelHost::from_prefs(&Prefs::default(), false).unwrap();
        let err = pop_out_panel(&mut host, PANEL_BROWSER, Point { x: 1800.0, y: 120.0 }).unwrap_err();
        assert_eq!(err.rule, Rule::Capability);

        // Saved Monitor 2 layout loads non-destructively when monitors are empty/unavailable.
        let mut desktop_host = PanelHost::from_prefs(&Prefs::default(), true).unwrap();
        pop_out(&mut desktop_host, PANEL_BROWSER, Point { x: 2400.0, y: 200.0 });
        let saved = desktop_host.to_prefs_value().unwrap();
        let prefs = Prefs { panel_layout: Some(saved), ..Prefs::default() };

        let mut app = ModelApp::new(Session::default());
        app.prefs = prefs.clone();
        app.panels = Some(PanelHost::from_prefs(&prefs, false).unwrap());
        // Empty monitor list keeps native_windows == false and preserves Monitor 2 in prefs.
        app.refresh_panel_monitors(Vec::new, true);
        let preserved = app.panels.as_ref().unwrap().to_prefs_value().unwrap();
        let decoded = solidcraft_panel_layout::decode(
            &serde_json::to_vec(&preserved).unwrap(),
            &crate::panels::build_registry().unwrap(),
            &crate::panels::seeded_default_layout(&crate::panels::build_registry().unwrap(), &Prefs::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(decoded.envelope.layout.containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap().location.zone(), Zone::Window);
        let _ = (HostCapabilities::default(), panel_kind(PANEL_BROWSER));
    }

    #[test]
    fn drop_preview_teardown_on_escape_drag_end_and_app_wide_focus_loss() {
        let mut app = desktop_app();
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();
        let ctx = egui::Context::default();
        let screen = Rect::from_min_size(pos2(0.0, 0.0), vec2(1280.0, 800.0));
        let root_raw = |events: Vec<egui::Event>| {
            let mut raw = egui::RawInput { screen_rect: Some(screen), events, ..Default::default() };
            let root = raw.viewports.entry(ViewportId::ROOT).or_default();
            root.inner_rect = Some(screen);
            root.outer_rect = Some(screen);
            root.focused = Some(true);
            root.native_pixels_per_point = Some(1.0);
            raw
        };

        // Warm up two frames so panel controls are laid out and locate Browser's draggable header rect.
        for _ in 0..2 {
            let mut out = ctx.run_ui(root_raw(Vec::new()), |ui| {
                app.ui(ui);
            });
            out.textures_delta.clear();
        }
        let header = app
            .panels
            .as_ref()
            .and_then(|h| h.host.frame())
            .and_then(|f| {
                f.layout.controls.iter().find_map(|c| {
                    matches!(&c.target, solidcraft_ui_kit::panel_runtime::HostPointerTarget::Header { instance } if instance.0 == PANEL_BROWSER)
                        .then(|| to_egui_rect(c.bounds))
                })
            })
            .expect("browser card header should be laid out");

        let drag_to_edge = |ctx: &egui::Context, app: &mut ModelApp, start: egui::Pos2, outside: egui::Pos2| {
            app.capturing = true;
            let mut out = ctx.run_ui(root_raw(vec![egui::Event::PointerMoved(start)]), |ui| app.ui(ui));
            out.textures_delta.clear();
            let mut out = ctx.run_ui(
                root_raw(vec![
                    egui::Event::PointerMoved(start),
                    egui::Event::PointerButton { pos: start, button: egui::PointerButton::Primary, pressed: true, modifiers: egui::Modifiers::NONE },
                ]),
                |ui| app.ui(ui),
            );
            out.textures_delta.clear();
            let mut out = ctx.run_ui(root_raw(vec![egui::Event::PointerMoved(outside)]), |ui| app.ui(ui));
            out.textures_delta.clear();
            app.capturing = false;
        };

        // 1. Drag out into monitor space and verify preview paints via MockDriver, then Esc tears down preview.
        drag_to_edge(&ctx, &mut app, header.center(), pos2(-180.0, header.center().y));
        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();
        assert!(host.host.window_landing_preview().is_some(), "drag past window edge must activate drop preview");
        let mut driver = MockDriver { drop_shown_return: true, root_focused: Some(true), ..Default::default() };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(driver.drop_previews_shown > 0, "active drag-out must show drop preview");

        // Escape cancels the drag gesture and stops showing drop preview.
        driver.escape = true;
        driver.drop_previews_shown = 0;
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.host.window_landing_preview().is_none());
        assert_eq!(driver.drop_previews_shown, 0);

        // 2. Releasing pointer (drag end) also clears window_landing_preview.
        app.panels = Some(host);
        drag_to_edge(&ctx, &mut app, header.center(), pos2(-180.0, header.center().y));
        app.capturing = true;
        let mut out = ctx.run_ui(
            root_raw(vec![egui::Event::PointerButton {
                pos: pos2(-180.0, header.center().y),
                button: egui::PointerButton::Primary,
                pressed: false,
                modifiers: egui::Modifiers::NONE,
            }]),
            |ui| app.ui(ui),
        );
        out.textures_delta.clear();
        app.capturing = false;
        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();
        let mut driver = MockDriver::default();
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert_eq!(driver.drop_previews_shown, 0);

        // 3. Root losing focus while a popped-out child window IS focused does NOT tear down preview;
        //    root losing focus with no child focused (app-wide focus loss) DOES tear down preview.
        pop_out(&mut host, PANEL_RESULTS, Point { x: 2000.0, y: 120.0 });
        host.shared.borrow_mut().results_available = true;
        host.reconcile_runtime().unwrap();
        driver.next_window_result = WindowFrameResult { focused: true, ..Default::default() };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.coordinator.any_child_focused);
        assert_eq!(host.coordinator.active_viewports().len(), 1);

        // With child focused, root_focused = Some(false) must keep gesture intact.
        app.panels = Some(host);
        drag_to_edge(&ctx, &mut app, header.center(), pos2(-180.0, header.center().y));
        let mut host = app.panels.take().unwrap();
        host.shared.borrow_mut().results_available = crate::results::shown(&app.session).is_some();
        host.shared.borrow_mut().inspector_available = crate::context::inspector_active(&app);
        host.reconcile_runtime().unwrap();
        assert!(host.host.window_landing_preview().is_some(), "drag past window edge must activate drop preview");
        driver.root_focused = Some(false);
        driver.drop_previews_shown = 0;
        driver.next_window_result = WindowFrameResult { focused: true, ..Default::default() };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.coordinator.any_child_focused);
        assert!(host.host.window_landing_preview().is_some());

        // Now lose child focus (app-wide focus loss) -> next frame cancels gesture.
        driver.next_window_result = WindowFrameResult { focused: false, ..Default::default() };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(!host.coordinator.any_child_focused);
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.host.window_landing_preview().is_none());
    }

    #[test]
    fn saved_native_origin_seeds_builder_position_and_falls_back_when_monitor_missing() {
        let mut saved_app = desktop_app();
        {
            let mut host = saved_app.panels.take().unwrap();
            pop_out(&mut host, PANEL_BROWSER, Point { x: 2000.0, y: 120.0 });
            saved_app.panels = Some(host);
        }
        let saved_layout = saved_app.panels.as_ref().unwrap().to_prefs_value().unwrap();
        let prefs = Prefs { panel_layout: Some(saved_layout), ..Prefs::default() };

        // 1. Fresh host with Secondary monitor connected places Browser on Monitor 2.
        let mut restored_app = ModelApp::new(Session::default());
        restored_app.prefs = prefs.clone();
        restored_app.panels = Some(PanelHost::from_prefs(&prefs, true).unwrap());
        let two_monitors =
            vec![monitor(Some("Primary"), 0.0, 0.0, 1920.0, 1080.0, 1.0, true), monitor(Some("Secondary"), 1920.0, 0.0, 2560.0, 1440.0, 1.0, false)];
        restored_app.refresh_panel_monitors(|| clone_monitors(&two_monitors), true);

        let ctx = egui::Context::default();
        let mut driver = MockDriver::default();
        let mut host = restored_app.panels.take().unwrap();
        host.show_desktop_windows_with(&mut driver, &mut restored_app, &ctx).unwrap();
        let browser_host = match &host.runtime.layout().layout().containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap().location {
            SectionLocation::Window { host, .. } => host.clone(),
            other => panic!("expected Window location, got {other:?}"),
        };
        let pos = driver.last_builders.get(&browser_host).and_then(|b| b.position).expect("seeded native builder must include position");
        assert!((pos.x - 2000.0).abs() < 1e-3 && (pos.y - 120.0).abs() < 1e-3, "got {pos:?}");

        // 2. When Secondary monitor is unavailable, builder falls back inside Primary monitor.
        let mut single_app = ModelApp::new(Session::default());
        single_app.prefs = prefs.clone();
        single_app.panels = Some(PanelHost::from_prefs(&prefs, true).unwrap());
        let primary_only = vec![monitor(Some("Primary"), 0.0, 0.0, 1920.0, 1080.0, 1.0, true)];
        single_app.refresh_panel_monitors(|| clone_monitors(&primary_only), true);
        let mut driver_single = MockDriver::default();
        let mut single_host = single_app.panels.take().unwrap();
        single_host.show_desktop_windows_with(&mut driver_single, &mut single_app, &ctx).unwrap();
        let single_pos =
            driver_single.last_builders.get(&browser_host).and_then(|b| b.position).expect("fallback native builder must include position");
        assert!((0.0..1920.0).contains(&single_pos.x), "expected fallback on primary monitor, got {single_pos:?}");
    }

    #[test]
    fn observed_native_window_geometry_updates_layout_and_marks_dirty() {
        let mut app = desktop_app();
        let mut host = app.panels.take().unwrap();
        pop_out(&mut host, PANEL_BROWSER, Point { x: 2000.0, y: 120.0 });
        host.ack_dirty();

        let ctx = egui::Context::default();
        // Minimized or unmeasured frames must not alter layout or mark dirty.
        let mut driver = MockDriver {
            next_window_result: WindowFrameResult { minimized_or_hidden: true, observed_geometry: None, ..Default::default() },
            ..Default::default()
        };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(!host.is_dirty());

        // Normal window move + resize updates layout + instance memory and marks store dirty.
        driver.next_window_result = WindowFrameResult {
            observed_geometry: Some((Some(Point { x: 2150.0, y: 240.0 }), Size { width: 520.0, height: 610.0 })),
            ..Default::default()
        };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.is_dirty());

        let saved = host.to_prefs_value().unwrap();
        let reg = crate::panels::build_registry().unwrap();
        let fallback = crate::panels::seeded_default_layout(&reg, &Prefs::default()).unwrap();
        let decoded = solidcraft_panel_layout::decode(saved.as_bytes(), &reg, &fallback).unwrap();
        let browser_sec = decoded.envelope.layout.containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap();
        assert!((browser_sec.size.width - 520.0).abs() < 1e-3);
        assert!((browser_sec.size.height - 610.0).abs() < 1e-3);
        match &browser_sec.location {
            SectionLocation::Window { origin, .. } => {
                assert!((origin.x - 2150.0).abs() < 1e-3 && (origin.y - 240.0).abs() < 1e-3);
            }
            other => panic!("expected Window location, got {other:?}"),
        }
    }

    #[test]
    fn undo_shortcut_is_ignored_while_prompt_is_open() {
        let mut app = desktop_app();
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        assert_eq!(app.session.model.nodes().len(), 1);

        let ctx = egui::Context::default();
        let screen = Rect::from_min_size(pos2(0.0, 0.0), vec2(1280.0, 800.0));
        app.prompt = Some(crate::files::Prompt::Unsaved(crate::files::Replace::New));
        let mut out = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(screen),
                events: vec![egui::Event::Key {
                    key: egui::Key::Z,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::COMMAND,
                }],
                ..Default::default()
            },
            |ui| app.shortcuts(ui.ctx()),
        );
        out.textures_delta.clear();
        assert_eq!(app.session.model.nodes().len(), 1, "Ctrl+Z must not undo while prompt is open");

        app.prompt = None;
        let mut out = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(screen),
                events: vec![egui::Event::Key {
                    key: egui::Key::Z,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::COMMAND,
                }],
                ..Default::default()
            },
            |ui| app.shortcuts(ui.ctx()),
        );
        out.textures_delta.clear();
        assert_eq!(app.session.model.nodes().len(), 0, "Ctrl+Z must undo once prompt is closed");
    }

    #[test]
    fn completed_native_close_redock_survives_capability_drop_and_return_during_deferral() {
        let mut app = desktop_app();
        let mut host = app.panels.take().unwrap();
        pop_out(&mut host, PANEL_BROWSER, Point { x: 2000.0, y: 120.0 });
        host.ack_dirty();

        let ctx = egui::Context::default();
        let mut driver = MockDriver {
            next_window_result: WindowFrameResult { action: Some(NativeWindowAction::Redock), ..Default::default() },
            ..Default::default()
        };
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        let browser_id = InstanceId::new(PANEL_BROWSER).unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&browser_id).unwrap().location.zone(), Zone::Left);
        assert_eq!(host.coordinator.pending_close_persist(), Some(2));

        // Native support drops then returns before the 2-frame countdown expires.
        assert!(!host.sync_desktop_monitors(&Prefs::default(), Vec::new, false));
        assert!(host.sync_desktop_monitors(&Prefs::default(), || vec![monitor(Some("Main"), 0.0, 0.0, 1920.0, 1080.0, 1.0, true)], true));
        assert_eq!(host.runtime.layout().layout().containing(&browser_id).unwrap().location.zone(), Zone::Left);

        // Countdown finishes cleanly and marks the store dirty.
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        host.show_desktop_windows_with(&mut driver, &mut app, &ctx).unwrap();
        assert!(host.is_dirty());
    }
}
