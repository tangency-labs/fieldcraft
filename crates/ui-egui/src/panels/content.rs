//! Panel content dispatch, availability/presentation controllers, and deferred commands.

use std::cell::RefCell;
use std::collections::VecDeque;
use std::rc::Rc;

use serde_json::Value;
use solidcraft_panel_layout::{Availability, ControlId, InstanceId, PanelResult, Point, Rule, error, valid_id};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::panel_runtime::{
    ContentCommandSink, ContentOwner, ContentServices, ContentSlot, FramePanelContent, LifecycleEvent, LocalRequest, LocalScope, PanelContext,
    PanelController, PanelLocalServices, PanelModels, PanelPresentation, RuntimeAction,
};
use solidcraft_ui_text::tr;

use super::{PANEL_BROWSER, PANEL_HISTORY, PANEL_INSPECTOR, PANEL_RESULTS};
use crate::ModelApp;
use crate::theme::Tokens;

const MAX_DEFERRED_COMMANDS: usize = 64;

/// Deferred UI/model commands emitted while panel controllers or bodies render.
/// Executed in Step 7 of `ModelApp::shell()` after `WorkspacePanelHost::finish`.
#[derive(Clone)]
pub(crate) enum PanelCommand {
    SetShowTree(bool),
    SetShowTimeline(bool),
    ClearResults,
    ClearInspectorTop,
    HistoryTransport(kit::HistoryTransport),
}

impl std::fmt::Debug for PanelCommand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SetShowTree(v) => f.debug_tuple("SetShowTree").field(v).finish(),
            Self::SetShowTimeline(v) => f.debug_tuple("SetShowTimeline").field(v).finish(),
            Self::ClearResults => f.write_str("ClearResults"),
            Self::ClearInspectorTop => f.write_str("ClearInspectorTop"),
            Self::HistoryTransport(_) => f.write_str("HistoryTransport(..)"),
        }
    }
}

/// Frame-synchronized snapshot read by `FieldCraftPanelController`s.
#[derive(Debug, Default)]
pub(crate) struct PanelSharedState {
    pub(crate) history_available: bool,
    pub(crate) results_available: bool,
    pub(crate) inspector_available: bool,
    pub(crate) inspector_heading: Option<(String, String)>,
    pub(crate) commands: VecDeque<PanelCommand>,
}

impl PanelSharedState {
    pub(crate) fn push_command(&mut self, command: PanelCommand) {
        if self.commands.len() < MAX_DEFERRED_COMMANDS {
            self.commands.push_back(command);
        } else {
            log::warn!("dropping panel command because queue reached limit ({MAX_DEFERRED_COMMANDS})");
        }
    }

    pub(crate) fn drain_commands(&mut self) -> Vec<PanelCommand> {
        self.commands.drain(..).collect()
    }
}

pub(crate) struct FieldCraftPanelController {
    id: &'static str,
    shared: Rc<RefCell<PanelSharedState>>,
}

impl FieldCraftPanelController {
    pub(crate) fn new(id: &'static str, shared: Rc<RefCell<PanelSharedState>>) -> Self {
        Self { id, shared }
    }
}

impl PanelController for FieldCraftPanelController {
    fn availability(&self, _owner: &ContentOwner, _models: &dyn PanelModels) -> Availability {
        let shared = self.shared.borrow();
        let available = match self.id {
            PANEL_BROWSER => true,
            PANEL_HISTORY => shared.history_available,
            PANEL_RESULTS => shared.results_available,
            PANEL_INSPECTOR => shared.inspector_available,
            _ => false,
        };
        if available { Availability::Available } else { Availability::Unavailable { reason: tr("ui.panel_host.rule.availability").into() } }
    }

    fn update(&mut self, _owner: &ContentOwner, _services: &mut ContentServices<'_>) -> PanelResult<()> {
        Ok(())
    }

    fn presentation(&self) -> PanelPresentation {
        match self.id {
            PANEL_BROWSER => {
                PanelPresentation { title: tr("panel.browser").into(), icon: "list".into(), status: None, attention: false, actions: Vec::new() }
            }
            PANEL_HISTORY => {
                let actions = [
                    ("first", tr("ui.history_first")),
                    ("back", tr("ui.history_back")),
                    ("forward", tr("ui.history_forward")),
                    ("last", tr("ui.history_last")),
                ]
                .into_iter()
                .filter_map(|(id, label)| {
                    Some(RuntimeAction { id: ControlId::new(id).ok()?, label: label.into(), enabled: true, disabled_reason: None, pressed: None })
                })
                .collect();
                PanelPresentation { title: tr("panel.history").into(), icon: "ui.views".into(), status: None, attention: false, actions }
            }
            PANEL_RESULTS => PanelPresentation {
                title: tr("panel.results").into(),
                icon: "results.deformed".into(),
                status: None,
                attention: false,
                actions: Vec::new(),
            },
            PANEL_INSPECTOR => {
                let shared = self.shared.borrow();
                let (icon, title) = shared.inspector_heading.clone().unwrap_or_else(|| ("list".into(), tr("panel.inspector").into()));
                PanelPresentation { title, icon, status: None, attention: false, actions: Vec::new() }
            }
            _ => PanelPresentation { title: self.id.into(), icon: "list".into(), status: None, attention: false, actions: Vec::new() },
        }
    }

    fn prepare_lifecycle(&self, _event: &LifecycleEvent) -> PanelResult<()> {
        Ok(())
    }

    fn lifecycle(&mut self, event: &LifecycleEvent) {
        let cmd = match (event, self.id) {
            (LifecycleEvent::Close, PANEL_BROWSER) => Some(PanelCommand::SetShowTree(false)),
            (LifecycleEvent::Close, PANEL_HISTORY) => Some(PanelCommand::SetShowTimeline(false)),
            (LifecycleEvent::Close, PANEL_RESULTS) => Some(PanelCommand::ClearResults),
            (LifecycleEvent::Close, PANEL_INSPECTOR) => Some(PanelCommand::ClearInspectorTop),
            (LifecycleEvent::Reopen, PANEL_BROWSER) => Some(PanelCommand::SetShowTree(true)),
            (LifecycleEvent::Reopen, PANEL_HISTORY) => Some(PanelCommand::SetShowTimeline(true)),
            _ => None,
        };
        if let Some(cmd) = cmd {
            self.shared.borrow_mut().push_command(cmd);
        }
    }

    fn render(&mut self, _slot: ContentSlot, _ui: &mut egui::Ui, _context: &PanelContext, _services: &mut ContentServices<'_>) -> PanelResult<()> {
        Err(error(Rule::Lifecycle, "FieldCraft panel bodies render through FieldCraftFrameContent"))
    }

    fn action(&mut self, _owner: &ContentOwner, action: &ControlId, _services: &mut ContentServices<'_>) -> PanelResult<()> {
        if self.id == PANEL_HISTORY {
            let transport = match action.0.as_str() {
                "first" => kit::HistoryTransport::First,
                "back" => kit::HistoryTransport::Back,
                "forward" => kit::HistoryTransport::Forward,
                "last" => kit::HistoryTransport::Last,
                _ => return Err(error(Rule::Capability, "unknown history header action")),
            };
            self.shared.borrow_mut().push_command(PanelCommand::HistoryTransport(transport));
        }
        Ok(())
    }
}

/// Frame-scoped renderer borrowing `ModelApp` cleanly during `WorkspacePanelHost::show_deferred_with_content`.
pub(crate) struct FieldCraftFrameContent<'a> {
    pub(crate) app: &'a mut ModelApp,
}

impl FramePanelContent for FieldCraftFrameContent<'_> {
    fn render(
        &mut self,
        slot: ContentSlot,
        ui: &mut egui::Ui,
        context: &PanelContext,
        _content: &mut dyn PanelController,
        _services: &mut ContentServices<'_>,
    ) -> PanelResult<()> {
        if !matches!(slot, ContentSlot::Body | ContentSlot::Inline) {
            return Ok(());
        }
        let outer = egui::Rect::from_min_size(
            egui::pos2(context.outer.origin.x, context.outer.origin.y),
            egui::vec2(context.outer.size.width, context.outer.size.height),
        );
        let t = Tokens::get(self.app.ui.dark);
        match context.instance.0.as_str() {
            PANEL_BROWSER => {
                if !context.preview_geometry {
                    self.app.place("browser", outer);
                }
                ui.push_id("fc_panel_browser", |ui| {
                    crate::browser::body(self.app, ui);
                });
            }
            PANEL_HISTORY => {
                if !context.preview_geometry {
                    self.app.place("timeline", outer);
                }
                let floating = slot != ContentSlot::Inline;
                ui.push_id("fc_panel_history", |ui| {
                    crate::timeline::body(self.app, ui, floating);
                });
            }
            PANEL_RESULTS => {
                if !context.preview_geometry {
                    self.app.place("context:results", outer);
                }
                ui.push_id("fc_panel_results", |ui| {
                    crate::results_panel::body(self.app, ui, &t);
                });
            }
            PANEL_INSPECTOR => {
                if !context.preview_geometry
                    && let Some(key) = crate::context::inspector_region_key(self.app)
                {
                    self.app.place(key, outer);
                }
                ui.push_id("fc_panel_inspector", |ui| {
                    crate::context::inspector_body(self.app, ui, &t);
                });
            }
            _ => {}
        }
        Ok(())
    }
}

/// Format an OS window title for the popped-out Inspector panel (used in PR 4b).
#[allow(dead_code)]
pub(crate) fn inspector_window_title(app: &ModelApp) -> String {
    match crate::context::inspector_heading(app) {
        Some((_, title)) => format!("{} — {title}", tr("panel.inspector")),
        None => tr("panel.inspector").into(),
    }
}

/// Execute deferred commands drained from `PanelSharedState` after the frame's panel pass.
pub(crate) fn execute_commands(app: &mut ModelApp, commands: Vec<PanelCommand>) {
    for cmd in commands {
        match cmd {
            PanelCommand::SetShowTree(open) => {
                app.ui.show_tree = open;
            }
            PanelCommand::SetShowTimeline(open) => {
                app.ui.show_timeline = open;
            }
            PanelCommand::ClearResults => {
                let _ = app.run("results.clear", serde_json::json!({}));
            }
            PanelCommand::ClearInspectorTop => {
                crate::context::clear_inspector_top(app);
            }
            PanelCommand::HistoryTransport(control) => {
                crate::timeline::run_transport(app, control);
            }
        }
    }
}

pub(crate) struct FieldCraftModels;
impl PanelModels for FieldCraftModels {
    fn read(&self, _owner: &ContentOwner, _key: &str) -> PanelResult<Value> {
        Err(error(Rule::MissingTarget, "FieldCraft panel bodies read ModelApp directly"))
    }
}

pub(crate) struct FieldCraftCommands;
impl ContentCommandSink for FieldCraftCommands {
    fn execute(&mut self, _owner: &ContentOwner, _command: &str, _input: Value) -> PanelResult<Value> {
        Err(error(Rule::Capability, "FieldCraft panel commands use PanelCommand queue"))
    }
}

#[derive(Default)]
pub(crate) struct FieldCraftLocal {
    keyboard: Option<(LocalScope, InstanceId)>,
}

impl PanelLocalServices for FieldCraftLocal {
    fn request(&mut self, scope: &LocalScope, request: LocalRequest) -> PanelResult<()> {
        scope.validate()?;
        request.validate()?;
        valid_id(&request.instance().0)?;
        match request {
            LocalRequest::KeyboardOwner { instance, .. } => {
                self.keyboard = Some((scope.clone(), instance));
            }
            LocalRequest::ReleaseKeyboard { instance } => {
                if self.keyboard.as_ref().is_some_and(|(_, owner)| *owner == instance) {
                    self.keyboard = None;
                }
            }
            LocalRequest::Overlay { region, .. } => {
                let _ = Point { x: region.bounds.origin.x, y: region.bounds.origin.y };
            }
            _ => {}
        }
        Ok(())
    }

    fn release(&mut self, scope: &LocalScope, instance: &InstanceId) {
        if self.keyboard.as_ref().is_some_and(|(s, id)| s == scope && id == instance) {
            self.keyboard = None;
        }
    }
}
