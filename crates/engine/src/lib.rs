//! The command engine.
//!
//! Everything a user can do is a command ([`CommandSpec`]) with a JSON form. The desktop UI,
//! scripts, the CLI and the control channel all run the same commands through
//! [`Session::execute`], which never lets a failure or a panic damage the session: the model and
//! the geometry are restored and an [`EngineError`] comes back.
//!
//! Commands that change the model are journalled and can be undone and redone
//! (`edit.undo`, `edit.redo`; [`history`]). The session is saved to and opened from `.fieldcraft`
//! files (`file.save`, `file.open`; [`fieldcraft_io`]) and autosaved for recovery after a crash
//! ([`autosave`]).
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod autosave;
pub mod bug_report;
pub mod catalogue;
pub mod cmd;
pub mod crash;
pub mod flow;
pub mod logfile;
pub mod machine;

pub use cmd::analytic_surfaces;
pub use cmd::replace::{MESHERS, Replacement, replacement, thousands};
pub use cmd::update::enable_checks as enable_update_checks;
/// Write an output (a picture, an export) where a user named, never over a model or design.
pub use fieldcraft_io::fs::write_output;
pub use fieldcraft_io::{Item as ViewItem, ItemKind};
pub mod history;
pub mod progress;
pub mod requirement;
pub mod results;
pub mod sample;
pub mod schema;
pub mod units;
pub mod usage;
pub mod view;

use std::collections::BTreeMap;
use std::sync::Arc;

use fieldcraft_explicit::ExplicitResult;
use fieldcraft_geometry::Design;
use fieldcraft_geometry::tagged_mesh::TaggedSurface;
use fieldcraft_model::Model;
use serde_json::Value;

use crate::history::History;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum EngineError {
    #[error("unknown command `{0}`")]
    UnknownCommand(String),
    #[error("{0}")]
    Params(String),
    /// One parameter is missing or wrong: the schema's [`schema::Problem`], at the field's JSON
    /// path (`elastic/young`). The message names the key ([`schema::Problem::for_agents`]); a
    /// form shows it beside the field, by its label ([`schema::Problem::for_form`]).
    #[error("{}", .0.for_agents())]
    Param(schema::Problem),
    #[error(transparent)]
    Model(#[from] fieldcraft_model::Error),
    #[error(transparent)]
    Geometry(#[from] fieldcraft_geometry::Error),
    #[error(transparent)]
    Mesh(#[from] fieldcraft_mesh::Error),
    #[error(transparent)]
    File(#[from] fieldcraft_io::Error),
    #[error(transparent)]
    Deck(#[from] fieldcraft_decks::Error),
    #[error(transparent)]
    Solver(#[from] fieldcraft_solver::Error),
    #[error(transparent)]
    Explicit(#[from] fieldcraft_explicit::Error),
    #[error(transparent)]
    Fluid(#[from] fieldcraft_lbm::Error),
    /// The session lacks what the command needs ([`requirement`]); nothing was done.
    #[error("{0}")]
    Unavailable(requirement::Unavailable),
    /// A command panicked; the model is unchanged.
    #[error("internal error in {command}: {message}")]
    Internal { command: String, message: String },
    /// An error of what `context` names (`body `Block``): kept whole, so its numbers are shown
    /// in the display units at the command's boundary ([`EngineError::shown`]).
    #[error("{context}: {source}")]
    Context { context: String, source: Box<EngineError> },
}

impl EngineError {
    /// Whether the command stopped because it was cancelled
    /// ([`fieldcraft_model::progress::Progress::cancel`]).
    pub fn is_cancelled(&self) -> bool {
        use fieldcraft_model::Error::Cancelled;
        if let EngineError::Context { source, .. } = self {
            return source.is_cancelled();
        }
        matches!(
            self,
            EngineError::Model(Cancelled)
                | EngineError::Deck(fieldcraft_decks::Error::Model(Cancelled))
                | EngineError::Geometry(fieldcraft_geometry::Error::Model(Cancelled))
                | EngineError::Mesh(fieldcraft_mesh::Error::Model(Cancelled))
                | EngineError::Solver(fieldcraft_solver::Error::Model(Cancelled))
                | EngineError::Explicit(fieldcraft_explicit::Error::Model(Cancelled))
                | EngineError::Fluid(fieldcraft_lbm::Error::Model(Cancelled) | fieldcraft_lbm::Error::Mesh(fieldcraft_mesh::Error::Model(Cancelled)))
        )
    }

    /// A body's mesher failure naming the body, or the cancellation as it is.
    pub(crate) fn of_body(body: &str, e: fieldcraft_mesh::Error) -> Self {
        match e {
            fieldcraft_mesh::Error::Model(fieldcraft_model::Error::Cancelled) => EngineError::Model(fieldcraft_model::Error::Cancelled),
            e => EngineError::context(format!("body `{body}`"), e),
        }
    }

    /// `e` as an error of what `context` names; a cancellation stays one.
    pub(crate) fn context(context: impl Into<String>, e: impl Into<EngineError>) -> Self {
        let e = e.into();
        if e.is_cancelled() { e } else { EngineError::Context { context: context.into(), source: Box::new(e) } }
    }

    /// The error with the numbers its meshers' messages name, which are in `from` (the command's
    /// frame), shown in `to` (the display system) with their units.
    pub fn shown(self, from: &fieldcraft_model::DisplaySystem, to: &fieldcraft_model::DisplaySystem) -> Self {
        let units = units::MeshUnits { from, to };
        match self {
            EngineError::Mesh(e) => EngineError::Mesh(e.shown(&units)),
            EngineError::Fluid(fieldcraft_lbm::Error::Mesh(e)) => EngineError::Fluid(fieldcraft_lbm::Error::Mesh(e.shown(&units))),
            EngineError::Context { context, source } => EngineError::Context { context, source: Box::new(source.shown(from, to)) },
            e => e,
        }
    }

    /// A problem with the parameter at `path` (`elastic/young`).
    pub fn param(path: impl Into<String>, kind: schema::ProblemKind) -> Self {
        Self::Param(schema::Problem::at(path, kind))
    }

    /// A parameter problem with something wrong that the schema does not describe: what the
    /// value must be or do instead.
    pub fn invalid(path: impl Into<String>, why: impl Into<String>) -> Self {
        Self::param(path, schema::ProblemKind::Invalid(why.into()))
    }
}

pub type Result<T> = std::result::Result<T, EngineError>;

/// The result of a solved step, from the implicit, the explicit or the lattice Boltzmann solver.
#[derive(Debug, Clone, PartialEq)]
pub enum StepResult {
    Implicit(fieldcraft_solver::StepResult),
    Explicit(ExplicitResult),
    Fluid(fieldcraft_lbm::FluidResult),
}

impl StepResult {
    pub fn step(&self) -> &str {
        match self {
            StepResult::Implicit(r) => r.step(),
            StepResult::Explicit(r) => &r.step,
            StepResult::Fluid(r) => &r.step,
        }
    }

    /// The node labels the result's per-node values follow, sorted.
    pub fn nodes(&self) -> &[u32] {
        match self {
            StepResult::Implicit(r) => r.nodes(),
            StepResult::Explicit(r) => &r.nodes,
            // A fluid result is per cell, not per node.
            StepResult::Fluid(_) => &[],
        }
    }

    /// The result as a file keeps it.
    pub fn as_stored(&self) -> fieldcraft_io::results::SolvedRef<'_> {
        use fieldcraft_io::results::SolvedRef;
        match self {
            StepResult::Implicit(r) => SolvedRef::Implicit(r),
            StepResult::Explicit(r) => SolvedRef::Explicit(r),
            StepResult::Fluid(r) => SolvedRef::Fluid(r),
        }
    }
}

impl From<fieldcraft_io::results::Solved> for StepResult {
    fn from(r: fieldcraft_io::results::Solved) -> Self {
        use fieldcraft_io::results::Solved;
        match r {
            Solved::Implicit(r) => StepResult::Implicit(r),
            Solved::Explicit(r) => StepResult::Explicit(r),
            Solved::Fluid(r) => StepResult::Fluid(r),
        }
    }
}

/// What a command works on, which decides how [`Session::execute`] records it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CommandKind {
    /// Reads or changes the model and its geometry. A run that changes something is journalled
    /// and becomes an undo step.
    Model,
    /// Works on the session itself: opens, saves, undoes and redoes. Never journalled; such a
    /// command keeps the journal and the undo history in step itself.
    Session,
}

/// A command: one user-visible action.
pub struct CommandSpec {
    /// Our own dotted, lower-case id: `<area>.<thing>[.<variant>]`, such as `mesh.node.add`.
    pub id: &'static str,
    pub label: &'static str,
    /// The JSON parameters, documented for scripts, `fieldcraft-cli commands` and agents.
    pub params: &'static str,
    /// The parameters field by field, for forms and JSON Schema; `None` when not described yet
    /// (the form is then the JSON text alone).
    pub form: Option<&'static [schema::Field]>,
    pub kind: CommandKind,
    /// What the session must hold before the command can run ([`requirement`]): checked by
    /// [`Session::execute`], and by the UI to grey the command with the missing pieces named.
    pub requires: &'static [requirement::Requirement],
    pub run: fn(&mut Session, &Value) -> Result<Value>,
}

/// All commands.
pub fn commands() -> &'static [CommandSpec] {
    cmd::ALL
}

impl CommandSpec {
    /// The command as JSON for listings (`fieldcraft-cli commands --json`, the control channel's
    /// `engine.commands`, MCP `list_commands`): id, label, the parameter doc and, when described,
    /// the parameters as JSON Schema, and what it requires.
    pub fn listing(&self) -> Value {
        let mut v = serde_json::json!({"id": self.id, "label": self.label, "params": self.params, "requires": self.requires});
        if let Some(form) = self.form {
            v["schema"] = schema::json_schema(form);
        }
        v
    }
}

pub fn find_command(id: &str) -> Option<&'static CommandSpec> {
    commands().iter().find(|c| c.id == id)
}

/// Refuse a parameter the command does not take, naming the ones it does: a script that sets
/// `deformed` on `results.contour` learns it belongs elsewhere instead of being ignored.
fn check_keys(spec: &CommandSpec, params: &Value) -> Result<()> {
    match spec.form.and_then(|form| schema::unknown_key(form, params)) {
        Some((key, accepted)) if accepted.is_empty() => Err(EngineError::Params(format!("{} takes no parameters (got `{key}`)", spec.id))),
        Some((key, accepted)) => Err(EngineError::Params(format!("{} has no parameter `{key}`; it takes {}", spec.id, accepted.join(", ")))),
        None => Ok(()),
    }
}

/// A new revision number, unique in the process: every state of every session has its own.
pub(crate) fn next_revision() -> u64 {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
    NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

/// The defaults of the work caps that commands take as parameters: what a command may spend
/// when its parameters name no cap of their own. A parameter still overrides them within the
/// command's documented range. The defaults are the production ones; a host may lower them (the
/// tests do, so hostile input is refused quickly).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Limits {
    /// `mesh.tet.generate`'s `max_steiner`.
    pub max_steiner: usize,
    /// `mesh.tet.generate`'s `max_tets`.
    pub max_tets: usize,
    /// `mesh.hex.octree`'s `max_elements`.
    pub max_hex_elements: usize,
    /// Threads the solvers use (the `solve.*` commands' `threads`): by default
    /// `FIELDCRAFT_THREADS` or the available cores, 1 to [`fieldcraft_parallel::MAX_THREADS`].
    /// Results do not depend on it.
    pub threads: usize,
    /// The memory in bytes a fluid solve may take once its caps (`max_cells`,
    /// `max_cell_updates`, `max_frame_bytes`) are raised above their defaults, which also bounds
    /// how far they may be raised: half the machine's physical memory ([`machine`]), 0 where it
    /// is unknown.
    pub memory: u64,
}

impl Default for Limits {
    fn default() -> Self {
        let cdt = fieldcraft_mesh::cdt::CdtOptions::default();
        Limits {
            max_steiner: cdt.max_steiner,
            max_tets: cdt.max_tets,
            max_hex_elements: cmd::octree::DEFAULT_MAX_ELEMENTS,
            threads: fieldcraft_parallel::default_threads(),
            memory: machine::memory_bound(),
        }
    }
}

/// The model being edited, the geometry it is built on and the commands that built it.
#[derive(Debug, Clone, Default)]
pub struct Session {
    pub model: Model,
    /// The open SolidCraft design (`geometry.open`); shared, so a snapshot costs nothing.
    pub design: Option<Arc<Design>>,
    /// Face-tagged boundary surfaces (`geometry.import.mesh`): the input to volume meshing.
    /// Each is shared, so the per-command snapshot is cheap.
    pub surfaces: Vec<Arc<TaggedSurface>>,
    /// The files the design and the surfaces were read from, with the hash of what was read:
    /// what a `.fieldcraft` file records to open them again.
    pub sources: fieldcraft_io::Geometry,
    /// Commands that changed the model, in order (a script that rebuilds it). Undo takes the
    /// last one off; redo puts it back. Its parameters are as they were given: bare numbers in
    /// `journal_header.units`, others with their units.
    pub journal: Vec<(String, Value)>,
    /// The unit system the journal's bare numbers are in (the input system of its first
    /// command), and, for a journal of a model that predates SI storage, the system the model was
    /// shown in before its first command.
    pub journal_header: fieldcraft_io::JournalHeader,
    /// The display system of a model that sets none of its own ([`Session::units`]): the user's
    /// setting in the app, mm-t-s-N-MPa elsewhere. A file replacing the session keeps it.
    pub default_units: fieldcraft_model::DisplaySystem,
    /// While a script with a units header replays: the system its bare numbers are in.
    pinned_input: Option<fieldcraft_model::DisplaySystem>,
    /// While a journal from before format 3 replays (its header names `model_units`).
    legacy_replay: bool,
    /// The `.fieldcraft` file the session was opened from or last saved to.
    pub path: Option<String>,
    /// The format of the file at `path` when it is older than the one this version writes: the model
    /// was migrated from an older file (`file.open`), or recovered from an autosave of one
    /// (`file.recover`). The first save over that file keeps a copy of it first
    /// (`<stem>.v<N>.fieldcraft`, [`fieldcraft_io::fs::save_model`]); then it is `None`.
    pub migrated_from: Option<u32>,
    /// What migrating that older file changed in what it means ([`fieldcraft_io::ModelFile::notes`]):
    /// shown with the migrated-from note until the first save.
    pub migration_notes: Vec<String>,
    /// The sample the session was opened as (`file.sample`), by its label; the document is named
    /// after it until it is saved under a name of its own.
    pub sample: Option<&'static str>,
    history: History,
    /// Identifies the current state: it changes with every edit, undo and redo, and comes back
    /// when an undo returns to an earlier state.
    revision: u64,
    /// The revision last saved to (or opened from) `path`; `None` when the content was never
    /// saved, where only the empty start (revision 0) is clean.
    saved_revision: Option<u64>,
    /// The latest result of each solved step, by step name, with the revision it was solved at;
    /// [`Session::result`] hands out only those that describe the current state.
    results: BTreeMap<String, (u64, Arc<StepResult>)>,
    /// The definition each solved step was solved from ([`fieldcraft_io::results::definition_hash`]): an edit
    /// that leaves it the same keeps the step's result current. Files keep it with the result.
    result_inputs: BTreeMap<String, [u8; 32]>,
    /// Counts the results stored (a solve, or a file opened): a result solved since the last save
    /// is unsaved too, since files keep results.
    results_stored: u64,
    /// [`Session::results_stored`] as last saved (or opened).
    saved_results: u64,
    /// What the results display shows (`results.*`); `None` shows the model.
    pub results_view: Option<results::ResultsView>,
    /// The default work caps of commands; a file replacing the session keeps them.
    pub limits: Limits,
    /// The folders commands may read and write files in (an automation client's roots,
    /// docs/mcp.md); unconfined by default. A file replacing the session keeps them.
    pub file_access: fieldcraft_io::roots::Roots,
    /// Identifies the open document: it changes when another model, sample or design is opened
    /// in place of the current one, never with an edit (see [`Session::opened`]).
    opened: u64,
    /// What the viewport hides (`view.*`): saved with the model, but no edit of it, so it is
    /// neither journalled nor undone, and it leaves results current.
    pub view: fieldcraft_io::View,
    /// Copies of the geometry files the model's file embeds, by content hash: opened where a file
    /// is missing, and embedded again when it saves packed.
    pub embedded: BTreeMap<String, Arc<Vec<u8>>>,
    /// The steps' preferred run settings (`solve.settings`): saved with the model, outside its
    /// definition, so neither journalled nor undone; the machine may override them
    /// (`FIELDCRAFT_SOLVE_BACKEND`, or a solve's own `backend`).
    pub run: fieldcraft_io::Run,
    /// The run settings as last saved (or opened): a change to them is unsaved too.
    saved_run: fieldcraft_io::Run,
    /// Whether saving embeds the geometry files (`file.save`'s `geometry`): a model opened from a
    /// packed file stays packed.
    pub pack_geometry: bool,
    /// The view as last saved (or opened): a change to it is unsaved too.
    saved_view: fieldcraft_io::View,
    /// A saved view waiting for its step's result to be shown (`view.restore` or opening a
    /// file before the step is solved).
    pending_view: Option<String>,
    /// The camera pictures are aimed with, in the display units: the one last aimed by a host
    /// (the app's viewport, an MCP picture; what `view.save` keeps), or asked for by a saved
    /// view ([`Session::request_camera`]). Session state, like the results view.
    pub camera: Option<view::Aim>,
    /// Counts the cameras asked for, so a host aims its own once per request.
    camera_requests: u64,
    /// This is the app's own model: its top-level commands go into the crash journal
    /// ([`crash`]), so an unclean exit can be replayed. Off for scratch sessions and tests.
    pub crash_journal: bool,
    /// Where this session's commands come from, for the usage analytics ([`usage`]): the app
    /// sets it on its own model (and the control channel and the command search while they
    /// run); `None` (scratch sessions, the CLI, tests) counts nothing.
    pub origin: Option<usage::Origin>,
    /// What the last `geometry.check` found, drawn as markers while the design it checked is
    /// open ([`view::geometry_marks`]). Not an edit: neither journalled, saved nor undone.
    pub geometry_marks: Option<Arc<view::GeometryMarks>>,
    /// The streamlines last traced for a fluid step's display ([`flow::TraceCache`]): the
    /// `results.streamlines` reply and the plot share them. Not an edit either.
    pub traced: flow::TraceCache,
}

impl Session {
    /// Run a command. On an error or a panic the model, the design and the surfaces are left as
    /// they were.
    pub fn execute(&mut self, id: &str, params: &Value) -> Result<Value> {
        // An id from hostile input is not echoed in full.
        let spec = find_command(id).ok_or_else(|| EngineError::UnknownCommand(id.chars().take(100).collect()))?;
        let timer = usage::Timer::start(self.origin);
        let elements = self.model.elements().len();
        let result = self
            .availability(spec)
            .map_err(EngineError::Unavailable)
            .and_then(|()| check_keys(spec, params))
            .and_then(|()| self.run_command(spec, params));
        if timer.counting() {
            let milestone = usage::milestone(spec.id, (elements, self.model.elements().len()));
            let sample =
                (spec.id == "file.sample").then(|| params.get("name").and_then(Value::as_str).and_then(sample::find)).flatten().map(|s| s.name);
            timer.done(spec.id, spec.kind, result.as_ref().err(), milestone, sample);
        }
        result
    }

    /// [`Session::execute`] once the command is found and can run.
    fn run_command(&mut self, spec: &'static CommandSpec, params: &Value) -> Result<Value> {
        let id = spec.id;
        // Session commands (open, save, undo, redo) change the session only in their last step,
        // after everything that can fail, so they need no copy of the model to restore.
        let before = (spec.kind == CommandKind::Model).then(|| (self.snapshot(), self.document()));
        let scope = crash::Scope::enter();
        let journaled = crash::command_begin(&scope, self, id, params);
        let (revision, opened) = (self.revision, self.opened);
        // Before format 3 the geometry and mesher commands read their numbers as millimetres in a
        // model built on a design or 3MF meshes, whatever it declared: a replayed journal from
        // then is read so. One that acted on a block typed in another unit cannot be replayed
        // exactly, and is refused rather than guessed.
        let legacy_geometry = self.legacy_replay && units::frame(id) == units::Frame::Geometry && self.on_geometry();
        if legacy_geometry && let Some(why) = self.legacy_block_target(id, params) {
            return Err(EngineError::Params(format!("this {id} from FieldCraft 0.7 acted on {why}; its numbers cannot be replayed exactly")));
        }
        let input = if legacy_geometry { units::geometry().clone() } else { self.input_units().clone() };
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let target = match units::frame(id) {
                units::Frame::Si => units::si().clone(),
                units::Frame::Geometry => units::geometry().clone(),
                units::Frame::Input => input.clone(),
            };
            let converted = match spec.form {
                Some(form) => units::convert_params(form, params, &input, &target)?,
                None => params.clone(),
            };
            let shown = self.units().clone();
            let mut reply = self
                .in_frame(&target, |s| (spec.run)(s, &converted).and_then(|v| cmd::relax::after_remesh(s, id, v)))
                .map_err(|e| e.shown(&target, &shown))?;
            units::convert_reply(id, &mut reply, &target, &shown)?;
            Ok(reply)
        }));
        let result = match outcome {
            Ok(Err(EngineError::Param(p))) => Err(EngineError::Param(match spec.form {
                Some(form) => schema::refine(form, params, p),
                None => p,
            })),
            Ok(r) => r,
            Err(p) => {
                let message = p.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| p.downcast_ref::<String>().cloned()).unwrap_or_default();
                log::error!("command {id} panicked: {message}");
                Err(EngineError::Internal { command: id.to_string(), message })
            }
        };
        if let Some((mut before, document)) = before {
            // A command that makes the session another document (a new model) keeps the old
            // one in its undo step.
            if self.document() != document {
                before.document = Some(document);
            }
            if result.is_err() {
                self.restore(before);
            } else if before.document.is_some() || self.differs_from(&before) {
                if self.journal.is_empty() {
                    self.journal_header.units = input.clone();
                }
                // Bare numbers the journal would read otherwise say their units.
                let entry = match spec.form {
                    Some(form) if input != *self.journal_reads(id) => units::with_units(form, params, &input),
                    _ => params.clone(),
                };
                self.journal.push((id.to_string(), entry));
                self.revision = next_revision();
                self.record_undo(before);
                // An item deleted (or a body cut into cells) is no longer hidden.
                cmd::view::forget_missing(self);
            }
        }
        if self.revision != revision {
            self.refresh_results();
        }
        drop(scope);
        crash::command_done(self, journaled, self.revision != revision || self.opened != opened, id, params, &result);
        result
    }

    /// Run a script: a JSON array of `{"command": id, "params": {...}}` (bare numbers in the
    /// session's display system), or a journal's form, `{"units": system, "model_units"?:
    /// system, "commands": [...]}`, whose bare numbers are all read in `units` (and whose model,
    /// when `model_units` is given and nothing was done yet, starts shown in that system).
    /// Returns each step's output, in order. Other fields of a step (such as an `id` naming it)
    /// are left to the caller.
    pub fn run_script(&mut self, script: &Value) -> Result<Vec<Value>> {
        let (steps, header) = match script {
            Value::Array(steps) => (steps, None),
            Value::Object(o) => {
                let steps = o
                    .get("commands")
                    .and_then(Value::as_array)
                    .ok_or_else(|| EngineError::Params("a script with a header lists its `commands`".into()))?;
                let mut head = o.clone();
                head.remove("commands");
                let header: fieldcraft_io::JournalHeader =
                    serde_json::from_value(Value::Object(head)).map_err(|e| EngineError::Params(format!("the script's header: {e}")))?;
                (steps, Some(header))
            }
            _ => return Err(EngineError::Params("a script is a JSON array of steps, or {\"units\", \"commands\"}".into())),
        };
        let pinned = self.pinned_input.clone();
        let was_legacy = self.legacy_replay;
        // A journal from before format 3 names the units its model was in; format 2's defaults
        // were in its own units.
        let legacy = header.as_ref().filter(|h| h.model_units.is_some()).map(|h| h.units.preset().label());
        self.legacy_replay = legacy.is_some() || was_legacy;
        if let Some(h) = &header {
            if self.journal.is_empty() && h.model_units.is_some() && self.model == fieldcraft_model::Model::default() {
                self.model.set_display_units(h.model_units.clone());
                self.journal_header.model_units = h.model_units.clone();
            }
            self.pinned_input = Some(h.units.clone());
        }
        let mut outputs = Vec::with_capacity(steps.len());
        let mut failed = None;
        for (i, step) in steps.iter().enumerate() {
            let Some(id) = step.get("command").and_then(Value::as_str) else {
                failed = Some(EngineError::Params(format!("step {}: no `command`", i + 1)));
                break;
            };
            let mut params = step.get("params").cloned().unwrap_or(Value::Object(Default::default()));
            // A journal from before format 3 (it names the units its model was in) left out what
            // format 2's defaults filled in, in its own units.
            if let Some(system) = legacy {
                fieldcraft_io::migrate::format_2_defaults(id, &mut params, system);
            }
            match self.execute(id, &params) {
                Ok(v) => outputs.push(v),
                Err(e) => {
                    failed = Some(EngineError::Params(format!("step {} ({id}): {e}", i + 1)));
                    break;
                }
            }
        }
        self.pinned_input = pinned;
        self.legacy_replay = was_legacy;
        match failed {
            Some(e) => Err(e),
            None => Ok(outputs),
        }
    }

    /// The display system: the model's own, else the session's default.
    pub fn units(&self) -> &fieldcraft_model::DisplaySystem {
        self.model.display_units().unwrap_or(&self.default_units)
    }

    /// The system bare numbers of commands are read in: a replaying script's, else the display
    /// system.
    /// A journal entry's parameters in the display units, for a form to show and run again: the
    /// journal's bare numbers are in its own units. Left as they are when they cannot be converted
    /// (a command without a form has no numbers that follow the units).
    pub fn journal_params_shown(&self, command: &str, params: &Value) -> Value {
        match find_command(command).and_then(|spec| spec.form) {
            Some(form) => units::convert_params(form, params, self.journal_reads(command), self.units()).unwrap_or_else(|_| params.clone()),
            None => params.clone(),
        }
    }

    pub fn input_units(&self) -> &fieldcraft_model::DisplaySystem {
        self.pinned_input.as_ref().unwrap_or_else(|| self.units())
    }

    /// `params` of `command` with every bare number given its unit in the input system: they
    /// mean the same when run again after the display units change (a reply the UI refreshes).
    pub fn params_with_units(&self, command: &str, params: &Value) -> Value {
        match find_command(command).and_then(|spec| spec.form) {
            Some(form) => units::with_units(form, params, self.input_units()),
            None => params.clone(),
        }
    }

    /// Whether the model is built on the geometry: a design, or imported 3MF meshes.
    fn on_geometry(&self) -> bool {
        self.design.is_some() || !self.surfaces.is_empty()
    }

    /// The system a replay of the journal reads the bare numbers of command `id` in: the journal's
    /// own, except in a journal from before format 3 of a model built on the geometry, whose
    /// geometry and mesher commands read them in millimetres, as they did then.
    fn journal_reads(&self, id: &str) -> &fieldcraft_model::DisplaySystem {
        if self.journal_header.model_units.is_some() && units::frame(id) == units::Frame::Geometry && self.on_geometry() {
            units::geometry()
        } else {
            &self.journal_header.units
        }
    }

    /// A block body (typed in the journal's own units, not millimetres) a replayed format 2
    /// geometry command `id` acts on: named in `body`/`bodies`, holding one of its `nodes`, with a
    /// node in its sphere (`center`, `radius`, in millimetres as it was then), or, for a command
    /// over the whole mesh given no target, any block beside a mesh of the design.
    fn legacy_block_target(&self, id: &str, params: &Value) -> Option<String> {
        let unit = self.input_units().symbol(fieldcraft_model::units::Quantity::Length).to_string();
        if unit == "mm" {
            return None;
        }
        let journal: Vec<Value> = self.journal.iter().map(|(c, p)| serde_json::json!({"command": c, "params": p})).collect();
        let blocks = fieldcraft_io::migrate::format_2_typed_bodies(&Value::Array(journal));
        let named: Vec<&str> = match (params.get("body"), params.get("bodies")) {
            (Some(Value::String(b)), _) => vec![b.as_str()],
            (_, Some(Value::Array(a))) => a.iter().filter_map(Value::as_str).collect(),
            _ => Vec::new(),
        };
        if let Some(b) = named.iter().find(|b| blocks.iter().any(|x| x == *b)) {
            return Some(format!("block body `{b}`, typed in {unit}"));
        }
        // A command over the whole mesh (no body, nodes, region or sphere) acts on the blocks too
        // when the mesh has both kinds.
        const WHOLE_MESH: [&str; 4] = ["mesh.node.relax", "mesh.node.relax.apply", "mesh.node.merge", "mesh.mirror"];
        const TARGETS: [&str; 10] = ["body", "bodies", "nodes", "elements", "faces", "edges", "set", "center", "around_element", "geometry"];
        if WHOLE_MESH.contains(&id)
            && !TARGETS.iter().any(|k| params.get(*k).is_some_and(|v| !v.is_null()))
            && let Some(b) = blocks.iter().find(|b| !self.model.body_elements(b).is_empty())
        {
            let meshed: Vec<String> = self.model.classification().entities().iter().map(|e| e.body.clone()).collect();
            if meshed.iter().any(|m| !blocks.contains(m) && !self.model.body_elements(m).is_empty()) {
                return Some(format!("every body of the mesh, the block body `{b}` (typed in {unit}) among them beside meshes of the design"));
            }
        }
        let given: Vec<u32> = params
            .get("nodes")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|n| n.as_u64())
            .filter_map(|n| u32::try_from(n).ok())
            .collect();
        let centre = params.get("center").and_then(|c| serde_json::from_value::<[f64; 3]>(c.clone()).ok());
        let radius = params.get("radius").and_then(Value::as_f64);
        for b in &blocks {
            for e in self.model.body_elements(b) {
                for n in self.model.element(e).map(|e| e.nodes.clone()).unwrap_or_default() {
                    let in_sphere = match (centre, radius, self.model.node(n)) {
                        (Some(c), Some(r), Some(p)) => {
                            let c = units::geometry_point(c);
                            let r = units::geometry_length(r);
                            (0..3).map(|k| (p.xyz[k] - c[k]).powi(2)).sum::<f64>() <= r * r
                        }
                        _ => false,
                    };
                    if given.contains(&n) || in_sphere {
                        return Some(format!("node {n} of block body `{b}`, typed in {unit}"));
                    }
                }
            }
        }
        None
    }

    /// The model expressed in `frame`, as a command run there sees it: a function of the stored
    /// model alone, never of what ran before (so a replay sees what the live session saw, and a
    /// preview what the command will).
    pub(crate) fn framed_model(&self, frame: &fieldcraft_model::DisplaySystem) -> Result<fieldcraft_model::Model> {
        Ok(self.model.converted(&|q, v| frame.express(q, v).ok())?)
    }

    /// Run `f` with the model expressed in `frame` (unchanged when it is SI): what the command
    /// leaves changed comes back to SI, what it left alone keeps its exact SI value.
    fn in_frame(&mut self, frame: &fieldcraft_model::DisplaySystem, f: impl FnOnce(&mut Session) -> Result<Value>) -> Result<Value> {
        if frame == units::si() {
            return f(self);
        }
        let stored = self.model.clone();
        let framed = self.framed_model(frame)?;
        self.model = framed.clone();
        let result = f(self);
        let after = std::mem::take(&mut self.model);
        match after.from_converted(&framed, &stored, &|q, v| frame.to_si(q, v).ok()) {
            Ok(m) => self.model = m,
            Err(e) => {
                self.model = stored;
                return Err(e.into());
            }
        }
        result
    }

    /// Identifies the open document: it changes when `file.open`, `file.recover`, `file.sample`
    /// or `deck.read` replaces the session and when `geometry.open` opens another design file or starts a new
    /// model with one; an edit, a solve, undo, redo or opening the same design again after an
    /// edit keep it. A view frames the model again when it changes.
    pub fn opened(&self) -> u64 {
        self.opened
    }

    /// Make the session another document (see [`Session::opened`]).
    pub(crate) fn mark_opened(&mut self) {
        self.opened = next_revision();
    }

    /// Identifies the current state (see [`Session::is_dirty`]).
    pub fn revision(&self) -> u64 {
        self.revision
    }

    /// Whether the session has changes its file does not hold. Undoing back to the saved state
    /// makes it clean again.
    pub fn is_dirty(&self) -> bool {
        self.saved_revision.unwrap_or(0) != self.revision
            || self.view != self.saved_view
            || self.results_stored != self.saved_results
            || self.run != self.saved_run
    }

    /// The stored result of a step, if it was solved in the current state (an edit makes it stale;
    /// undoing back to the solved state makes it current again).
    pub fn result(&self, step: &str) -> Option<&Arc<StepResult>> {
        self.results.get(step).filter(|(revision, _)| *revision == self.revision).map(|(_, r)| r)
    }

    /// The steps with a current result.
    pub fn solved_steps(&self) -> impl Iterator<Item = &str> {
        self.results.iter().filter(|(_, (revision, _))| *revision == self.revision).map(|(name, _)| name.as_str())
    }

    /// Keep a step's result for the current state. With no results shown, the display shows it;
    /// a view of the same step stays as it is (on mode 1 and the last frame if its own are gone).
    pub(crate) fn store_result(&mut self, result: StepResult) {
        match &mut self.results_view {
            None => self.results_view = Some(results::ResultsView::for_result(&result)),
            Some(v) if v.step == result.step() && results::displacements(&result, v.pick()).is_err() => {
                v.mode = 1;
                v.frame = None;
            }
            Some(_) => {}
        }
        let step = result.step().to_string();
        match fieldcraft_io::results::definition_hash(&self.model, &step, &self.sources) {
            Some(f) => self.result_inputs.insert(step.clone(), f),
            None => self.result_inputs.remove(&step),
        };
        self.results.insert(step.clone(), (self.revision, Arc::new(result)));
        self.results_stored = self.results_stored.wrapping_add(1);
        cmd::views::show_pending(self, &step);
    }

    /// Keep a result read from a file, solved from the definition `definition`: current when that
    /// is the model's, else stale (kept, and current again should an edit bring the model back to
    /// it). Opening is not solving: nothing is unsaved.
    pub(crate) fn keep_stored_result(&mut self, result: StepResult, definition: [u8; 32]) -> bool {
        let step = result.step().to_string();
        let current = fieldcraft_io::results::definition_hash(&self.model, &step, &self.sources) == Some(definition);
        let revision = if current { self.revision } else { next_revision() };
        self.result_inputs.insert(step.clone(), definition);
        self.results.insert(step, (revision, Arc::new(result)));
        current
    }

    /// Every step's latest result with the definition it was solved from, as files keep them (a
    /// result whose step is gone is kept with a definition no model has).
    pub(crate) fn stored_results(&self) -> Vec<(&str, fieldcraft_io::results::StoredRef<'_>)> {
        self.results
            .iter()
            .map(|(step, (_, r))| {
                let definition = self.result_inputs.get(step).copied().unwrap_or([0; 32]);
                (step.as_str(), fieldcraft_io::results::StoredRef { definition, result: r.as_stored() })
            })
            .collect()
    }

    /// After the state changed (an edit, undo or redo): a result whose step reads what it read
    /// when it was solved describes the new state too, so it stays (or is again) current; the
    /// others go stale.
    fn refresh_results(&mut self) {
        let revision = self.revision;
        for (step, (solved, _)) in self.results.iter_mut().filter(|(_, (r, _))| *r != revision) {
            let read = self.result_inputs.get(step).copied();
            if read.is_some() && read == fieldcraft_io::results::definition_hash(&self.model, step, &self.sources) {
                *solved = revision;
            }
        }
    }

    /// The path to use for `access` to the file at `path`, checked against [`Session::file_access`].
    pub fn checked_path(&self, path: &str, access: fieldcraft_io::roots::Access) -> Result<String> {
        let checked = self.file_access.check(std::path::Path::new(path), access)?;
        match checked.to_str() {
            Some(p) => Ok(p.to_string()),
            None => Err(EngineError::Params(format!("`{path}` resolves to a path that is not valid UTF-8"))),
        }
    }

    /// Ask the hosts to aim their camera as `aim` says (a saved view shown): it becomes the
    /// session's camera and [`Session::camera_requests`] counts one more.
    pub fn request_camera(&mut self, aim: view::Aim) {
        self.camera = Some(aim);
        self.camera_requests = self.camera_requests.wrapping_add(1);
    }

    /// An empty session to open a whole model into off the host's thread: this one's work caps,
    /// file access, crash journal, command origin and display units, and its count of camera
    /// requests, so a view the opened model shows is a request the host has not seen yet.
    pub fn empty_like(&self) -> Session {
        Session {
            limits: self.limits,
            file_access: self.file_access.clone(),
            crash_journal: self.crash_journal,
            origin: self.origin,
            default_units: self.default_units.clone(),
            camera_requests: self.camera_requests,
            ..Session::default()
        }
    }

    /// How many cameras were asked for: a host that saw fewer aims at [`Session::camera`].
    pub fn camera_requests(&self) -> u64 {
        self.camera_requests
    }

    /// The current state is what `path` now holds.
    pub(crate) fn mark_saved(&mut self, path: String) {
        self.path = Some(path);
        self.migrated_from = None;
        self.migration_notes.clear();
        self.saved_revision = Some(self.revision);
        self.saved_view = self.view.clone();
        self.saved_results = self.results_stored;
        self.saved_run = self.run.clone();
    }
}

#[cfg(test)]
mod cancel_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_geometry::sample::{Plate, plate_design};
    use serde_json::json;

    #[test]
    fn command_ids_are_unique_and_well_formed() {
        let mut ids: Vec<&str> = commands().iter().map(|c| c.id).collect();
        for id in &ids {
            assert!(id.split('.').count() >= 2 && id.chars().all(|c| c.is_ascii_lowercase() || c == '.' || c == '_'), "{id}");
        }
        ids.sort_unstable();
        let n = ids.len();
        ids.dedup();
        assert_eq!(ids.len(), n);
    }

    #[test]
    fn script_builds_a_hex_and_failures_leave_the_model_alone() {
        let mut s = Session::default();
        let mut script = Vec::new();
        for (i, c) in [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]].iter().enumerate() {
            script.push(json!({"command": "mesh.node.add", "params": {"label": i + 1, "xyz": c}}));
        }
        script.push(json!({"command": "mesh.element.add", "params": {"label": 1, "kind": "hex8", "nodes": [1,2,3,4,5,6,7,8]}}));
        s.run_script(&Value::Array(script)).unwrap();
        assert_eq!(s.journal.len(), 9);
        let info = s.execute("model.info", &json!({})).unwrap();
        assert_eq!(info["nodes"], 8);
        assert_eq!(info["elements"]["hex8"], 1);

        let model = s.model.clone();
        assert!(s.execute("mesh.element.add", &json!({"label": 2, "kind": "tet4", "nodes": [1, 2, 3, 99]})).is_err());
        assert!(s.execute("mesh.node.remove", &json!({"label": 1})).is_err());
        assert_eq!(s.model, model);
        assert_eq!(s.journal.len(), 9);
    }

    #[test]
    fn hostile_params_never_panic() {
        let hostile = [
            json!(null),
            json!([]),
            json!("x"),
            json!({}),
            json!({"label": -1}),
            json!({"label": 1e300}),
            json!({"label": u64::MAX}),
            json!({"label": 1, "xyz": [1, 2]}),
            json!({"label": 1, "xyz": ["a", null, {}]}),
            json!({"label": 1, "xyz": [1e308, -1e308, 0]}),
            json!({"label": 1, "kind": "hex8", "nodes": null}),
            json!({"label": 1, "kind": "nope", "nodes": [1]}),
            json!({"label": 1, "kind": "tet4", "nodes": [1, 1, 1, 1]}),
            json!({"label": 1, "kind": "hex20", "nodes": vec![u32::MAX; 20]}),
            json!({"path": ""}),
            json!({"path": "/"}),
            json!({"path": "/dev/zero"}),
            json!({"path": "/no/such/design.solidcraft"}),
            json!({"path": "x".repeat(100_000)}),
            json!({"format": "abaqus"}),
            json!({"format": "abaqus", "path": ""}),
            json!({"format": "abaqus", "path": "/"}),
            json!({"format": "abaqus", "path": "/dev/zero"}),
            json!({"format": "abaqus", "path": "/no/such/part.inp", "units": "guess"}),
            json!({"format": "abaqus", "path": "/no/such/dir/part.inp"}),
            json!({"format": "abaqus", "path": format!("/no/such/{}", "x".repeat(100_000))}),
            json!({"format": "abaqus", "path": 7}),
            json!({"format": "x".repeat(100_000), "path": "/no/such/part.inp"}),
            json!({"format": "lsdyna", "path": "/no/such/dir/part.k"}),
            json!({"format": "lsdyna", "path": "/"}),
            json!({"format": "lsdyna", "path": "x\u{0}.k"}),
            json!({"format": null, "path": "/no/such/part.inp"}),
            json!({"format": "nastran", "path": "/no/such/dir/part.bdf"}),
            json!({"format": "nastran", "path": "/no/such/dir/part.bdf", "steps": ["Bend"]}),
            json!({"format": "nastran", "path": "/no/such/dir/part.bdf", "steps": [null, 3, "Modes"]}),
            json!({"format": "nastran", "path": "/no/such/dir/part.bdf", "steps": vec!["x".repeat(1000); 2000]}),
            json!({"format": "nastran", "path": "/no/such/dir/part.bdf", "steps": {}}),
            json!({"path": 7}),
            json!({"path": "/no/such/dir/out.vtu", "step": "Nope"}),
            json!({"path": "/no/such/dir/out.vtu", "mode": 0}),
            json!({"path": "/no/such/dir/out.pvd", "mode": 1_000_000_000}),
            json!({"path": "/no/such/dir/out.pvd", "frames": 3, "encoding": []}),
            json!({"path": "x\u{0}.vtu"}),
            json!({"path": format!("/no/such/{}.vtu", "x".repeat(100_000))}),
            json!({"path": "x\u{0}.stl"}),
            json!({"path": "/no/such/dir/skin.stl", "set": 7}),
            json!({"path": "/no/such/dir/skin.stl", "set": "x".repeat(100_000)}),
            json!({"path": "/no/such/dir/skin.stl", "format": ["ascii"]}),
            // (`max_points` and `max_passes` keep mesh.tet.generate's quality stage small.)
            json!({"body": "Body1", "tolerance": 0.05, "max_points": 3000, "max_passes": 1}),
            json!({"body": null, "tolerance": "fine"}),
            json!({"elements": [u32::MAX, -1, 1.5]}),
            json!({"elements": {}, "bins": -5}),
            json!({"elements": [1], "bins": 1e300}),
            json!({"points": null}),
            json!({"points": [[1, 2]]}),
            json!({"points": [[1e308, -1e308, 0], [0, 0, 0]]}),
            json!({"points": [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]}),
            json!({"points": vec![[1, 1, 1]; 30_000]}),
            json!({"surface": 0}),
            json!({"surface": -1}),
            json!({"surface": 1e300}),
            json!({"surface": "0", "weld": -1, "feature_angle": 1e300}),
            json!({"polygon": vec![[0.0, 0.0]; 200_000], "size": 1}),
            json!({"polygon": [[0, 0], [1e308, 0], [0, 1e308]], "size": 1e-300}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "holes": vec![[[0.1, 0.1], [0.2, 0.1], [0.1, 0.2]]; 20_000], "size": 0.5}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "holes": [[[0.1, 0.1], [0.2, 0.1], [0.1, 0.2]], [[0.1, 0.1], [0.2, 0.1], [0.1, 0.2]]], "size": 0.2}),
            json!({"polygon": [["a", 0], [1, null], {}], "size": "big", "kind": 7}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 1e300, "kind": "tri6", "max_elements": 1e300}),
            json!({"body": "Body1", "face": "", "tolerance": 0.05, "size": 1, "max_points": 3000, "max_passes": 1}),
            // Structured hex meshing.
            json!({"size": -1, "order": 9, "method": "x"}),
            json!({"size": 1e-300}),
            json!({"order": 2, "method": "mapped", "tolerance": 1e300}),
            json!({"surface": 0, "size": 1e300, "order": 2}),
            // Model definitions.
            json!({"name": "x".repeat(100_000)}),
            json!({"name": "A\u{0}"}),
            json!({"name": 7, "nodes": [1]}),
            json!({"name": "A", "nodes": [0, -1, 1e300]}),
            json!({"name": "A", "nodes": vec![u32::MAX; 3], "faces": [{"element": 1, "face": 0}]}),
            json!({"name": "A", "geometry": [{"body": ""}, {"body": "B", "face": "f", "edge": "e"}]}),
            json!({"name": "A", "geometry": [{"body": "x".repeat(10_000)}]}),
            json!({"name": "A", "elastic": {"type": "isotropic", "young": 1e308, "poisson": 0.49999999}}),
            json!({"name": "A", "elastic": {"type": "engineering_constants", "young": [1e-300, 1e300, 1], "poisson": [1e300, 0, 0], "shear": [0, 0, 0]}}),
            json!({"name": "A", "elastic": {"type": "isotropic", "young": "stiff"}, "density": -0.0}),
            json!({"name": "A", "set": "A", "material": "A", "kind": {"type": "solid"}}),
            json!({"name": "Steel", "plastic": {"type": "bilinear", "yield_stress": -1e308, "hardening_modulus": f64::MAX}}),
            json!({"name": "Steel", "plastic": {"type": "table", "points": vec![json!({"strain": 0, "stress": 1}); 2_000]}}),
            json!({"name": "Steel", "plastic": {"type": "table", "points": [{"strain": 0, "stress": 1}, {"strain": 0, "stress": 2}]}}),
            json!({"name": "Steel", "plastic": {"type": "johnson_cook", "a": 1, "b": 1, "n": -1, "room_temperature": 1e308, "melt_temperature": -1e308}}),
            json!({"name": "Steel", "plastic": {"type": "johnson_cook", "a": 1, "b": 1e308, "n": 1e308, "c": 1e308, "m": 1e-308, "room_temperature": 0, "melt_temperature": 1, "taylor_quinney": 2}}),
            json!({"name": "Steel", "plastic": "hard"}),
            json!({"name": "Steel", "damage": {"type": "johnson_cook", "d1": 1e308, "d2": -1e308, "d3": 1e308}, "erosion": {"plastic_strain": -1}}),
            json!({"name": "Steel", "damage": null, "erosion": {"plastic_strain": 1e-308}}),
            json!({"name": "A", "set": "all", "material": "Steel", "kind": {"type": "beam", "profile": {"type": "pipe", "radius": 1, "thickness": 2}, "orientation": [0, 0, 0]}}),
            json!({"name": "A", "set": "all", "material": "Steel", "kind": {"type": "shell", "thickness": 1e308, "offset": -0.5, "points": 4294967295u64}}),
            json!({"name": "A", "xyz": [1e308, 0, 0]}),
            json!({"name": "A", "set": "all", "material": "Steel", "kind": {"type": "plane_stress", "thickness": -1}}),
            json!({"name": "A", "set": "all", "material": "Steel", "kind": {"type": "plane_strain", "thickness": 1e308}}),
            json!({"name": "A", "set": "all", "material": "Steel", "kind": {"type": "axisymmetric", "axis": "x"}}),
            json!({"name": "A", "type": "centrifugal", "angular_velocity": 1e300, "point": [0, 0, 0], "axis": [0, 1, 0]}),
            json!({"name": "A", "type": "centrifugal", "angular_velocity": 1, "point": [0, 0, 0], "axis": [0, 0, 0]}),
            json!({"name": "A", "edges": [{"element": u32::MAX, "edge": 0}]}),
            json!({"name": "A", "reference_point": "RP", "set": "end", "kind": "rigid", "dofs": [0, 7]}),
            json!({"name": "A", "type": "fixed", "target": {"set": "nope"}}),
            json!({"name": "A", "type": "displacement", "target": {"reference_point": "RP"}, "values": [null, null, null, null, null, null]}),
            json!({"name": "A", "type": "displacement", "target": {"set": "root"}, "values": [1e308, -1e308, 0, 0, 0, 0, 0]}),
            json!({"name": "A", "type": "pressure", "set": "end", "magnitude": -1e308}),
            json!({"name": "A", "type": "gravity", "acceleration": [0, 0, 0]}),
            json!({"name": "A", "type": "force", "target": {"set": "end"}, "vector": [1, 0, 0], "distribution": "lumpy"}),
            json!({"name": "A", "type": "traction", "set": "root", "vector": [1e308, 1e308, 1e308]}),
            json!({"name": "A", "type": "traction", "set": "nope", "vector": [1, 0, 0]}),
            json!({"name": "RP", "xyz": [0, 0, 0], "label": 4294967295u64}),
            json!({"name": "A", "xyz": [0, 0, 0], "label": 1}),
            json!({"name": "A", "type": "body_force", "set": "all", "vector": [1, 2]}),
            json!({"name": "A", "procedure": {"type": "frequency", "modes": 4294967295u64, "min_hz": -1}}),
            json!({"name": "A", "procedure": {"type": "frequency", "modes": 5}, "output": ["everything"], "before": 3}),
            json!({"name": "Bend", "procedure": {"type": "static"}, "loads": "push"}),
            json!({"step": "Bend", "loads": ["nope"], "active": "yes"}),
            json!({"step": "Bend", "boundary_conditions": [], "loads": []}),
            json!({"name": "A", "target": {"set": "end"}, "mass": [1e308, -1, f64::MAX], "inertia": [1, 1, 1, 1e300, 0, 0]}),
            json!({"name": "A", "target": {"reference_point": "RP"}, "inertia": [1, 2]}),
            json!({"name": "A", "target": {"set": "end"}, "mass": "heavy"}),
            json!({"name": "A", "type": "spring", "first": {"set": "end"}, "second": {"set": "end"}, "direction": "line", "stiffness": 1e308}),
            json!({"name": "A", "type": "spring", "first": {"set": "end"}, "direction": {"dof": 255}, "stiffness": 1}),
            json!({"name": "A", "type": "dashpot", "first": {"set": "end"}, "direction": {"vector": [1e-320, 0, 0]}, "coefficient": -1}),
            json!({"name": "A", "type": "dashpot", "first": {"reference_point": "RP"}, "direction": {"dof": 4}, "coefficient": 1e300}),
            json!({"system": "furlongs"}),
            // Fluid domains, boundaries and steps.
            json!({"name": "A", "material": "Steel", "min": [0, 0, 0], "max": [1e308, 1, 1], "spacing": 1e-300}),
            json!({"name": "A", "material": "Steel", "min": [0, 0, 0], "max": [1, 1, 1], "spacing": 0.5, "faces": {"x_min": {"type": "periodic"}}}),
            json!({"name": "A", "material": "nope", "min": [0, 0, 0], "max": [1, 1, 1], "spacing": 0.5, "obstacles": ["", ""]}),
            json!({"name": "A", "density": 1, "fluid": {"type": "newtonian", "viscosity": -1}}),
            json!({"name": "A", "fluid": {"type": "newtonian", "viscosity": 1}}),
            json!({"domain": "A", "at": "x_mid", "condition": {"type": "wall"}}),
            json!({"domain": "A", "at": {"body": "B", "face": "f"}, "condition": {"type": "pressure", "pressure": 1e308}}),
            json!({"domain": null, "at": null, "condition": null}),
            json!({"name": "A", "procedure": {"type": "fluid", "domain": "A", "steps": 1e300}}),
            json!({"name": "A", "procedure": {"type": "fluid", "domain": "A", "end_time": -1, "steps": 3}}),
            json!({"name": "Bend", "procedure": {"type": "fluid", "domain": "nope", "end_time": 1}}),
            json!({"step": "Bend", "max_cells": 1e300, "max_steps": -3}),
            // Heat transfer.
            json!({"name": "A", "thermal": {"conductivity": -1}}),
            json!({"name": "A", "thermal": {"conductivity": [1e308, f64::MIN_POSITIVE, 0]}, "density": 1}),
            json!({"name": "A", "thermal": {"conductivity": "hot", "specific_heat": 0, "expansion": [1, 2]}}),
            json!({"name": "A", "thermal": {}}),
            json!({"name": "A", "thermal": {"expansion": 1e308, "reference_temperature": -1e308}}),
            json!({"name": "A", "type": "temperature", "target": {"reference_point": "RP"}, "value": 1}),
            json!({"name": "A", "type": "temperature", "target": {"set": "root"}, "value": "warm"}),
            json!({"name": "A", "type": "surface_flux", "set": "end", "flux": 0}),
            json!({"name": "A", "type": "convection", "set": "end", "film": -1, "sink": 1e308}),
            json!({"name": "A", "type": "heat_source", "set": "all", "power": -1e308}),
            json!({"name": "A", "procedure": {"type": "heat_transfer", "transient": {"end_time": 1e308, "increment": 1e-308}}}),
            json!({"name": "A", "procedure": {"type": "heat_transfer", "transient": {"end_time": 1, "increment": 0.1, "theta": 0.1}}}),
            json!({"name": "A", "procedure": {"type": "heat_transfer", "transient": {"end_time": -1, "increment": 0}}}),
            json!({"name": "A", "procedure": {"type": "coupled_temperature_displacement", "extra": 1}}),
            json!({"step": "Heat", "frame": u64::MAX, "nodes": [1]}),
            // Contact.
            json!({"name": "C", "main": 7, "secondary": null}),
            json!({"name": "C", "main": "end", "secondary": "root", "friction": {"static": 1e308, "kinetic": -1, "decay": -1e308}}),
            json!({"name": "C", "main": "end", "secondary": "end"}),
            json!({"name": "C", "main": "x".repeat(10_000), "secondary": "root", "penalty_scale": 1e308}),
            json!({"name": "C", "surface": "root", "penalty_scale": 0, "initial_penetration": "ignore"}),
            json!({"name": "C", "surface": 3}),
            json!({"name": "C", "surfaces": {"type": "pair", "main": "end"}}),
            json!({"name": "C", "surfaces": {"type": "general", "surface": "nope"}, "friction": "sticky"}),
            json!({"tolerance": 1e300, "create": true, "friction": {"static": 0.1}}),
            json!({"tolerance": -0.0, "create": true}),
            json!({"tolerance": 1e-300, "create": 1}),
            json!({"geometry": {"body": "Bar"}, "elements": [1, 99], "faces": [{"element": 1, "face": 200}]}),
            json!({"geometry": {"body": "Bar", "vertex": "v"}, "nodes": [1], "edges": [{"element": 1, "edge": 0}]}),
        ];
        // Small work caps, so input that would mesh up to the production caps (needle triangles
        // need Steiner points all along them, an octree hex mesh is fitted and untangled) is
        // refused at once. The caps themselves are what
        // `tet_generate_honours_the_session_caps` checks.
        let capped =
            || Session { limits: Limits { max_steiner: 64, max_tets: 100_000, max_hex_elements: 500, ..Limits::default() }, ..Session::default() };
        for spec in commands() {
            for p in &hostile {
                let mut s = capped();
                let _ = s.execute(spec.id, &p.clone());
            }
        }
        assert!(capped().execute("no.such", &json!({})).is_err());

        // The same against a defined model, so edits and deletes reach existing definitions.
        let mut defined = capped();
        crate::cmd::define::tests::defined_bar(&mut defined);
        for spec in commands() {
            for p in &hostile {
                let mut s = defined.clone();
                let before = s.model.clone();
                if s.execute(spec.id, p).is_err() {
                    assert_eq!(s.model, before, "{} {p}", spec.id);
                }
            }
        }

        // The same against an open design.
        let (_dir, path) = sample_design_file(&Plate::default());
        let mut opened = capped();
        opened.execute("geometry.open", &json!({"path": path})).unwrap();
        let tess = [
            json!({"body": "Body1"}),
            json!({"body": "Body1", "tolerance": 0}),
            json!({"body": "Body1", "tolerance": -1}),
            json!({"body": "Body1", "tolerance": 1e-300}),
            json!({"body": "Body1", "tolerance": 1e300}),
            json!({"body": "Body2", "tolerance": 0.1}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "surface": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "weld": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "weld": 1e-3}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "weld": f64::MIN_POSITIVE}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "weld": 0, "weld_distance": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "weld_distance": 1e300}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "weld_distance": 0.05}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "feature_angle": 1e-300}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1, "feature_angle": null}),
            // The quality stage's own parameters, hostile, with small caps so this stays fast.
            json!({"body": "Body1", "tolerance": 0.5, "size": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": -1}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 1e300}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 1e-300, "max_points": 2000, "max_tets": 5000}),
            json!({"body": "Body1", "tolerance": 0.5, "size": "fine"}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "order": "cubic"}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "order": "quadratic", "quality": false}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "quality": "yes"}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "quality_target": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "quality_target": 1e300}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "curvature_angle": -1}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "curvature_angle": 1e-300, "max_points": 2000}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "grading": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "grading": f64::MAX}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_passes": 1e9}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_points": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 20, "max_points": 30, "max_passes": 1}),
            json!({"body": "Body1", "tolerance": 0.5, "size": 0.5, "max_tets": 3000, "max_passes": 1}),
            json!({"body": "Body1", "size": 1e-300}),
            json!({"body": "Body1", "size": 1e300, "order": 2}),
            json!({"body": "Body1", "size": 4, "method": "sweep", "order": 2}),
            json!({"body": "Body1", "method": "revolve"}),
        ];
        for p in tess.iter().chain(&hostile) {
            for spec in commands() {
                let mut s = opened.clone();
                let _ = s.execute(spec.id, p);
            }
        }
    }

    /// A design file in a fresh temporary directory (removed when the guard drops).
    pub(crate) fn sample_design_file(plate: &Plate) -> (TempDir, String) {
        design_file(&plate_design(plate).unwrap())
    }

    /// Design file bytes written to a temporary directory.
    pub(crate) fn design_file(bytes: &[u8]) -> (TempDir, String) {
        let dir = TempDir::new();
        let path = dir.0.join("design.solidcraft");
        std::fs::write(&path, bytes).unwrap();
        (dir, path.to_string_lossy().into_owned())
    }

    pub(crate) struct TempDir(std::path::PathBuf);
    impl TempDir {
        fn new() -> Self {
            static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let dir = std::env::temp_dir().join(format!("fieldcraft-engine-{}-{n}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            TempDir(dir)
        }
    }
    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn open_a_design_list_its_bodies_and_tessellate_one() {
        let (_dir, path) = sample_design_file(&Plate::default());
        let mut s = Session::default();
        let opened = s.execute("geometry.open", &json!({"path": path})).unwrap();
        assert_eq!(opened["bodies"], json!(["Body1"]));
        assert_eq!(opened["errors"], json!([]));
        assert_eq!(opened["warnings"], json!([]));
        assert_eq!(s.journal.len(), 1);

        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let volume = bodies[0]["volume"].as_f64().unwrap();
        let expect = (40.0 * 30.0 - std::f64::consts::PI * 25.0) * 10.0;
        assert!((volume - expect).abs() < 1e-3 * expect, "{volume}");
        let faces = bodies[0]["faces"].as_array().unwrap();
        assert_eq!(faces.len(), 8);
        let edges = bodies[0]["edges"].as_array().unwrap();
        assert_eq!(edges.len(), 18);
        // Through the thickness: four box corners, the two hole seams.
        assert_eq!(edges.iter().filter(|e| (e["length"].as_f64().unwrap() - 10.0).abs() < 1e-6).count(), 6);

        let t = s.execute("geometry.tessellate", &json!({"body": "Body1", "tolerance": 0.05})).unwrap();
        let per_face: Vec<u64> = t["faces"].as_array().unwrap().iter().map(|f| f["triangles"].as_u64().unwrap()).collect();
        assert!(per_face.iter().all(|&n| n > 0), "{t}");
        assert_eq!(per_face.iter().sum::<u64>(), t["triangles"].as_u64().unwrap());

        // A failed open keeps the design that was open.
        assert!(s.execute("geometry.open", &json!({"path": "/no/such/design.solidcraft"})).is_err());
        assert!(s.design.is_some());
        assert_eq!(s.journal.len(), 1);

        // Re-opening the edited design rebuilds it with the same face names.
        std::fs::write(&path, plate_design(&Plate { width: 60.0, ..Plate::default() }).unwrap()).unwrap();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let wider = s.execute("geometry.bodies", &json!({})).unwrap();
        assert!(wider[0]["volume"].as_f64().unwrap() > volume);
        assert_eq!(wider[0]["faces"], bodies[0]["faces"]);
        let names = |b: &Value| b[0]["edges"].as_array().unwrap().iter().map(|e| e["name"].clone()).collect::<Vec<_>>();
        assert_eq!(names(&wider), names(&bodies));
        assert_eq!(s.journal.len(), 2);
    }

    /// A meshed, defined model of body `body` (a hex block whose bottom face is named `bottom`)
    /// with a set on it and one on its face `face`.
    fn defined_model(body: &str, face: &str) -> Session {
        let mut s = Session::default();
        let cube = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]];
        let script = json!([
            {"command": "mesh.block", "params": {"body": body, "kind": "hex8", "corners": cube, "divisions": [1, 1, 1], "faces": {"1": "bottom"}}},
            {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": body}]}},
            {"command": "model.set.create", "params": {"name": "hole", "geometry": [{"body": body, "face": face}]}},
            {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}},
        ]);
        s.run_script(&script).unwrap();
        s
    }

    /// `opened` changes with another document (a sample, another design file, a new model) and
    /// stays with an edit, a solve or the same design opened again: views frame on that change.
    #[test]
    fn the_open_document_changes_only_when_another_is_opened() {
        let mut s = Session::default();
        let start = s.opened();
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        let cantilever = s.opened();
        assert_ne!(cantilever, start);
        s.execute("solve.run", &json!({})).unwrap();
        s.execute("model.material.create", &json!({"name": "Al", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        assert_eq!(s.opened(), cantilever, "edits and solves keep the document");
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        assert_ne!(s.opened(), cantilever, "a sample opened again is a fresh start");

        let (_dir, path) = sample_design_file(&Plate::default());
        let before = s.opened();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let design = s.opened();
        assert_ne!(design, before, "another design file");
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        assert_eq!(s.opened(), design, "the same design opened again after an edit");
        s.execute("geometry.open", &json!({"path": path, "model": "new"})).unwrap();
        assert_ne!(s.opened(), design, "a new model");
    }

    #[test]
    fn reopening_the_same_design_keeps_the_mesh_and_definitions() {
        let (_dir, path) = sample_design_file(&Plate::default());
        let mut probe = Session::default();
        probe.execute("geometry.open", &json!({"path": path})).unwrap();
        let faces = probe.execute("geometry.bodies", &json!({})).unwrap()[0]["faces"].clone();
        let mut s = defined_model("Body1", faces[0]["name"].as_str().unwrap());
        let before = s.model.clone();
        let r = s.execute("geometry.open", &json!({"path": path})).unwrap();
        assert_eq!((r["model"].as_str(), r["related"].as_bool()), (Some("kept"), Some(true)), "{r}");
        assert_eq!(r["stale"], json!({"mesh_bodies": [], "sets": [], "seeds": [], "references": []}));
        assert_eq!(s.model, before);
    }

    #[test]
    fn a_different_design_starts_a_new_model_or_reports_what_is_stale() {
        let (_dir, path) = sample_design_file(&Plate::default());
        let defined = defined_model("Seal", "no-such-face");
        let before = defined.model.clone();
        assert!(!before.elements().is_empty() && before.sets().len() == 2);

        // Kept (the default): nothing is lost, and the reply names what no longer resolves.
        let mut kept = defined.clone();
        let r = kept.execute("geometry.open", &json!({"path": path})).unwrap();
        assert_eq!((r["model"].as_str(), r["related"].as_bool()), (Some("kept"), Some(false)), "{r}");
        assert_eq!(r["stale"], json!({"mesh_bodies": ["Seal"], "sets": ["all", "hole"], "seeds": [], "references": []}));
        assert_eq!(kept.model, before);

        // New: the mesh and every definition go, the units stay; one undo brings them back.
        let mut fresh = defined.clone();
        fresh.execute("model.units", &json!({"system": "m-kg-s-N-Pa"})).unwrap();
        let units = fresh.model.display_units().cloned();
        let before = fresh.model.clone();
        let r = fresh.execute("geometry.open", &json!({"path": path, "model": "new"})).unwrap();
        assert_eq!((r["model"].as_str(), r["related"].as_bool()), (Some("new"), Some(false)), "{r}");
        assert_eq!(r["stale"], json!({"mesh_bodies": [], "sets": [], "seeds": [], "references": []}));
        assert!(fresh.model.elements().is_empty() && fresh.model.nodes().is_empty());
        assert!(fresh.model.definition_counts().iter().all(|(_, n)| *n == 0));
        assert_eq!(fresh.model.display_units().cloned(), units);
        assert!(fresh.design.is_some());
        fresh.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(fresh.model, before);
        assert!(fresh.design.is_none());

        // A bad choice changes nothing.
        let mut s = defined.clone();
        assert!(s.execute("geometry.open", &json!({"path": path, "model": "discard"})).is_err());
        assert!(s.design.is_none() && s.model == defined.model);
    }

    #[test]
    fn a_new_model_is_untitled_and_undo_brings_the_file_back() {
        let (dir, design) = sample_design_file(&Plate::default());
        let file = dir.0.join("a.fieldcraft").to_string_lossy().into_owned();
        let mut s = defined_model("Seal", "no-such-face");
        s.execute("file.save_as", &json!({"path": file})).unwrap();
        assert!(!s.is_dirty());

        s.execute("geometry.open", &json!({"path": design, "model": "new"})).unwrap();
        assert_eq!(s.path, None, "a new model must not save over the old model's file");
        assert!(s.is_dirty());
        assert!(s.execute("file.save", &json!({})).is_err(), "save asks for a file");

        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.path.as_deref(), Some(file.as_str()));
        assert!(!s.is_dirty(), "undo returns to the saved state of the old file");

        s.execute("edit.redo", &json!({})).unwrap();
        assert_eq!(s.path, None);
        assert!(s.is_dirty());

        // Saving the new model and then editing it: undoing the edit keeps the new file.
        let b = dir.0.join("b.fieldcraft").to_string_lossy().into_owned();
        s.execute("file.save_as", &json!({"path": b})).unwrap();
        s.execute("model.material.create", &json!({"name": "Al", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.path.as_deref(), Some(b.as_str()));
        assert!(!s.is_dirty());
        // Undoing the open after that save returns to the old file; redo to the saved new one.
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.path.as_deref(), Some(file.as_str()));
        s.execute("edit.redo", &json!({})).unwrap();
        assert_eq!(s.path.as_deref(), Some(b.as_str()));
        assert!(!s.is_dirty());
    }

    #[test]
    fn keeping_the_model_keeps_its_file_and_a_sample_has_none() {
        let (dir, design) = sample_design_file(&Plate::default());
        let file = dir.0.join("a.fieldcraft").to_string_lossy().into_owned();
        let mut s = defined_model("Seal", "no-such-face");
        s.execute("file.save_as", &json!({"path": file})).unwrap();
        s.execute("geometry.open", &json!({"path": design})).unwrap();
        assert_eq!(s.path.as_deref(), Some(file.as_str()));
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        assert_eq!(s.path, None, "a sample must not save over the last opened file");
        assert!(s.execute("file.save", &json!({})).is_err());
    }
}
