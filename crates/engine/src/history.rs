//! Undo and redo on the command journal.
//!
//! Every journalled command leaves an undo step: the state before it (model, design, surfaces
//! and their sources). `edit.undo` puts that state back and moves the command from the end of the
//! journal onto the redo stack with the state it led to; `edit.redo` does the opposite. So the
//! journal always replays to the current model, and the undo steps are the journal's last
//! entries, newest last. A new command clears the redo stack.
//!
//! A step also carries the session's [`Document`] (its file and the revision saved there) when
//! its command changed which document the session is, as `geometry.open` with `"model": "new"`
//! does: undo makes the previous file current again, redo the new untitled model. Other steps
//! leave the document alone, so saving between edits and then undoing them keeps the file.
//!
//! History is bounded by step count ([`MAX_UNDO_STEPS`]) and by an estimate of the memory the
//! kept states hold ([`DEFAULT_UNDO_MEMORY`], [`Session::set_undo_memory`]): when either is
//! exceeded, the oldest undo steps go first, then the redo steps farthest from the current state.
//! Designs and surfaces are shared ([`Arc`]), so a step only pays for those it alone keeps alive.

use std::collections::VecDeque;
use std::sync::Arc;

use fieldcraft_geometry::Design;
use fieldcraft_geometry::tagged_mesh::TaggedSurface;
use fieldcraft_model::Model;
use serde_json::{Value, json};

use crate::{EngineError, Result, Session};

/// Most undo steps kept.
pub const MAX_UNDO_STEPS: usize = 500;

/// Memory the undo and redo steps may hold by default (estimated).
pub const DEFAULT_UNDO_MEMORY: usize = 1 << 30;

/// Estimated bytes of one node (label, coordinates, map entry overhead).
const NODE_BYTES: usize = 64;
/// Estimated bytes of one element without its node list.
const ELEMENT_BYTES: usize = 72;
const TAGGED_ELEMENT_BYTES: usize = 160;
const TAGGED_NODE_BYTES: usize = 32;
const DEFINITION_BYTES: usize = 512;
/// A rebuilt design's B-rep is not measurable through its API; each design a step alone keeps
/// alive counts as this much.
const DESIGN_BYTES: usize = 16 << 20;

/// What an undo step restores.
#[derive(Debug, Clone, Default)]
pub(crate) struct State {
    model: Model,
    design: Option<Arc<Design>>,
    surfaces: Vec<Arc<TaggedSurface>>,
    sources: fieldcraft_io::Geometry,
    revision: u64,
    /// The document before the step, when the step changed it.
    pub(crate) document: Option<Document>,
}

/// Which document the session is: the `.fieldcraft` file it belongs to and the revision saved
/// there ([`Session::is_dirty`]).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Document {
    pub(crate) path: Option<String>,
    pub(crate) saved_revision: Option<u64>,
    pub(crate) migrated_from: Option<u32>,
}

#[derive(Debug, Clone)]
struct UndoStep {
    /// The state before the journal's matching entry ran.
    state: State,
    bytes: usize,
}

#[derive(Debug, Clone)]
struct RedoStep {
    /// The state the entry led to.
    state: State,
    entry: (String, Value),
    bytes: usize,
}

#[derive(Debug, Clone)]
pub(crate) struct History {
    undo: VecDeque<UndoStep>,
    redo: Vec<RedoStep>,
    /// Estimated bytes held by both stacks.
    bytes: usize,
    limit: usize,
}

impl Default for History {
    fn default() -> Self {
        History { undo: VecDeque::new(), redo: Vec::new(), bytes: 0, limit: DEFAULT_UNDO_MEMORY }
    }
}

impl History {
    fn push_undo(&mut self, step: UndoStep) {
        self.bytes += step.bytes;
        self.undo.push_back(step);
        self.trim();
    }

    fn pop_undo(&mut self) -> Option<UndoStep> {
        let step = self.undo.pop_back()?;
        self.bytes -= step.bytes;
        Some(step)
    }

    fn pop_redo(&mut self) -> Option<RedoStep> {
        let step = self.redo.pop()?;
        self.bytes -= step.bytes;
        Some(step)
    }

    fn push_redo(&mut self, step: RedoStep) {
        self.bytes += step.bytes;
        self.redo.push(step);
        self.trim();
    }

    fn clear_redo(&mut self) {
        self.bytes -= self.redo.iter().map(|r| r.bytes).sum::<usize>();
        self.redo.clear();
    }

    /// Drop the oldest undo steps, then the farthest redo steps, until within the limits.
    fn trim(&mut self) {
        while self.undo.len() > MAX_UNDO_STEPS || (self.bytes > self.limit && !self.undo.is_empty()) {
            if let Some(old) = self.undo.pop_front() {
                self.bytes -= old.bytes;
            }
        }
        let mut cut = 0;
        while self.bytes > self.limit && cut < self.redo.len() {
            self.bytes -= self.redo.get(cut).map_or(0, |r| r.bytes);
            cut += 1;
        }
        self.redo.drain(..cut);
    }
}

fn model_bytes(m: &Model) -> usize {
    let elements: usize = m.elements().values().map(|e| ELEMENT_BYTES + 4 * e.nodes.len()).sum();
    let c = m.classification();
    let definitions: usize = m.definition_counts().iter().map(|(_, n)| n).sum();
    m.nodes().len() * NODE_BYTES
        + elements
        + c.elements().len() * TAGGED_ELEMENT_BYTES
        + c.nodes().len() * TAGGED_NODE_BYTES
        + definitions * DEFINITION_BYTES
}

fn surface_bytes(s: &TaggedSurface) -> usize {
    s.mesh.positions.len() * 24 + s.mesh.triangles.len() * 16
}

impl Session {
    pub(crate) fn snapshot(&self) -> State {
        State {
            model: self.model.clone(),
            design: self.design.clone(),
            surfaces: self.surfaces.clone(),
            sources: self.sources.clone(),
            revision: self.revision,
            document: None,
        }
    }

    pub(crate) fn document(&self) -> Document {
        Document { path: self.path.clone(), saved_revision: self.saved_revision, migrated_from: self.migrated_from }
    }

    /// Make `document` current; returns the one it replaces.
    pub(crate) fn set_document(&mut self, document: Document) -> Document {
        Document {
            path: std::mem::replace(&mut self.path, document.path),
            saved_revision: std::mem::replace(&mut self.saved_revision, document.saved_revision),
            migrated_from: std::mem::replace(&mut self.migrated_from, document.migrated_from),
        }
    }

    /// Make `state` current; returns the state it replaces.
    pub(crate) fn restore(&mut self, state: State) -> State {
        State {
            model: std::mem::replace(&mut self.model, state.model),
            design: std::mem::replace(&mut self.design, state.design),
            surfaces: std::mem::replace(&mut self.surfaces, state.surfaces),
            sources: std::mem::replace(&mut self.sources, state.sources),
            revision: std::mem::replace(&mut self.revision, state.revision),
            document: state.document.map(|d| self.set_document(d)),
        }
    }

    /// Whether the session's state differs from `state`.
    pub(crate) fn differs_from(&self, state: &State) -> bool {
        let shared = |a: &[Arc<TaggedSurface>], b: &[Arc<TaggedSurface>]| a.len() == b.len() && a.iter().zip(b).all(|(x, y)| Arc::ptr_eq(x, y));
        self.model != state.model
            || self.design.as_ref().map(Arc::as_ptr) != state.design.as_ref().map(Arc::as_ptr)
            || !shared(&self.surfaces, &state.surfaces)
            || self.sources != state.sources
    }

    /// Estimated memory `state` holds beyond what the current state shares with it.
    fn kept_bytes(&self, state: &State) -> usize {
        let design = match (&state.design, &self.design) {
            (Some(d), Some(current)) if Arc::ptr_eq(d, current) => 0,
            (Some(_), _) => DESIGN_BYTES,
            (None, _) => 0,
        };
        let surfaces: usize = state.surfaces.iter().filter(|s| !self.surfaces.iter().any(|c| Arc::ptr_eq(c, s))).map(|s| surface_bytes(s)).sum();
        model_bytes(&state.model) + design + surfaces
    }

    /// A journalled command ran from `before`: keep it as an undo step.
    pub(crate) fn record_undo(&mut self, before: State) {
        let bytes = self.kept_bytes(&before);
        self.history.clear_redo();
        self.history.push_undo(UndoStep { state: before, bytes });
    }

    /// Forget every undo and redo step (a new model was opened).
    pub(crate) fn clear_history(&mut self) {
        let limit = self.history.limit;
        self.history = History { limit, ..History::default() };
    }

    /// Bound the estimated memory of the undo history (bytes); older steps beyond it go.
    pub fn set_undo_memory(&mut self, bytes: usize) {
        self.history.limit = bytes;
        self.history.trim();
    }

    /// The number of steps that can be undone and redone.
    pub fn undo_depth(&self) -> (usize, usize) {
        (self.history.undo.len(), self.history.redo.len())
    }

    /// The undone commands redo would run, the next one first.
    pub fn redo_entries(&self) -> impl Iterator<Item = &(String, Value)> {
        self.history.redo.iter().rev().map(|r| &r.entry)
    }

    fn history_report(&self, done: &str, command: &str) -> Value {
        let (undo, redo) = self.undo_depth();
        json!({done: command, "undo": undo, "redo": redo})
    }

    /// Undo the last journalled command.
    pub fn undo(&mut self) -> Result<Value> {
        if self.history.undo.is_empty() || self.journal.is_empty() {
            return Err(EngineError::Params("nothing to undo".into()));
        }
        let (Some(step), Some(entry)) = (self.history.pop_undo(), self.journal.pop()) else {
            return Err(EngineError::Params("nothing to undo".into()));
        };
        let after = self.restore(step.state);
        let bytes = self.kept_bytes(&after);
        let command = entry.0.clone();
        self.history.push_redo(RedoStep { state: after, entry, bytes });
        Ok(self.history_report("undone", &command))
    }

    /// Redo the last undone command.
    pub fn redo(&mut self) -> Result<Value> {
        let step = self.history.pop_redo().ok_or_else(|| EngineError::Params("nothing to redo".into()))?;
        let before = self.restore(step.state);
        let bytes = self.kept_bytes(&before);
        let command = step.entry.0.clone();
        self.journal.push(step.entry);
        self.history.push_undo(UndoStep { state: before, bytes });
        Ok(self.history_report("redone", &command))
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    fn add_node(s: &mut Session, label: u32) {
        s.execute("mesh.node.add", &json!({"label": label, "xyz": [label, 0, 0]})).unwrap();
    }

    #[test]
    fn undo_and_redo_walk_the_journal() {
        let mut s = Session::default();
        assert!(s.execute("edit.undo", &json!({})).is_err());
        assert!(s.execute("edit.redo", &json!({})).is_err());
        for label in 1..=3 {
            add_node(&mut s, label);
        }
        s.execute("model.info", &json!({})).unwrap();
        assert_eq!(s.undo_depth(), (3, 0));
        let three = (s.model.clone(), s.revision());

        let r = s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(r, json!({"undone": "mesh.node.add", "undo": 2, "redo": 1}));
        assert_eq!(s.model.nodes().len(), 2);
        assert_eq!(s.journal.len(), 2);
        s.execute("edit.undo", &json!({})).unwrap();
        s.execute("edit.undo", &json!({})).unwrap();
        assert!(s.model.nodes().is_empty() && s.journal.is_empty());
        assert!(!s.is_dirty(), "back to the start");
        assert!(s.execute("edit.undo", &json!({})).is_err());

        for _ in 0..3 {
            s.execute("edit.redo", &json!({})).unwrap();
        }
        assert_eq!((s.model.clone(), s.revision()), three);
        assert_eq!(s.journal.len(), 3);
        assert!(s.execute("edit.redo", &json!({})).is_err());

        // A new command after an undo drops what could be redone.
        s.execute("edit.undo", &json!({})).unwrap();
        add_node(&mut s, 9);
        assert_eq!(s.undo_depth(), (3, 0));
        assert_eq!(s.journal.last().unwrap().1["label"], 9);
        // The journal replays to the current model.
        let mut replay = Session::default();
        let script: Vec<Value> = s.journal.iter().map(|(c, p)| json!({"command": c, "params": p})).collect();
        replay.run_script(&Value::Array(script)).unwrap();
        assert_eq!(replay.model, s.model);
    }

    #[test]
    fn failed_commands_and_queries_are_not_undo_steps() {
        let mut s = Session::default();
        add_node(&mut s, 1);
        assert!(s.execute("mesh.node.add", &json!({"label": 0, "xyz": [0, 0, 0]})).is_err());
        s.execute("model.check", &json!({})).unwrap();
        assert_eq!(s.undo_depth(), (1, 0));
    }

    #[test]
    fn history_is_bounded_by_steps_and_memory() {
        let mut s = Session::default();
        for label in 1..=(MAX_UNDO_STEPS as u32 + 20) {
            add_node(&mut s, label);
        }
        assert_eq!(s.undo_depth().0, MAX_UNDO_STEPS);
        // The oldest steps go; undo still works back to there and stops.
        for _ in 0..MAX_UNDO_STEPS {
            s.execute("edit.undo", &json!({})).unwrap();
        }
        assert_eq!(s.model.nodes().len(), 20);
        assert!(s.execute("edit.undo", &json!({})).is_err());
        assert_eq!(s.journal.len(), 20, "the journal keeps what can no longer be undone");

        // A memory cap keeps only the steps that fit.
        let mut s = Session::default();
        s.set_undo_memory(40 * NODE_BYTES);
        for label in 1..=20 {
            add_node(&mut s, label);
        }
        // The steps before the 20th and 19th node hold 19 and 18 nodes.
        assert_eq!(s.undo_depth(), (2, 0));
        assert!(s.history.bytes <= 40 * NODE_BYTES);
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.undo_depth(), (1, 1), "the redo step fits too");
        s.set_undo_memory(0);
        assert_eq!(s.undo_depth(), (0, 0));
        assert_eq!(s.history.bytes, 0);
    }
}
