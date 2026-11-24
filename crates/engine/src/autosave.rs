//! Autosave for recovery after a crash.
//!
//! The desktop app keeps an [`Autosaver`] and calls [`Autosaver::save`] on a timer: while the
//! session has unsaved changes not yet autosaved, the whole model is written to the app's entry in
//! the recovery folder ([`fieldcraft_io::recovery`]); once nothing is unsaved the entry goes. After
//! a crash the next launch lists the entry (`file.recovery.list`) and opens it (`file.recover`).

use std::path::Path;

use fieldcraft_io::recovery::Slot;

use crate::external_view::ExternalView;
use crate::{Result, Session};

pub use fieldcraft_io::recovery::default_dir;

/// One app's autosaver.
#[derive(Debug)]
pub struct Autosaver {
    slot: Slot,
    /// Imported frames and their display can change without a model edit.
    saved: Option<(u64, u64, Option<ExternalView>)>,
}

impl Autosaver {
    /// Start autosaving into `dir` (created if needed).
    pub fn new(dir: &Path) -> Result<Autosaver> {
        Ok(Autosaver { slot: Slot::new(dir)?, saved: None })
    }

    /// The id of this app's entry.
    pub fn id(&self) -> &str {
        self.slot.id()
    }

    /// Autosave if the session has unsaved changes not written yet, or drop the entry when
    /// nothing is unsaved. Returns whether it wrote.
    pub fn save(&mut self, s: &Session) -> Result<bool> {
        if !s.is_dirty() {
            self.slot.clear();
            self.saved = None;
            return Ok(false);
        }
        let current = (s.revision(), s.external_generation(), s.external_view.clone());
        if self.saved.as_ref() == Some(&current) {
            return Ok(false);
        }
        // Geometry paths stay absolute: the entry is not next to the model's own file.
        let bytes = s.recovery_bytes()?;
        self.slot.write(&bytes, s.path.as_deref(), s.model.nodes().len(), s.model.elements().len())?;
        self.saved = Some(current);
        crate::crash::autosaved(s, &self.slot.file());
        Ok(true)
    }

    /// The app is closing: a last autosave, then the lock is released. Unsaved work stays in the
    /// recovery folder for the next launch to offer.
    pub fn close(mut self, s: &Session) -> Result<()> {
        let saved = self.save(s);
        self.slot.close();
        saved.map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    #[test]
    fn autosaves_unsaved_work_and_recovers_it_after_a_crash() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-autosave-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let mut s = Session::default();
        let mut a = Autosaver::new(&dir).unwrap();
        assert!(!a.save(&s).unwrap(), "nothing unsaved");
        s.execute("mesh.node.add", &json!({"label": 1, "xyz": [1, 2, 3]})).unwrap();
        assert!(a.save(&s).unwrap());
        assert!(!a.save(&s).unwrap(), "unchanged since the last autosave");
        let id = a.id().to_string();
        let listed = s.execute("file.recovery.list", &json!({"dir": dir})).unwrap();
        assert_eq!(listed, json!([]), "a running app's autosave is not offered");

        // The app dies.
        drop(a);
        let listed = s.execute("file.recovery.list", &json!({"dir": dir})).unwrap();
        assert_eq!(listed[0]["id"], id);
        assert_eq!(listed[0]["nodes"], 1);
        assert!(listed[0]["age"].as_u64().is_some_and(|a| a < 60), "{listed}");
        let mut next = Session::default();
        let r = next.execute("file.recover", &json!({"id": id, "dir": dir})).unwrap();
        assert_eq!(r["nodes"], 1);
        assert_eq!(next.model, s.model);
        assert_eq!(next.journal, s.journal);
        assert!(next.is_dirty(), "recovered work is unsaved");
        next.execute("file.recovery.discard", &json!({"all": true, "dir": dir})).unwrap();
        assert_eq!(next.execute("file.recovery.list", &json!({"dir": dir})).unwrap(), json!([]));

        // Undoing back to the saved state drops the entry; a clean close leaves nothing.
        let mut a = Autosaver::new(&dir).unwrap();
        let mut s = Session::default();
        s.execute("mesh.node.add", &json!({"label": 1, "xyz": [1, 2, 3]})).unwrap();
        assert!(a.save(&s).unwrap());
        s.execute("edit.undo", &json!({})).unwrap();
        assert!(!a.save(&s).unwrap());
        a.close(&s).unwrap();
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
