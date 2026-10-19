//! Files: open and save `.fieldcraft` models, and recover autosaved ones.
//!
//! A saved model refers to its design and imported meshes by path and content hash; opening it
//! opens them again. Paths inside the model file's folder are stored relative to it, so a folder
//! of a model and its design can be moved together. A geometry file that is missing, unreadable
//! or changed since the save does not stop the model from opening: it is reported in the result's
//! `warnings` (sets follow geometry names, so a remesh brings a changed design's mesh up to date).

use std::path::{Component, Path, PathBuf};
use std::sync::Arc;

use fieldcraft_geometry::tagged_mesh;
use fieldcraft_geometry::{Design, MAX_DESIGN_BYTES, open_design};
use fieldcraft_io::{Geometry, ModelFile, SourceFile};
use serde_json::{Value, json};

use fieldcraft_io::roots::{Access, Roots};

use super::{file_param, req};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, next_revision};

pub const NEW: CommandSpec = CommandSpec {
    id: "file.new",
    label: "New Model",
    params: r#"{} — starts an empty model in place of the current one: no design, mesh, definitions, results or journal; unsaved changes are discarded (the window asks first) and the undo history starts empty"#,
    form: Some(&[]),
    kind: CommandKind::Session,
    requires: &[],
    run: new_model,
};

pub const OPEN: CommandSpec = CommandSpec {
    id: "file.open",
    label: "Open Model",
    params: r#"{"path": "bracket.fieldcraft"} — opens the model, its mesh and its journal, and opens its design and imported meshes again from their files (one that is missing or changed is listed in `warnings`; the model still opens); unsaved changes of the current model are discarded and the undo history starts empty. A file of an older format is migrated as it is read (`migrated_from` in the reply) and nothing is written; a later format is refused"#,
    form: Some(super::forms::FILE_OPEN),
    kind: CommandKind::Session,
    requires: &[],
    run: open,
};

pub const SAVE: CommandSpec = CommandSpec {
    id: "file.save",
    label: "Save Model",
    params: r#"{"path"?: "bracket.fieldcraft"} — saves to `path`, else to the file the model was opened from or last saved to; written atomically, the previous version kept as `<file>.bak`. A file there that is not a model of the current format (one this model was migrated from, a later or damaged one, any other file) is first copied to `<stem>.v<N>.fieldcraft` (`<stem>.backup.<ext>` when it is not a model), never replaced later (`backup` in the reply)"#,
    form: Some(super::forms::FILE_SAVE),
    kind: CommandKind::Session,
    requires: &[],
    run: save,
};

pub const SAVE_AS: CommandSpec = CommandSpec {
    id: "file.save_as",
    label: "Save Model As",
    params: r#"{"path": "bracket.fieldcraft"} — saves to a new file, which becomes the model's file (a file there that is not a model of the current format is first copied, as file.save does)"#,
    form: Some(super::forms::FILE_SAVE_AS),
    kind: CommandKind::Session,
    requires: &[],
    run: save_as,
};

pub const RECOVERY_LIST: CommandSpec = CommandSpec {
    id: "file.recovery.list",
    label: "Recoverable Models",
    params: r#"{"dir"?: "…"} — models autosaved by apps that crashed or closed with unsaved changes, newest first, each with its id, path, saved_at, age (seconds) and node and element counts (default folder: FIELDCRAFT_RECOVERY_DIR, a temporary folder in runs driven by automation, else the per-user data folder)"#,
    form: Some(super::forms::RECOVERY_LIST),
    kind: CommandKind::Session,
    requires: &[],
    run: recovery_list,
};

pub const RECOVER: CommandSpec = CommandSpec {
    id: "file.recover",
    label: "Recover Model",
    params: r#"{"id": "…", "dir"?: "…"} — opens an autosaved model (an id from file.recovery.list) with its original file path, as unsaved changes; the autosave stays until file.recovery.discard"#,
    form: Some(super::forms::RECOVER),
    kind: CommandKind::Session,
    requires: &[],
    run: recover,
};

pub const RECOVERY_DISCARD: CommandSpec = CommandSpec {
    id: "file.recovery.discard",
    label: "Discard Recovered Model",
    params: r#"{"id": "…"} or {"all": true}; "dir"?: "…" — deletes autosaves of apps that are no longer running"#,
    form: Some(super::forms::RECOVERY_DISCARD),
    kind: CommandKind::Session,
    requires: &[],
    run: recovery_discard,
};

pub const SAMPLE: CommandSpec = CommandSpec {
    id: "file.sample",
    label: "Open Sample",
    params: r#"{"name": "cantilever"} — opens a sample model (fieldcraft_engine::sample; the list is in docs/examples.md) in place of the current one, unsaved changes discarded; every sample but `parts` is ready to solve (solve.run)"#,
    form: Some(super::forms::SAMPLE),
    kind: CommandKind::Session,
    requires: &[],
    run: open_sample,
};

fn new_model(s: &mut Session, _: &Value) -> Result<Value> {
    let mut empty = Session::default();
    // Nothing to save until a command changes it.
    empty.saved_revision = Some(empty.revision);
    empty.saved_view = empty.view.clone();
    s.replace_with(empty);
    Ok(summary(s, &[]))
}

fn open_sample(s: &mut Session, p: &Value) -> Result<Value> {
    let name: &str = req(p, "name")?;
    let names = || crate::sample::SAMPLES.iter().map(|x| x.name).collect::<Vec<_>>().join(", ");
    let sample = crate::sample::find(name).ok_or_else(|| EngineError::Params(format!("no sample `{name}` (samples: {})", names())))?;
    let mut opened = Session::default();
    if let Some(recipe) = sample.design {
        opened.design = Some(Arc::new(open_design(&fieldcraft_geometry::script_design(recipe)?)?));
    }
    if let Some(prepare) = sample.prepare {
        opened.run_script(&prepare())?;
    }
    let bodies = opened.design.as_ref().map_or(&[][..], |d| d.bodies.as_slice());
    let script = (sample.script)(bodies);
    opened.run_script(&script)?;
    // A sample just opened has nothing to save: it opens again the same way.
    opened.saved_revision = Some(opened.revision);
    opened.saved_view = opened.view.clone();
    opened.sample = Some(sample.label);
    s.replace_with(opened);
    let mut reply = summary(s, &[]);
    reply["sample"] = json!(sample.name);
    Ok(reply)
}

/// A file of at most `max` bytes, through FieldCraft's file edge ([`fieldcraft_io::fs`]: the disk,
/// or the page's files in the browser).
pub(crate) fn read_file(path: &str, max: u64) -> Result<Vec<u8>> {
    Ok(fieldcraft_io::fs::read_capped(Path::new(path), max)?)
}

/// `path` made absolute against the working folder (as given where there is none).
pub(crate) fn absolute(path: &str) -> String {
    std::path::absolute(path).map(|p| p.to_string_lossy().into_owned()).unwrap_or_else(|_| path.to_string())
}

/// How a geometry path is stored in a file saved in `base`: relative (with `/`) when inside it.
fn stored_path(path: &str, base: Option<&Path>) -> String {
    let relative = base.and_then(|b| Path::new(path).strip_prefix(b).ok()).and_then(|rel| {
        let parts: Option<Vec<&str>> = rel.components().map(|c| if let Component::Normal(s) = c { s.to_str() } else { None }).collect();
        parts.filter(|p| !p.is_empty()).map(|p| p.join("/"))
    });
    relative.unwrap_or_else(|| path.to_string())
}

/// A stored geometry path resolved against the folder of the file that holds it.
fn resolved_path(path: &str, base: Option<&Path>) -> String {
    match base {
        Some(b) if !Path::new(path).is_absolute() => b.join(path).to_string_lossy().into_owned(),
        _ => path.to_string(),
    }
}

/// Open one geometry file again: what it gives, and the source as now read (or as stored, with
/// its path resolved, when it could not be). A file outside the read root of `roots` is not
/// opened.
fn reopen<T>(
    what: &str,
    stored: &SourceFile,
    base: Option<&Path>,
    roots: &Roots,
    max: u64,
    open: impl FnOnce(&[u8]) -> fieldcraft_geometry::Result<T>,
    warnings: &mut Vec<String>,
) -> (Option<T>, SourceFile) {
    let path = resolved_path(stored.path(), base);
    let kept = stored.with_path(path.clone()).unwrap_or_else(|_| stored.clone());
    let checked = match roots.check(Path::new(&path), Access::Read) {
        Ok(p) => p.to_string_lossy().into_owned(),
        Err(e) => {
            warnings.push(format!("{what} `{path}` was not opened ({e}); the model opens without it"));
            return (None, kept);
        }
    };
    let bytes = match read_file(&checked, max) {
        Ok(b) => b,
        Err(e) => {
            warnings.push(format!("{what} `{path}` was not opened ({e}); the model opens without it"));
            return (None, kept);
        }
    };
    let Ok(read) = SourceFile::of(path.clone(), &bytes) else { return (None, kept) };
    if read.sha256() != stored.sha256() {
        warnings
            .push(format!("{what} `{path}` changed since the model was saved; sets follow geometry names, so remeshing brings the mesh up to date"));
    }
    match open(&bytes) {
        Ok(v) => (Some(v), read),
        Err(e) => {
            warnings.push(format!("{what} `{path}` could not be opened ({e}); the model opens without it"));
            (None, kept)
        }
    }
}

impl Session {
    /// The bytes of the `.fieldcraft` file holding this session, with geometry paths inside
    /// `base` (the folder the file goes to) stored relative to it.
    pub fn file_bytes(&self, base: Option<&Path>) -> Result<Vec<u8>> {
        let stored = |s: &SourceFile| s.with_path(stored_path(s.path(), base));
        let geometry = Geometry {
            design: self.sources.design.as_ref().map(stored).transpose()?,
            meshes: self.sources.meshes.iter().map(stored).collect::<fieldcraft_io::Result<_>>()?,
        };
        Ok(fieldcraft_io::write(&self.model, &geometry, &self.journal, &self.view)?)
    }

    /// A session holding what a file holds, its geometry opened again from the paths it names
    /// (relative ones against `base`) where `roots` lets it read them; problems with the
    /// geometry come back as warnings. The session's own file access is unconfined.
    pub fn from_file(file: ModelFile, base: Option<&Path>, roots: &Roots) -> (Session, Vec<String>) {
        let mut warnings = Vec::new();
        let mut sources = Geometry::default();
        let mut design: Option<Arc<Design>> = None;
        if let Some(stored) = &file.geometry.design {
            let (opened, source) = reopen("design", stored, base, roots, MAX_DESIGN_BYTES as u64, open_design, &mut warnings);
            design = opened.and_then(|d| match super::datum::partitioned(&d, &file.model) {
                Ok(cut) => Some(Arc::new(cut)),
                Err(e) => {
                    warnings.push(format!("the design's partitions could not be applied ({e}); the model opens without the design"));
                    None
                }
            });
            if let Some(d) = &design {
                warnings.extend(
                    d.errors
                        .iter()
                        .filter(|e| e.feature.starts_with("partition ") || e.feature.starts_with("merge "))
                        .map(|e| format!("{}: {}", e.feature, e.message)),
                );
            }
            sources.design = Some(source);
        }
        let mut surfaces = Vec::new();
        for stored in &file.geometry.meshes {
            let (opened, source) = reopen("mesh", stored, base, roots, tagged_mesh::MAX_FILE_BYTES as u64, tagged_mesh::read_3mf, &mut warnings);
            surfaces.extend(opened.into_iter().flatten().map(Arc::new));
            sources.meshes.push(source);
        }
        let mut session = Session {
            model: file.model,
            design,
            surfaces,
            sources,
            journal: file.journal,
            revision: next_revision(),
            saved_view: file.view.clone(),
            view: file.view,
            ..Session::default()
        };
        let (_, failed) = super::stitch::restitch(&mut session);
        warnings.extend(failed.into_iter().map(|(n, why)| format!("stitch `{n}`: {why}")));
        (session, warnings)
    }

    /// Replace this session with `opened`, keeping its undo memory limit, its work caps and its
    /// file access (and where its commands come from); the history starts empty.
    fn replace_with(&mut self, mut opened: Session) {
        opened.history = std::mem::take(&mut self.history);
        opened.limits = self.limits;
        opened.file_access = std::mem::take(&mut self.file_access);
        opened.crash_journal = self.crash_journal;
        opened.origin = self.origin;
        opened.clear_history();
        opened.mark_opened();
        *self = opened;
    }
}

fn summary(s: &Session, warnings: &[String]) -> Value {
    json!({
        "path": s.path,
        "nodes": s.model.nodes().len(),
        "elements": s.model.elements().len(),
        "journal": s.journal.len(),
        "design": s.design.as_ref().map(|d| d.name.clone()),
        "surfaces": s.surfaces.len(),
        "migrated_from": s.migrated_from,
        "warnings": warnings,
    })
}

/// Replace the session with a `.fieldcraft` file's content; `path` is its own file (geometry
/// paths are relative to its folder), `saved` whether that file holds this content (else the
/// content is an autosave of it). Nothing is written.
fn load(s: &mut Session, bytes: &[u8], path: Option<String>, saved: bool) -> Result<Value> {
    let (file, format) = fieldcraft_io::read_versioned(bytes)?;
    let (mut opened, warnings) = Session::from_file(file, path.as_deref().and_then(|p| Path::new(p).parent()), &s.file_access);
    // The format of the file at `path`: the one just read, or (an autosave) the one on disk,
    // read only where the session may read. Only an older one is "migrated from"; a save over a
    // later one keeps a copy all the same (fs::save_model decides from the file itself).
    let on_disk = match &path {
        Some(_) if saved => Some(format),
        Some(p) => s.checked_path(p, Access::Read).ok().and_then(|p| fieldcraft_io::fs::format_of_file(Path::new(&p))),
        None => None,
    };
    opened.migrated_from = on_disk.filter(|&f| f < fieldcraft_io::FORMAT_VERSION);
    opened.path = path;
    if saved {
        opened.saved_revision = Some(opened.revision);
    }
    s.replace_with(opened);
    Ok(summary(s, &warnings))
}

fn open(s: &mut Session, p: &Value) -> Result<Value> {
    let path = absolute(&file_param(s, p, "path", Access::Read)?);
    let bytes = read_file(&path, fieldcraft_io::MAX_FILE_BYTES)?;
    load(s, &bytes, Some(path.clone()), true).map_err(|e| EngineError::Params(format!("{path}: {e}")))
}

/// Save to `path` (as given by the caller).
fn save_to(s: &mut Session, path: &str) -> Result<Value> {
    let target = absolute(path);
    if std::fs::metadata(&target).is_ok_and(|m| !m.is_file()) {
        return Err(EngineError::Params(format!("`{target}` exists and is not a regular file")));
    }
    let bytes = s.file_bytes(Path::new(&target).parent())?;
    let kept = fieldcraft_io::fs::save_model(Path::new(&target), &bytes)?;
    s.mark_saved(target.clone());
    let mut reply = json!({"path": target, "bytes": bytes.len()});
    if let Some(kept) = kept {
        let name = |p: &Path| p.file_name().map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into_owned());
        reply["backup"] = json!(kept.to_string_lossy());
        reply["summary"] = json!(format!(
            "saved {} in format {}; the file it replaced is kept as {}",
            name(Path::new(&target)),
            fieldcraft_io::FORMAT_VERSION,
            name(&kept)
        ));
    }
    Ok(reply)
}

fn save(s: &mut Session, p: &Value) -> Result<Value> {
    let path = match p.get("path") {
        Some(_) => file_param(s, p, "path", Access::Write)?,
        None => {
            let own = s.path.as_deref().ok_or_else(|| EngineError::Params("the model has no file yet: give `path` (or use file.save_as)".into()))?;
            s.checked_path(own, Access::Write)?
        }
    };
    save_to(s, &path)
}

fn save_as(s: &mut Session, p: &Value) -> Result<Value> {
    let path = file_param(s, p, "path", Access::Write)?;
    save_to(s, &path)
}

/// The recovery folder: `dir` (checked for `access`), else the default one (the app's own folder,
/// which no client names, so it is not confined).
fn recovery_dir(s: &Session, p: &Value, access: Access) -> Result<PathBuf> {
    match p.get("dir") {
        Some(_) => Ok(PathBuf::from(file_param(s, p, "dir", access)?)),
        None => {
            crate::autosave::default_dir().ok_or_else(|| EngineError::Params("no recovery folder: set FIELDCRAFT_RECOVERY_DIR or give `dir`".into()))
        }
    }
}

fn recovery_list(s: &mut Session, p: &Value) -> Result<Value> {
    let mut listed = Vec::new();
    for e in fieldcraft_io::recovery::orphans(&recovery_dir(s, p, Access::Read)?) {
        let mut v = serde_json::to_value(&e).map_err(|e| EngineError::Params(e.to_string()))?;
        if let Some(o) = v.as_object_mut() {
            o.insert("age".into(), json!(e.age()));
        }
        listed.push(v);
    }
    Ok(Value::Array(listed))
}

fn recover(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = recovery_dir(s, p, Access::Read)?;
    let (bytes, entry) = fieldcraft_io::recovery::load(&dir, req::<&str>(p, "id")?)?;
    // Recovered work is unsaved until it is saved. Autosaves keep geometry paths absolute.
    load(s, &bytes, entry.path, false)
}

fn recovery_discard(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = recovery_dir(s, p, Access::Write)?;
    if p.get("all").and_then(Value::as_bool) == Some(true) {
        let ids: Vec<String> = fieldcraft_io::recovery::orphans(&dir).into_iter().map(|e| e.id).collect();
        for id in &ids {
            fieldcraft_io::recovery::remove(&dir, id)?;
        }
        return Ok(json!({"discarded": ids}));
    }
    let id = req::<&str>(p, "id")?;
    fieldcraft_io::recovery::remove(&dir, id)?;
    Ok(json!({"discarded": [id]}))
}

#[cfg(test)]
mod tests {
    use fieldcraft_geometry::sample::{Plate, plate_design};
    use serde_json::json;

    use super::*;

    fn dir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("fieldcraft-file-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn defined(s: &mut Session) {
        for (i, c) in [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]].iter().enumerate() {
            s.execute("mesh.node.add", &json!({"label": i + 1, "xyz": c})).unwrap();
        }
        s.execute("mesh.element.add", &json!({"label": 1, "kind": "tet4", "nodes": [1, 2, 3, 4]})).unwrap();
        s.execute("model.material.create", &json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}})).unwrap();
    }

    /// `file.new` empties the model, its journal and its undo history, and has nothing to save.
    #[test]
    fn new_starts_an_empty_model() {
        let mut s = Session::default();
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        s.execute("mesh.node.add", &json!({"label": 10_000, "xyz": [0, 0, 0]})).unwrap();
        let r = s.execute("file.new", &json!({})).unwrap();
        assert_eq!((r["nodes"].as_u64(), r["elements"].as_u64(), r["journal"].as_u64()), (Some(0), Some(0), Some(0)), "{r}");
        assert!(s.model.materials().is_empty() && s.model.steps().is_empty() && s.design.is_none() && s.sample.is_none());
        assert!(!s.is_dirty(), "an empty model has nothing to save");
        assert!(s.execute("edit.undo", &json!({})).is_err(), "the undo history starts empty");
    }

    #[test]
    fn save_and_open_keep_the_model_journal_and_design() {
        let d = dir("roundtrip");
        let design = d.join("plate.solidcraft");
        std::fs::write(&design, plate_design(&Plate::default()).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": design})).unwrap();
        defined(&mut s);
        assert!(s.execute("file.save", &json!({})).is_err(), "no file yet, and no dialog");
        let path = d.join("bracket.fieldcraft");
        s.execute("file.save_as", &json!({"path": path})).unwrap();
        assert!(!s.is_dirty());
        let bytes = std::fs::read(&path).unwrap();
        let file = fieldcraft_io::read(&bytes).unwrap();
        assert_eq!(file.geometry.design.as_ref().unwrap().path(), "plate.solidcraft", "stored next to the model, relative");

        // Saving again keeps the previous version and gives the same bytes.
        s.execute("file.save", &json!({})).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        assert_eq!(std::fs::read(fieldcraft_io::fs::backup_path(&path)).unwrap(), bytes);

        let mut t = Session::default();
        t.execute("mesh.node.add", &json!({"label": 9, "xyz": [0, 0, 0]})).unwrap();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert_eq!(r["warnings"], json!([]));
        assert_eq!(r["design"], json!(s.design.as_ref().unwrap().name));
        assert_eq!(t.model, s.model);
        assert_eq!(t.journal, s.journal);
        assert_eq!(t.undo_depth(), (0, 0), "a new model's history starts empty");
        assert!(!t.is_dirty());
        assert!(t.design.is_some());

        // The design changed since: the model opens, with a warning.
        std::fs::write(&design, plate_design(&Plate { width: 60.0, ..Plate::default() }).unwrap()).unwrap();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert!(r["warnings"][0].as_str().unwrap().contains("changed since"), "{r}");
        assert!(t.design.is_some());

        // The design is gone: the model opens without it and keeps the reference.
        std::fs::remove_file(&design).unwrap();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert!(r["warnings"][0].as_str().unwrap().contains("not opened"), "{r}");
        assert!(t.design.is_none());
        assert_eq!(t.model, s.model);
        assert_eq!(t.sources.design.as_ref().unwrap().path(), design.to_string_lossy());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_confined_session_reads_and_writes_only_inside_its_roots() {
        let d = dir("roots");
        let (input, output, outside) = (d.join("in"), d.join("out"), d.join("outside"));
        for f in [&input, &output, &outside] {
            std::fs::create_dir_all(f).unwrap();
        }
        let design = plate_design(&Plate::default()).unwrap();
        std::fs::write(input.join("plate.solidcraft"), &design).unwrap();
        std::fs::write(outside.join("plate.solidcraft"), &design).unwrap();
        let mut s = Session { file_access: Roots::new(Some(&input), Some(&output)).unwrap(), ..Session::default() };
        let refused = |r: Result<Value>| {
            let e = r.unwrap_err().to_string();
            assert!(e.starts_with("automation"), "{e}");
        };
        refused(s.execute("geometry.open", &json!({"path": outside.join("plate.solidcraft")})));
        refused(s.execute("geometry.import.mesh", &json!({"path": outside.join("x.3mf")})));
        refused(s.execute("file.open", &json!({"path": outside.join("x.fieldcraft")})));
        refused(s.execute("file.recovery.list", &json!({"dir": outside})));
        refused(s.execute("file.recovery.discard", &json!({"all": true, "dir": input})));
        s.execute("geometry.open", &json!({"path": input.join("plate.solidcraft")})).unwrap();
        defined(&mut s);
        // A section, so deck.write has all it needs and reaches its path.
        s.execute("model.set.create", &json!({"name": "all", "elements": [1]})).unwrap();
        s.execute("model.section.create", &json!({"name": "solid", "set": "all", "material": "Steel", "kind": {"type": "solid"}})).unwrap();
        // Writes go to the write root only: not the read root, not elsewhere.
        refused(s.execute("file.save_as", &json!({"path": input.join("m.fieldcraft")})));
        refused(s.execute("file.save_as", &json!({"path": outside.join("m.fieldcraft")})));
        refused(s.execute("file.save_as", &json!({"path": output.join("../outside/m.fieldcraft")})));
        refused(s.execute("deck.write", &json!({"format": "abaqus", "path": outside.join("m.inp")})));
        let saved = output.join("m.fieldcraft");
        s.execute("file.save_as", &json!({"path": saved})).unwrap();
        assert!(saved.is_file());
        assert_eq!(std::fs::read_dir(&outside).unwrap().count(), 1, "nothing written outside");

        // The roots survive opening a file; a model whose design lies outside the read root
        // opens without it.
        let mut other = Session::default();
        other.execute("geometry.open", &json!({"path": outside.join("plate.solidcraft")})).unwrap();
        other.execute("file.save_as", &json!({"path": input.join("other.fieldcraft")})).unwrap();
        let r = s.execute("file.open", &json!({"path": input.join("other.fieldcraft")})).unwrap();
        assert!(r["warnings"][0].as_str().unwrap().contains("automation path rejected"), "{r}");
        assert!(s.design.is_none());
        assert!(s.file_access.is_confined());
        refused(s.execute("file.save", &json!({})));
        let _ = std::fs::remove_dir_all(&d);
    }

    /// The writers of the format 1 fixtures of `fieldcraft-io` (crates/io/tests/formats/v1).
    const V1_WRITERS: [&str; 5] = ["0.4.0", "0.5.0", "0.6.0", "0.7.0", "develop-2eb3515c"];

    /// The format 1 fixtures `writer` made, copied to a fresh folder.
    fn v1_fixtures(tag: &str, writer: &str) -> PathBuf {
        let d = dir(&format!("{tag}-{writer}"));
        let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../io/tests/formats/v1").join(writer);
        for e in std::fs::read_dir(src).unwrap() {
            let p = e.unwrap().path();
            if p.extension().is_some_and(|x| x == "fieldcraft" || x == "solidcraft") {
                std::fs::copy(&p, d.join(p.file_name().unwrap())).unwrap();
            }
        }
        d
    }

    /// Every file in `d` with its bytes and modification time.
    fn contents(d: &Path) -> Vec<(String, Vec<u8>, std::time::SystemTime)> {
        let mut all: Vec<_> = std::fs::read_dir(d)
            .unwrap()
            .map(|e| {
                let p = e.unwrap().path();
                (
                    p.file_name().unwrap().to_string_lossy().into_owned(),
                    std::fs::read(&p).unwrap(),
                    std::fs::metadata(&p).unwrap().modified().unwrap(),
                )
            })
            .collect();
        all.sort();
        all
    }

    /// Every format 1 fixture opens, migrated, as the reader gives it, with its geometry; opening
    /// writes nothing.
    #[test]
    fn format_1_models_open_migrated_and_opening_writes_nothing() {
        for writer in V1_WRITERS {
            format_1_models_of(writer);
        }
    }

    fn format_1_models_of(writer: &str) {
        let d = v1_fixtures("v1-open", writer);
        let before = contents(&d);
        let models: Vec<PathBuf> = before.iter().filter(|f| f.0.ends_with(".fieldcraft")).map(|f| d.join(&f.0)).collect();
        assert!(models.len() >= 6, "{writer}");
        for path in &models {
            let mut s = Session::default();
            let r = s.execute("file.open", &json!({"path": path})).unwrap();
            assert_eq!(r["migrated_from"], 1, "{r}");
            assert_eq!(r["warnings"], json!([]), "{}: {r}", path.display());
            assert_eq!(s.migrated_from, Some(1));
            assert!(!s.is_dirty(), "a migrated model has nothing unsaved: opening it again gives the same");
            let file = fieldcraft_io::read(&std::fs::read(path).unwrap()).unwrap();
            assert_eq!((&s.model, &s.journal, &s.view), (&file.model, &file.journal, &file.view));
        }
        assert_eq!(contents(&d), before, "opening wrote nothing");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// The first save over a format 1 file keeps a copy of it, says so, and is the only one.
    #[test]
    fn the_first_save_over_a_format_1_file_keeps_it() {
        let d = v1_fixtures("v1-save", "0.7.0");
        let path = d.join("hex-static.fieldcraft");
        let original = std::fs::read(&path).unwrap();
        let mut s = Session::default();
        s.execute("file.open", &json!({"path": path})).unwrap();
        s.execute("model.material.create", &json!({"name": "Alu", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        // Undo back to the migrated state keeps the mark.
        s.execute("edit.undo", &json!({})).unwrap();
        s.execute("edit.redo", &json!({})).unwrap();
        assert_eq!(s.migrated_from, Some(1));
        let r = s.execute("file.save", &json!({})).unwrap();
        let kept = d.join("hex-static.v1.fieldcraft");
        assert_eq!(r["backup"], json!(kept.to_string_lossy()), "{r}");
        assert!(r["summary"].as_str().unwrap().contains("kept as hex-static.v1.fieldcraft"), "{r}");
        assert_eq!(std::fs::read(&kept).unwrap(), original);
        assert_eq!(s.migrated_from, None);
        assert_eq!(fieldcraft_io::fs::format_of_file(&path), Some(fieldcraft_io::FORMAT_VERSION));
        let r = s.execute("file.save", &json!({})).unwrap();
        assert!(r.get("backup").is_none(), "{r}");
        assert_eq!(std::fs::read(&kept).unwrap(), original, "the copy stays as it was");

        // Save As elsewhere leaves the format 1 original alone and drops the mark.
        let other = d.join("every-kind.fieldcraft");
        let other_bytes = std::fs::read(&other).unwrap();
        s.execute("file.open", &json!({"path": other})).unwrap();
        s.execute("file.save_as", &json!({"path": d.join("copy.fieldcraft")})).unwrap();
        assert_eq!((s.migrated_from, std::fs::read(&other).unwrap()), (None, other_bytes));
        assert!(!d.join("every-kind.v1.fieldcraft").exists());
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Work recovered from an autosave of a migrated model still counts as migrated: the autosave
    /// is format 2, but the file it saves over is format 1.
    #[test]
    fn recovered_work_over_a_format_1_file_keeps_it_too() {
        let d = v1_fixtures("v1-recover", "0.6.0");
        let recovery = d.join("recovery");
        let path = d.join("tet-frequency.fieldcraft");
        let original = std::fs::read(&path).unwrap();
        let mut s = Session::default();
        s.execute("file.open", &json!({"path": path})).unwrap();
        s.execute("model.material.create", &json!({"name": "Alu", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        let mut autosave = crate::autosave::Autosaver::new(&recovery).unwrap();
        assert!(autosave.save(&s).unwrap());
        let id = autosave.id().to_string();
        drop(autosave);
        assert_eq!(std::fs::read(&path).unwrap(), original, "autosaving never touches the model's own file");

        let mut t = Session::default();
        let r = t.execute("file.recover", &json!({"id": id, "dir": recovery})).unwrap();
        assert_eq!(r["migrated_from"], 1, "{r}");
        assert!(t.is_dirty());
        assert_eq!(std::fs::read(&path).unwrap(), original, "recovering writes nothing");
        t.execute("file.save", &json!({})).unwrap();
        assert_eq!(std::fs::read(d.join("tet-frequency.v1.fieldcraft")).unwrap(), original);

        // The file was replaced by a later FieldCraft since the autosave: the work is not
        // "migrated from" it, but a save still keeps it. A confined session reads no file
        // outside its roots to find out.
        let later = {
            use std::io::{Read, Write};
            let mut zip = zip::ZipArchive::new(std::io::Cursor::new(&original)).unwrap();
            let mut out = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
            for i in 0..zip.len() {
                let mut f = zip.by_index(i).unwrap();
                let (name, mut b) = (f.name().unwrap().to_string(), Vec::new());
                f.read_to_end(&mut b).unwrap();
                if name == "fieldcraft.json" {
                    b = String::from_utf8(b)
                        .unwrap()
                        .replace("fieldcraft/1", &format!("fieldcraft/{}", fieldcraft_io::FORMAT_VERSION + 1))
                        .into_bytes();
                }
                out.start_file(name, zip::write::SimpleFileOptions::default()).unwrap();
                out.write_all(&b).unwrap();
            }
            out.finish().unwrap().into_inner()
        };
        std::fs::write(&path, &later).unwrap();
        let mut u = Session::default();
        let r = u.execute("file.recover", &json!({"id": id, "dir": recovery})).unwrap();
        assert_eq!(r["migrated_from"], Value::Null, "{r}");
        let r = u.execute("file.save", &json!({})).unwrap();
        assert!(r["backup"].as_str().unwrap().ends_with(&format!("tet-frequency.v{}.fieldcraft", fieldcraft_io::FORMAT_VERSION + 1)), "{r}");
        let elsewhere = d.join("elsewhere");
        std::fs::create_dir_all(&elsewhere).unwrap();
        let mut confined = Session { file_access: Roots::new(Some(&recovery), Some(&elsewhere)).unwrap(), ..Session::default() };
        std::fs::write(&path, &original).unwrap();
        let r = confined.execute("file.recover", &json!({"id": id, "dir": recovery})).unwrap();
        assert_eq!(r["migrated_from"], Value::Null, "the file outside the read root is not read: {r}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// A deck, an STL, a VTU or a CSV aimed at the model's own file is refused; the model stays.
    #[test]
    fn outputs_never_replace_the_model_file() {
        let d = dir("outputs");
        let mut s = Session::default();
        defined(&mut s);
        s.execute("model.set.create", &json!({"name": "all", "elements": [1]})).unwrap();
        s.execute("model.section.create", &json!({"name": "solid", "set": "all", "material": "Steel", "kind": {"type": "solid"}})).unwrap();
        s.execute("model.set.create", &json!({"name": "base", "nodes": [1, 2, 3]})).unwrap();
        s.execute("model.boundary_condition.create", &json!({"name": "fix", "type": "fixed", "target": {"set": "base"}})).unwrap();
        s.execute("model.step.create", &json!({"name": "Load", "procedure": {"type": "static"}})).unwrap();
        s.execute("model.step.activate", &json!({"step": "Load", "boundary_conditions": ["fix"]})).unwrap();
        let path = d.join("bracket.fieldcraft");
        s.execute("file.save_as", &json!({"path": path})).unwrap();
        let saved = std::fs::read(&path).unwrap();
        for (id, p) in [("deck.write", json!({"format": "abaqus", "path": path})), ("mesh.export.stl", json!({"path": path}))] {
            let e = s.execute(id, &p).unwrap_err().to_string();
            assert!(e.contains("FieldCraft model") && e.contains("choose another file name"), "{id}: {e}");
            assert_eq!(std::fs::read(&path).unwrap(), saved, "{id}");
        }
        s.execute("deck.write", &json!({"format": "abaqus", "path": d.join("bracket.inp")})).unwrap();
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn hostile_file_params_and_files_change_nothing() {
        let d = dir("hostile");
        let junk = d.join("junk.fieldcraft");
        std::fs::write(&junk, b"PK\x03\x04 not a model").unwrap();
        let mut s = Session::default();
        defined(&mut s);
        let before = (s.model.clone(), s.journal.clone(), s.revision());
        for (id, p) in [
            ("file.open", json!({"path": junk})),
            ("file.open", json!({"path": d})),
            ("file.open", json!({"path": d.join("missing.fieldcraft")})),
            ("file.open", json!({"path": "\u{0}"})),
            ("file.open", json!({"path": 5})),
            ("file.save", json!({"path": ""})),
            ("file.save", json!({"path": d})),
            ("file.save_as", json!({"path": d.join("no/such/dir/a.fieldcraft")})),
            ("file.save_as", json!({"path": "x".repeat(5000)})),
            ("file.save_as", json!({})),
            ("file.recover", json!({"id": "../../etc/passwd", "dir": d})),
            ("file.recover", json!({"dir": d})),
            ("file.recovery.discard", json!({"dir": d})),
            ("file.recovery.list", json!({"dir": 7})),
        ] {
            assert!(s.execute(id, &p).is_err(), "{id} {p}");
            assert_eq!((s.model.clone(), s.journal.clone(), s.revision()), before, "{id} {p}");
        }
        assert_eq!(std::fs::read_dir(&d).unwrap().count(), 1, "nothing written");
        let _ = std::fs::remove_dir_all(&d);
    }
}
