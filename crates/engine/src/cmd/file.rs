//! Files: open and save `.fieldcraft` models, and recover autosaved ones.
//!
//! A saved model refers to its design and imported meshes by path and content hash; opening it
//! opens them again. Paths inside the model file's folder are stored relative to it, so a folder
//! of a model and its design can be moved together. A geometry file that is missing, unreadable
//! or changed since the save does not stop the model from opening: it is reported in the result's
//! `warnings` (sets follow geometry names, so a remesh brings a changed design's mesh up to date).

use std::collections::BTreeMap;
use std::path::{Component, Path, PathBuf};
use std::sync::Arc;

use fieldcraft_geometry::tagged_mesh;
use fieldcraft_geometry::{Design, MAX_DESIGN_BYTES, open_design};
use fieldcraft_io::{Geometry, ModelFile, SourceFile};
use serde_json::{Value, json};

use fieldcraft_io::roots::{Access, Roots};

use super::{file_param, opt, req};
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
    params: r#"{"path": "bracket.fieldcraft"} — opens the model, its mesh and its journal, and opens its design and imported meshes again from their files (one that is missing or changed is listed in `warnings`; the model still opens), shows the saved view it was saved showing (`view`: its name, the camera the hosts aim at, `pending` until its step is solved) and its saved results (`results` lists each step's, `stale` when it was solved from another definition than the model's: kept and shown as stale, with a warning each); unsaved changes of the current model are discarded and the undo history starts empty. A file of an older format is migrated as it is read (`migrated_from` in the reply; what migrating changed in what it means, such as a mesh of the design in a model of metres, leads `warnings`) and nothing is written; a later format is refused. Names the design no longer has (a SolidCraft that names a face otherwise): a face the mesh lies on is found again by where the mesh is, its edges and vertices through their faces, and every definition follows (`renamed`, and a warning each; the model is then changed, to save); what is not found again is listed in `warnings`, to pick again, and in `stale` (as geometry.open's: `mesh_bodies`, `sets`, `seeds` and `references` [{owner, reference, why}]); `mismatched` lists mesh faces that lie on another face of the design than their name's (a name SolidCraft now gives another face)"#,
    form: Some(super::forms::FILE_OPEN),
    kind: CommandKind::Session,
    requires: &[],
    run: open,
};

pub const SAVE: CommandSpec = CommandSpec {
    id: "file.save",
    label: "Save Model",
    params: r#"{"path"?: "bracket.fieldcraft", "results"?: true, "geometry"?: "refer" | "embed"} — saves to `path`, else to the file the model was opened from or last saved to; written atomically, the previous version kept as `<file>.bak`. The file keeps every solved step's latest result with the definition it was solved from (`results: false` leaves them out: a smaller file, and one that diffs by its text parts). `"geometry": "embed"` also keeps a copy of the design and imported meshes inside the file (`geometry/<file>`, opened when a referenced file is missing: whoever gets the file gets the geometry too), `"refer"` only the references (the default, unless the model's file embeds them); the reply's `contents` says what the file holds. A file there that is not a model of the current format (one this model was migrated from, a later or damaged one, any other file) is first copied to `<stem>.v<N>.fieldcraft` (`<stem>.backup.<ext>` when it is not a model), never replaced later (`backup` in the reply)"#,
    form: Some(super::forms::FILE_SAVE),
    kind: CommandKind::Session,
    requires: &[],
    run: save,
};

pub const SAVE_AS: CommandSpec = CommandSpec {
    id: "file.save_as",
    label: "Save Model As",
    params: r#"{"path": "bracket.fieldcraft", "results"?: true, "geometry"?: "refer" | "embed"} — saves to a new file, which becomes the model's file (`results`, `geometry` and a file there that is not a model of the current format as file.save)"#,
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
    opened.saved_results = opened.results_stored;
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
/// opened. When the file cannot be read and the model embeds a copy of it (`embedded`), the copy
/// is opened instead, with a warning saying so.
fn reopen<T>(
    what: &str,
    stored: &SourceFile,
    base: Option<&Path>,
    roots: &Roots,
    max: u64,
    open: impl Fn(&[u8]) -> fieldcraft_geometry::Result<T>,
    embedded: Option<&[u8]>,
    warnings: &mut Vec<String>,
) -> (Option<T>, SourceFile) {
    let path = resolved_path(stored.path(), base);
    let kept = stored.with_path(path.clone()).unwrap_or_else(|_| stored.clone());
    let read = roots.check(Path::new(&path), Access::Read).map_err(EngineError::from).and_then(|checked| read_file(&checked.to_string_lossy(), max));
    let bytes = match (read, embedded) {
        (Ok(b), _) => b,
        (Err(e), Some(copy)) => {
            warnings.push(format!("{what} `{path}` was not opened ({e}); the copy embedded in the model was opened instead"));
            return match open(copy) {
                Ok(v) => (Some(v), kept),
                Err(e) => {
                    warnings.push(format!("the embedded copy of {what} `{path}` could not be opened ({e}); the model opens without it"));
                    (None, kept)
                }
            };
        }
        (Err(e), None) => {
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

/// The model expressed in millimetres, as the geometry is (the stored model is in SI).
pub(crate) fn frame_geometry(model: &fieldcraft_model::Model) -> std::result::Result<fieldcraft_model::Model, String> {
    let mm = crate::units::geometry();
    model.converted(&|q, v| mm.express(q, v).ok()).map_err(|e| e.to_string())
}

/// What a saved file keeps besides the model, its mesh, journal and view.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Keep {
    /// The solved results (`results.bin`).
    pub results: bool,
    /// Copies of the geometry files (`geometry/…`), besides the references to them.
    pub geometry: bool,
}

impl Session {
    /// The bytes of the `.fieldcraft` file holding this session's model, mesh, journal and view
    /// (no results: what autosaves, crash and bug reports keep), with geometry paths inside
    /// `base` (the folder the file goes to) stored relative to it.
    pub fn file_bytes(&self, base: Option<&Path>) -> Result<Vec<u8>> {
        self.file_bytes_with(base, Keep::default())
    }

    /// [`Session::file_bytes`], keeping what `keep` says.
    pub fn file_bytes_with(&self, base: Option<&Path>, keep: Keep) -> Result<Vec<u8>> {
        let stored = |s: &SourceFile| s.with_path(stored_path(s.path(), base));
        let geometry = Geometry {
            design: self.sources.design.as_ref().map(stored).transpose()?,
            meshes: self.sources.meshes.iter().map(stored).collect::<fieldcraft_io::Result<_>>()?,
        };
        let copies = if keep.geometry { self.geometry_copies()? } else { Vec::new() };
        let extras = fieldcraft_io::Extras {
            results: if keep.results { self.stored_results() } else { Vec::new() },
            embed: copies.iter().map(|(sha, b)| (sha.clone(), b.as_slice())).collect(),
            // The settings of steps the model still has.
            run: self.run.iter().filter(|(step, _)| self.model.step(step).is_some()).map(|(k, v)| (k.clone(), *v)).collect(),
        };
        Ok(fieldcraft_io::write_with(&self.model, &geometry, &self.journal, &self.journal_header, &self.view, &extras)?)
    }

    /// The bytes of every geometry file the model is built on, by content hash, to embed: the
    /// copy the model was opened with, else the file read again (it must still be what was
    /// opened).
    fn geometry_copies(&self) -> Result<Vec<(String, Arc<Vec<u8>>)>> {
        let mut out: Vec<(String, Arc<Vec<u8>>)> = Vec::new();
        for source in self.sources.design.iter().chain(&self.sources.meshes) {
            let sha = source.sha256().to_string();
            if out.iter().any(|(s, _)| *s == sha) {
                continue;
            }
            let bytes = match self.embedded.get(&sha) {
                Some(b) => b.clone(),
                None => {
                    let path = self.checked_path(source.path(), Access::Read)?;
                    let bytes = read_file(&path, fieldcraft_io::MAX_GEOMETRY_PART_BYTES)
                        .map_err(|e| EngineError::Params(format!("cannot embed `{}`: {e}", source.path())))?;
                    if fieldcraft_io::sha256_hex(&bytes) != sha {
                        return Err(EngineError::Params(format!(
                            "cannot embed `{}`: it changed since the model opened it (open the geometry again, then save)",
                            source.path()
                        )));
                    }
                    Arc::new(bytes)
                }
            };
            out.push((sha, bytes));
        }
        Ok(out)
    }

    /// A session holding what a file holds, its geometry opened again from the paths it names
    /// (relative ones against `base`) where `roots` lets it read them; problems with the
    /// geometry come back as warnings. The session's own file access is unconfined.
    /// The model's names are brought to the design's ([`super::names::reconcile`]): what was
    /// renamed and what is still stale come back in the report, and as warnings.
    pub(crate) fn from_file(file: ModelFile, base: Option<&Path>, roots: &Roots) -> (Session, Vec<String>, super::names::Report) {
        // What migrating an older file changed in what it means comes first.
        let mut warnings = file.notes.clone();
        // The geometry is in millimetres: it is cut and stitched by the model expressed in them.
        let mut framed = match frame_geometry(&file.model) {
            Ok(m) => m,
            Err(e) => {
                warnings.push(format!("the model's lengths cannot be expressed in millimetres ({e}); its partitions and stitches are not applied"));
                fieldcraft_model::Model::default()
            }
        };
        let mut sources = Geometry::default();
        let mut design: Option<Arc<Design>> = None;
        // A packed file's copies of its geometry files, used where a file is missing.
        let embedded: BTreeMap<String, Arc<Vec<u8>>> = file.embedded.into_iter().map(|(sha, b)| (sha, Arc::new(b))).collect();
        let copy = |stored: &SourceFile| stored.embedded().and_then(|_| embedded.get(stored.sha256())).map(|b| b.as_slice());
        let packed = file.geometry.design.iter().chain(&file.geometry.meshes).any(|s| s.embedded().is_some());
        let mut model = file.model;
        let mut names = super::names::Report::default();
        if let Some(stored) = &file.geometry.design {
            let (opened, source) = reopen("design", stored, base, roots, MAX_DESIGN_BYTES as u64, open_design, copy(stored), &mut warnings);
            design = opened.and_then(|d| {
                match super::datum::partitioned(&d, &framed).and_then(|cut| super::names::reconcile(&mut model, Some(&mut framed), &d, cut)) {
                    Ok((cut, report)) => {
                        names = report;
                        Some(Arc::new(cut))
                    }
                    Err(e) => {
                        warnings.push(format!("the design's partitions could not be applied ({e}); the model opens without the design"));
                        None
                    }
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
            warnings.extend(names.warnings());
            sources.design = Some(source.without_embedded());
        }
        let mut surfaces = Vec::new();
        for stored in &file.geometry.meshes {
            let (opened, source) =
                reopen("mesh", stored, base, roots, tagged_mesh::MAX_FILE_BYTES as u64, tagged_mesh::read_3mf, copy(stored), &mut warnings);
            surfaces.extend(opened.into_iter().flatten().map(Arc::new));
            sources.meshes.push(source.without_embedded());
        }
        let mut session = Session {
            model,
            design,
            surfaces,
            sources,
            journal: file.journal,
            journal_header: file.journal_header,
            revision: next_revision(),
            saved_view: file.view.clone(),
            view: file.view,
            migration_notes: file.notes,
            embedded,
            pack_geometry: packed,
            saved_run: file.run.clone(),
            run: file.run,
            ..Session::default()
        };
        let (_, failed) = super::stitch::restitch_surfaces(&mut session.surfaces, &framed);
        warnings.extend(failed.into_iter().map(|(n, why)| format!("stitch `{n}`: {why}")));
        // The results come back current when the model is still the one they were solved from.
        warnings.extend(file.results_note);
        for stored in file.results {
            let step = stored.result.step().to_string();
            if !session.keep_stored_result(stored.result.into(), stored.definition) {
                warnings.push(format!(
                    "the results of step `{step}` were solved from another definition than the model's: they are kept, shown as stale; solve again to bring them up to date"
                ));
            }
        }
        (session, warnings, names)
    }

    /// Replace this session with `opened`, keeping its undo memory limit, its work caps and its
    /// file access (and where its commands come from); the history starts empty.
    fn replace_with(&mut self, mut opened: Session) {
        opened.history = std::mem::take(&mut self.history);
        opened.limits = self.limits;
        opened.file_access = std::mem::take(&mut self.file_access);
        opened.crash_journal = self.crash_journal;
        opened.solver_diagnostics = self.solver_diagnostics;
        opened.origin = self.origin;
        opened.default_units = self.default_units.clone();
        // The hosts count camera requests across documents.
        opened.camera_requests = self.camera_requests;
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
        "results": s.result_steps().map(|(step, current)| json!({"step": step, "stale": !current})).collect::<Vec<_>>(),
        "migrated_from": s.migrated_from,
        "warnings": warnings,
    })
}

/// Replace the session with a `.fieldcraft` file's content; `path` is its own file (geometry
/// paths are relative to its folder), `saved` whether that file holds this content (else the
/// content is an autosave of it). Nothing is written.
fn load(s: &mut Session, bytes: &[u8], path: Option<String>, saved: bool) -> Result<Value> {
    let (file, format) = fieldcraft_io::read_versioned(bytes)?;
    let (mut opened, warnings, names) = Session::from_file(file, path.as_deref().and_then(|p| Path::new(p).parent()), &s.file_access);
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
    // Names rewritten to the design's are a change: the person saves it, opening writes nothing.
    if names.changed_model() {
        opened.revision = next_revision();
    }
    s.replace_with(opened);
    // The view the file was saved showing comes back (once its step is solved, when it is not).
    let mut warnings = warnings;
    let shown = match s.view.current.clone() {
        Some(name) => match super::views::show(s, &name) {
            Ok(r) => {
                warnings.extend(r["warnings"].as_array().into_iter().flatten().filter_map(|w| w.as_str().map(str::to_string)));
                json!({"name": name, "pending": r["pending"], "camera": r["camera"]})
            }
            Err(e) => {
                warnings.push(e.to_string());
                Value::Null
            }
        },
        None => Value::Null,
    };
    let mut reply = summary(s, &warnings);
    reply["view"] = shown;
    reply["stale"] = super::geometry::opened_stale(s, &names);
    let names = names.json();
    reply["renamed"] = names["renamed"].clone();
    reply["mismatched"] = names["mismatched"].clone();
    Ok(reply)
}

fn open(s: &mut Session, p: &Value) -> Result<Value> {
    let path = absolute(&file_param(s, p, "path", Access::Read)?);
    let bytes = read_file(&path, fieldcraft_io::MAX_FILE_BYTES)?;
    load(s, &bytes, Some(path.clone()), true).map_err(|e| EngineError::Params(format!("{path}: {e}")))
}

/// What a save of `s` keeps, from its parameters: `results` (by default every step's latest
/// result) and `geometry` (`embed` or `refer`; by default as the model's file was).
fn keep(s: &Session, p: &Value) -> Result<Keep> {
    let geometry = match opt::<String>(p, "geometry")?.as_deref() {
        None => s.pack_geometry,
        Some("embed") => true,
        Some("refer") => false,
        Some(other) => {
            return Err(EngineError::invalid(
                "geometry",
                format!("must be `embed` or `refer`, not `{}`", other.chars().take(40).collect::<String>()),
            ));
        }
    };
    Ok(Keep { results: opt::<bool>(p, "results")?.unwrap_or(true), geometry })
}

/// What a file of `s` saved with `keep` holds, in words: what a person sending it shares.
fn contents(s: &Session, keep: Keep) -> Vec<String> {
    let mut out = vec!["the model definition, its mesh and the commands that built it".to_string()];
    let steps: Vec<&str> = s.result_steps().map(|(n, _)| n).collect();
    if keep.results && !steps.is_empty() {
        let count = if steps.len() == 1 { "1 step".to_string() } else { format!("{} steps", steps.len()) };
        out.push(format!("the results of {count} ({})", steps.join(", ")));
    }
    let named: Vec<&str> = s.sources.design.iter().chain(&s.sources.meshes).map(|f| f.path()).collect();
    if !named.is_empty() && keep.geometry {
        out.push(format!(
            "a copy of its geometry files ({}): whoever gets the file gets the design and meshes too",
            named.iter().map(|p| file_name(p)).collect::<Vec<_>>().join(", ")
        ));
    } else if !named.is_empty() {
        out.push(format!("references to its geometry files, not the files: {}", named.join(", ")));
    }
    out
}

/// Save to `path` (as given by the caller), keeping what `keep` says.
fn save_to(s: &mut Session, path: &str, keep: Keep) -> Result<Value> {
    let target = absolute(path);
    if std::fs::metadata(&target).is_ok_and(|m| !m.is_file()) {
        return Err(EngineError::Params(format!("`{target}` exists and is not a regular file")));
    }
    let bytes = s.file_bytes_with(Path::new(&target).parent(), keep)?;
    let kept = fieldcraft_io::fs::save_model(Path::new(&target), &bytes)?;
    let holds = contents(s, keep);
    s.mark_saved(target.clone());
    s.pack_geometry = keep.geometry;
    let mut reply = json!({"path": target, "bytes": bytes.len(), "contents": holds});
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
    let keep = keep(s, p)?;
    save_to(s, &path, keep)
}

fn save_as(s: &mut Session, p: &Value) -> Result<Value> {
    let path = file_param(s, p, "path", Access::Write)?;
    let keep = keep(s, p)?;
    save_to(s, &path, keep)
}

/// The last part of `path`.
fn file_name(path: &str) -> &str {
    path.rsplit(['/', '\\']).next().unwrap_or(path)
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
            // A mesh of the millimetre design in a model of metres or inches opens with a note
            // saying what changed, kept with the migrated-from note until the first save.
            if ["design-metres", "design-inches", "mixed-metres"].iter().any(|n| path.file_stem().is_some_and(|s| s == *n)) {
                let w = r["warnings"].to_string();
                assert!(w.contains("was made from a SolidCraft design") && w.contains("results will differ"), "{w}");
                assert_eq!(s.migration_notes.len(), 1, "{}", path.display());
            } else {
                assert_eq!(r["warnings"], json!([]), "{}: {r}", path.display());
            }
            assert_eq!(s.migrated_from, Some(1));
            assert!(!s.is_dirty(), "a migrated model has nothing unsaved: opening it again gives the same");
            let file = fieldcraft_io::read(&std::fs::read(path).unwrap()).unwrap();
            assert_eq!((&s.model, &s.journal, &s.view), (&file.model, &file.journal, &file.view));
        }
        assert_eq!(contents(&d), before, "opening wrote nothing");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// The journal of every migrated fixture replays to its migrated model: the step keeps the
    /// journal's numbers as they were typed and names their units, and replay reads them so.
    #[test]
    fn every_migrated_journal_replays_to_its_migrated_model() {
        // As written, and relabelled into every system (the defaults format 2 filled in mean
        // other numbers in each): bit for bit.
        const SYSTEMS: [Option<&str>; 6] =
            [None, Some("mm-t-s-N-MPa"), Some("m-kg-s-N-Pa"), Some("mm-kg-ms-kN-GPa"), Some("mm-g-ms-N-MPa"), Some("in-lbf-s-psi")];
        let (mut checked, mut worst) = (0, 0.0f64);
        for writer in V1_WRITERS {
            let d = v1_fixtures("v1-replay", writer);
            let mut models: Vec<PathBuf> =
                std::fs::read_dir(&d).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|x| x == "fieldcraft")).collect();
            models.sort();
            for original in &models {
                for system in SYSTEMS {
                    let path = match system {
                        None => original.clone(),
                        Some(system) => {
                            let p = d.join(format!("relabelled-{system}.fieldcraft"));
                            std::fs::write(&p, relabelled(&std::fs::read(original).unwrap(), system)).unwrap();
                            p
                        }
                    };
                    let what = format!("{writer} {} as {system:?}", original.display());
                    let mut opened = Session::default();
                    opened.execute("file.open", &json!({"path": path})).unwrap_or_else(|e| panic!("{what}: {e}"));
                    let file = fieldcraft_io::read(&std::fs::read(&path).unwrap()).unwrap();
                    let h = &file.journal_header;
                    let commands: Vec<Value> = file.journal.iter().map(|(c, p)| json!({"command": c, "params": p})).collect();
                    let mut script = json!({"units": h.units, "model_units": h.model_units, "commands": commands});
                    // A design is named relative to the model's folder.
                    for c in script["commands"].as_array_mut().unwrap() {
                        if c["command"] == "geometry.open"
                            && let Some(p) = c["params"]["path"].as_str()
                        {
                            c["params"]["path"] = json!(d.join(p));
                        }
                    }
                    let mut replayed = Session::default();
                    replayed.run_script(&script).unwrap_or_else(|e| panic!("{what}: {e}"));
                    // Bit for bit, but for coordinates a mesher computed after another one on
                    // the same mesh (a relaxation after a tet mesh): it sees the millimetre image
                    // of the stored metres, within an ulp of what 0.7.0 saw.
                    let (mut a, mut b) = (serde_json::to_value(&replayed.model).unwrap(), serde_json::to_value(&opened.model).unwrap());
                    worst = worst.max(settle_coordinates(&mut a, &mut b, ""));
                    settle_prose(&mut a, &mut b, "");
                    if a != b {
                        panic!("{what}: the replay differs at {}", difference(&a, &b, ""));
                    }
                    checked += 1;
                }
            }
            let _ = std::fs::remove_dir_all(&d);
        }
        assert!(checked >= 100, "{checked}");
        // Of models about 0.1 m across: a few ulps.
        assert!(worst < 1e-15, "{worst}");
    }

    /// Where `a` and `b` hold the same node or relaxation coordinate up to 1e-12 (relative to
    /// 1 m), make them equal; the largest such difference.
    fn settle_coordinates(a: &mut Value, b: &mut Value, at: &str) -> f64 {
        let coordinate = at.starts_with("/nodes/") || (at.starts_with("/relaxations/") && (at.contains("/at/") || at.contains("/moved/")));
        match (a, b) {
            (Value::Number(x), Value::Number(y)) if coordinate => {
                let (p, q) = (x.as_f64().unwrap_or(f64::NAN), y.as_f64().unwrap_or(f64::NAN));
                let d = (p - q).abs();
                if d > 0.0 && d < 1e-12 {
                    *x = y.clone();
                    return d;
                }
                0.0
            }
            (Value::Object(x), Value::Object(y)) => {
                x.iter_mut().filter_map(|(k, v)| y.get_mut(k).map(|w| settle_coordinates(v, w, &format!("{at}/{k}")))).fold(0.0, f64::max)
            }
            (Value::Array(x), Value::Array(y)) => {
                x.iter_mut().zip(y.iter_mut()).enumerate().map(|(i, (v, w))| settle_coordinates(v, w, &format!("{at}/{i}"))).fold(0.0, f64::max)
            }
            _ => 0.0,
        }
    }

    /// Where `a` and `b` hold a mesh plan's prose (a feature's note, a cut's or unmapped
    /// feature's reason, the plan's notes), make them equal: the text is this release's wording
    /// and number format, the release that wrote the file had its own; the plan's data is
    /// compared as it is.
    fn settle_prose(a: &mut Value, b: &mut Value, at: &str) {
        let prose = at.starts_with("/mesh_plans/") && (at.ends_with("/note") || at.ends_with("/reason") || at.rsplit('/').nth(1) == Some("notes"));
        match (a, b) {
            (Value::String(x), Value::String(y)) if prose => *x = y.clone(),
            (Value::Object(x), Value::Object(y)) => x.iter_mut().for_each(|(k, v)| {
                if let Some(w) = y.get_mut(k) {
                    settle_prose(v, w, &format!("{at}/{k}"));
                }
            }),
            (Value::Array(x), Value::Array(y)) => {
                x.iter_mut().zip(y.iter_mut()).enumerate().for_each(|(i, (v, w))| settle_prose(v, w, &format!("{at}/{i}")))
            }
            _ => {}
        }
    }

    /// Where `a` and `b` first differ, with both values there.
    fn difference(a: &Value, b: &Value, at: &str) -> String {
        match (a, b) {
            (Value::Object(x), Value::Object(y)) => x.keys().chain(y.keys()).find(|k| x.get(*k) != y.get(*k)).map_or(format!("{at}: (same)"), |k| {
                difference(x.get(k).unwrap_or(&Value::Null), y.get(k).unwrap_or(&Value::Null), &format!("{at}/{k}"))
            }),
            (Value::Array(x), Value::Array(y)) if x.len() == y.len() => x
                .iter()
                .zip(y)
                .enumerate()
                .find(|(_, (p, q))| p != q)
                .map_or(format!("{at}: (same)"), |(i, (p, q))| difference(p, q, &format!("{at}/{i}"))),
            _ => format!("{at}: {a} vs {b}"),
        }
    }

    /// The format 1 file `bytes` declaring `system`: its model's `units` and the journal's
    /// `model.units` say so, the manifest's sizes and hashes follow.
    fn relabelled(bytes: &[u8], system: &str) -> Vec<u8> {
        use std::io::{Read, Write};
        let mut zip = zip::ZipArchive::new(std::io::Cursor::new(bytes)).unwrap();
        let mut parts: Vec<(String, Vec<u8>)> = (0..zip.len())
            .map(|i| {
                let mut f = zip.by_index(i).unwrap();
                let mut b = Vec::new();
                f.read_to_end(&mut b).unwrap();
                (f.name().unwrap().to_string(), b)
            })
            .collect();
        let mut hashes = Vec::new();
        for (name, b) in &mut parts {
            let edited = match name.as_str() {
                "model.json" => {
                    let mut m: Value = serde_json::from_slice(b).unwrap();
                    m["units"] = json!(system);
                    serde_json::to_vec(&m).unwrap()
                }
                "journal.json" => {
                    // A model gets its system from a command (a file of a system names it in
                    // its journal): changed where there is one, the first command otherwise.
                    let mut j: Value = serde_json::from_slice(b).unwrap();
                    let list = j.as_array_mut().unwrap();
                    let mut named = false;
                    for c in list.iter_mut() {
                        if c["command"] == "model.units" {
                            c["params"]["system"] = json!(system);
                            named = true;
                        }
                    }
                    if !named {
                        list.insert(0, json!({"command": "model.units", "params": {"system": system}}));
                    }
                    serde_json::to_vec(&j).unwrap()
                }
                _ => continue,
            };
            hashes.push((name.clone(), json!({"bytes": edited.len(), "sha256": fieldcraft_io::sha256_hex(&edited)})));
            *b = edited;
        }
        let mut manifest: Value = serde_json::from_slice(&parts[0].1).unwrap();
        for (name, h) in hashes {
            manifest["parts"][&name] = h;
        }
        parts[0].1 = serde_json::to_vec(&manifest).unwrap();
        let mut out = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        for (name, b) in &parts {
            out.start_file(name, zip::write::SimpleFileOptions::default().compression_method(zip::CompressionMethod::Stored)).unwrap();
            out.write_all(b).unwrap();
        }
        out.finish().unwrap().into_inner()
    }

    /// Replaying a 0.7.0 journal reads a geometry command's numbers as millimetres only in a
    /// model built on a design: a block model of metres moves its nodes in metres, and a command
    /// that acted on a block beside a mesh of the design is refused rather than guessed.
    #[test]
    fn legacy_geometry_commands_replay_in_the_unit_they_were_typed_in() {
        let d = v1_fixtures("v1-legacy", "0.7.0");
        let replay = |name: &str, system: &str, extra: Value| -> Result<Session> {
            let path = d.join(format!("{name}-{system}.fieldcraft"));
            std::fs::write(&path, relabelled(&std::fs::read(d.join(format!("{name}.fieldcraft"))).unwrap(), system)).unwrap();
            let file = fieldcraft_io::read(&std::fs::read(&path).unwrap()).unwrap();
            let mut commands: Vec<Value> = file.journal.iter().map(|(c, p)| json!({"command": c, "params": p})).collect();
            for c in &mut commands {
                if c["command"] == "geometry.open" {
                    c["params"]["path"] = json!(d.join(c["params"]["path"].as_str().unwrap()));
                }
            }
            commands.push(extra);
            let h = &file.journal_header;
            let mut s = Session::default();
            s.run_script(&json!({"units": h.units, "model_units": h.model_units, "commands": commands}))?;
            Ok(s)
        };
        let s = replay(
            "hex-static",
            "m-kg-s-N-Pa",
            json!({"command": "mesh.node.move", "params": {"nodes": [20], "by": [1, 0, 0], "on_geometry": false, "guard": "none"}}),
        )
        .unwrap();
        assert_eq!(s.model.node(20).unwrap().xyz[0], 41.0);
        // Node 1 is the block `Bar`'s (typed in metres) beside the design's `Case`.
        let e = replay(
            "mixed-metres",
            "m-kg-s-N-Pa",
            json!({"command": "mesh.node.move", "params": {"nodes": [1], "by": [1, 0, 0], "on_geometry": false, "guard": "none"}}),
        )
        .unwrap_err()
        .to_string();
        assert!(e.contains("from FieldCraft 0.7 acted on node 1 of block body `Bar`, typed in m") && e.contains("cannot be replayed exactly"), "{e}");
        // So is one over the whole mesh, given no target.
        let e = replay("mixed-metres", "m-kg-s-N-Pa", json!({"command": "mesh.node.relax", "params": {"mode": "strong"}})).unwrap_err().to_string();
        assert!(e.contains("acted on every body of the mesh, the block body `Bar` (typed in m)") && e.contains("cannot be replayed exactly"), "{e}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// A fluid model of 0.7.0 typed in kelvin keeps its lattice counting from 0 K (the scheme is
    /// not invariant under a shift of the temperature scale), whatever it is shown in: the same
    /// heat in SI shown in K or in °C.
    #[test]
    fn a_kelvin_fluid_model_solves_as_it_did_whatever_it_is_shown_in() {
        let d = v1_fixtures("v1-fluid", "0.7.0");
        let mut s = Session::default();
        s.execute("file.open", &json!({"path": d.join("fluid-metres.fieldcraft")})).unwrap();
        assert_eq!(s.model.fluid_domain("Duct").and_then(|f| f.thermal.as_ref()).map(|t| t.temperature_zero), Some(0.0));
        let heat = |s: &mut Session| -> Vec<f64> {
            s.execute("solve.fluid", &json!({})).unwrap();
            let Some(crate::StepResult::Fluid(r)) = s.result("Flow").map(|r| r.as_ref()) else { panic!("no fluid result") };
            r.heat.iter().map(|h| *h.power.last().unwrap()).collect()
        };
        let in_kelvin = heat(&mut s);
        s.execute("model.units", &json!({"system": "mm-t-s-N-MPa"})).unwrap();
        assert_eq!(heat(&mut s), in_kelvin, "how it is shown changes nothing solved");
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

    /// `bytes` with part `name` edited by `edit` and the manifest's size and hash following it.
    fn with_part_edited(bytes: &[u8], name: &str, edit: impl FnOnce(&mut Vec<u8>)) -> Vec<u8> {
        use std::io::{Read, Write};
        let mut zip = zip::ZipArchive::new(std::io::Cursor::new(bytes)).unwrap();
        let mut parts: Vec<(String, Vec<u8>)> = (0..zip.len())
            .map(|i| {
                let mut f = zip.by_index(i).unwrap();
                let mut b = Vec::new();
                f.read_to_end(&mut b).unwrap();
                (f.name().unwrap().to_string(), b)
            })
            .collect();
        let part = parts.iter_mut().find(|(n, _)| n == name).unwrap();
        edit(&mut part.1);
        let entry = json!({"bytes": part.1.len(), "sha256": fieldcraft_io::sha256_hex(&part.1)});
        let mut manifest: Value = serde_json::from_slice(&parts[0].1).unwrap();
        manifest["parts"][name] = entry;
        parts[0].1 = serde_json::to_vec(&manifest).unwrap();
        let mut out = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        for (name, b) in &parts {
            out.start_file(name, zip::write::SimpleFileOptions::default().compression_method(zip::CompressionMethod::Stored)).unwrap();
            out.write_all(b).unwrap();
        }
        out.finish().unwrap().into_inner()
    }

    /// Every bit of a step's result (Debug prints every float exactly).
    fn bits(s: &Session, step: &str) -> String {
        format!("{:?}", s.latest_result(step).map(|(r, _)| r))
    }

    /// Solved results are kept in the file (#644): a fresh session opens them current and shows
    /// them without solving, bit for bit; an edit marks them stale and undo makes them current
    /// again; a file saved with stale results opens them stale; `results: false` leaves them out.
    #[test]
    fn solved_results_reopen_without_solving() {
        let d = dir("results");
        let mut s = Session::default();
        crate::cmd::solve::tests::cantilever(&mut s);
        s.execute("solve.static", &json!({})).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        let path = d.join("beam.fieldcraft");
        s.execute("file.save_as", &json!({"path": path, "results": false})).unwrap();
        let lean = std::fs::read(&path).unwrap();
        assert!(fieldcraft_io::read(&lean).unwrap().results.is_empty(), "`results: false` keeps none");
        s.execute("solve.static", &json!({})).unwrap();
        assert!(s.is_dirty(), "a solve is unsaved: files keep results");
        let r = s.execute("file.save", &json!({})).unwrap();
        assert!(r["contents"].to_string().contains("the results of 2 steps (Bend, Modes)"), "{r}");
        assert!(!s.is_dirty());
        assert!(std::fs::read(&path).unwrap().len() > lean.len());

        let mut t = Session::default();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert_eq!(r["results"], json!([{"step": "Bend", "stale": false}, {"step": "Modes", "stale": false}]), "{r}");
        assert_eq!(r["warnings"], json!([]));
        assert!(!t.is_dirty(), "opening is not solving");
        for step in ["Bend", "Modes"] {
            assert_eq!(bits(&t, step), bits(&s, step), "{step}");
        }
        // Shown without solving again.
        let shown = t.execute("results.contour", &json!({"step": "Bend", "field": "von_mises"})).unwrap();
        assert_eq!(shown, s.execute("results.contour", &json!({"step": "Bend", "field": "von_mises"})).unwrap());
        assert_eq!(t.execute("results.modes", &json!({"step": "Modes", "mode": 2})).unwrap()["stale"], false);
        assert!(t.execute("results.reactions", &json!({"step": "Bend"})).is_ok());
        // An edit marks them stale; undo brings them back.
        t.execute("model.material.edit", &json!({"name": "Mat", "elastic": {"type": "isotropic", "young": 2000, "poisson": 0}, "density": 1}))
            .unwrap();
        assert_eq!(t.result_steps().collect::<Vec<_>>(), vec![("Bend", false), ("Modes", false)]);
        t.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(t.result_steps().collect::<Vec<_>>(), vec![("Bend", true), ("Modes", true)]);
        // Saved with the model changed since: they open stale, kept and said so.
        t.execute("edit.redo", &json!({})).unwrap();
        t.execute("file.save", &json!({})).unwrap();
        let mut u = Session::default();
        let r = u.execute("file.open", &json!({"path": path})).unwrap();
        assert_eq!(r["results"], json!([{"step": "Bend", "stale": true}, {"step": "Modes", "stale": true}]), "{r}");
        assert!(r["warnings"].to_string().contains("results of step `Bend` were solved from another definition"), "{r}");
        assert_eq!(bits(&u, "Bend"), bits(&s, "Bend"), "kept as they were");
        assert_eq!(u.execute("results.contour", &json!({"step": "Bend", "field": "von_mises"})).unwrap()["stale"], true);
        // Editing back to the definition they were solved from makes them current.
        u.execute("model.material.edit", &json!({"name": "Mat", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0}, "density": 1}))
            .unwrap();
        assert_eq!(u.result_steps().collect::<Vec<_>>(), vec![("Bend", true), ("Modes", true)]);
        // Autosaves keep no results.
        assert!(fieldcraft_io::read(&u.file_bytes(None).unwrap()).unwrap().results.is_empty());

        // Results of another results encoding: the model opens without them, saying why.
        let later = with_part_edited(&std::fs::read(&path).unwrap(), "results.bin", |b| b[6..8].copy_from_slice(b"77"));
        let other = d.join("later.fieldcraft");
        std::fs::write(&other, later).unwrap();
        let mut v = Session::default();
        let r = v.execute("file.open", &json!({"path": other})).unwrap();
        assert_eq!(r["results"], json!([]));
        assert!(r["warnings"][0].as_str().unwrap().contains("another results encoding (MCRSLT77"), "{r}");
        assert_eq!(v.model, t.model);
        // A damaged results part is a damaged file: refused, the session unchanged.
        let damaged = with_part_edited(&std::fs::read(&path).unwrap(), "results.bin", |b| b.truncate(b.len() - 3));
        std::fs::write(&other, damaged).unwrap();
        let before = v.model.clone();
        assert!(v.execute("file.open", &json!({"path": other})).unwrap_err().to_string().contains("results.bin"));
        assert_eq!(v.model, before);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// A fluid step (cell fields, histories, heat) and an explicit step (frames, histories,
    /// energies) reopen bit for bit, and their frames, slices and streamlines show without
    /// solving again.
    #[test]
    fn fluid_and_explicit_results_reopen_bit_for_bit() {
        let d = dir("results-kinds");
        let mut fluid = Session::default();
        crate::cmd::solve::tests::heated_couette(&mut fluid);
        fluid.execute("solve.fluid", &json!({})).unwrap();
        let mut explicit = Session::default();
        crate::cmd::solve::tests::explicit_impact(&mut explicit);
        explicit.execute("solve.explicit", &json!({"step": "Impact"})).unwrap();
        for (s, step, name) in [(&mut fluid, "Shear", "couette"), (&mut explicit, "Impact", "impact")] {
            let path = d.join(format!("{name}.fieldcraft"));
            s.execute("file.save_as", &json!({"path": path})).unwrap();
            // Shown in the same units (the explicit model is in SI).
            let mut t = Session { default_units: s.default_units.clone(), ..Session::default() };
            let r = t.execute("file.open", &json!({"path": path})).unwrap();
            assert_eq!(r["results"], json!([{"step": step, "stale": false}]), "{r}");
            assert_eq!(bits(&t, step), bits(s, step), "{step}");
            let frame = json!({"step": step, "frame": 1});
            assert_eq!(t.execute("results.frame", &frame).unwrap(), s.execute("results.frame", &frame).unwrap());
        }
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": d.join("couette.fieldcraft")})).unwrap();
        for (id, p) in [
            ("results.frame", json!({"step": "Shear", "frame": 1})),
            ("results.slice", json!({"step": "Shear", "field": "temperature", "axis": "z"})),
            ("results.streamlines", json!({"step": "Shear", "count": 4})),
            ("solve.probe", json!({"step": "Shear", "points": [[0.06, 0.5, 0.06]]})),
        ] {
            assert_eq!(t.execute(id, &p).unwrap(), fluid.execute(id, &p).unwrap(), "{id}");
        }
        let _ = std::fs::remove_dir_all(&d);
    }

    /// The results fixtures: files whose results.bin this version's encoding wrote, one per kind
    /// of result (crates/engine/tests/results/).
    fn results_fixtures() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/results").join(String::from_utf8_lossy(fieldcraft_io::results::MAGIC).as_ref())
    }

    /// Writes the results fixtures: run once for each new results encoding (a new magic), then
    /// the files are never edited.
    #[test]
    #[ignore = "writes the results fixtures for a new results encoding"]
    fn write_the_results_fixtures() {
        let d = results_fixtures();
        std::fs::create_dir_all(&d).unwrap();
        let mut implicit = Session::default();
        crate::cmd::solve::tests::cantilever(&mut implicit);
        implicit.execute("solve.run", &json!({})).unwrap();
        let mut fluid = Session::default();
        crate::cmd::solve::tests::heated_couette(&mut fluid);
        fluid.execute("solve.fluid", &json!({})).unwrap();
        let mut explicit = Session::default();
        crate::cmd::solve::tests::explicit_impact(&mut explicit);
        explicit.execute("solve.explicit", &json!({"step": "Impact"})).unwrap();
        for (s, name) in [(&mut implicit, "cantilever"), (&mut fluid, "couette"), (&mut explicit, "impact")] {
            s.execute("file.save_as", &json!({"path": d.join(format!("{name}.fieldcraft"))})).unwrap();
        }
    }

    /// Results this encoding wrote are read by it: a change to the solvers' result types that
    /// would misread them fails here, and is a new results encoding (a new magic, new fixtures).
    #[test]
    fn results_of_this_encoding_read_back() {
        let d = results_fixtures();
        for (name, steps) in [("cantilever", &["Bend", "Modes"][..]), ("couette", &["Shear"]), ("impact", &["Impact"])] {
            let bytes = std::fs::read(d.join(format!("{name}.fieldcraft"))).unwrap();
            let file = fieldcraft_io::read(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(file.results_note.is_none(), "{name}: {:?}", file.results_note);
            let read: Vec<&str> = file.results.iter().map(|r| r.result.step()).collect();
            assert_eq!(read, steps, "{name}");
            let mut s = Session::default();
            s.execute("file.open", &json!({"path": d.join(format!("{name}.fieldcraft"))})).unwrap();
            assert_eq!(s.result_steps().count(), steps.len(), "{name}");
        }
    }

    /// #647: a model saved with its geometry embedded opens on a machine with nothing else, with
    /// its design, and stays packed; references stay the default, and an embedded copy is used
    /// only where the file a reference names is missing.
    #[test]
    fn a_packed_model_opens_with_its_geometry_and_nothing_else() {
        let d = dir("packed");
        let design = d.join("plate.solidcraft");
        std::fs::write(&design, plate_design(&Plate::default()).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": design})).unwrap();
        defined(&mut s);
        let referring = d.join("refer.fieldcraft");
        let r = s.execute("file.save_as", &json!({"path": referring})).unwrap();
        assert!(r["contents"].to_string().contains("references to its geometry files, not the files"), "{r}");
        assert!(fieldcraft_io::read(&std::fs::read(&referring).unwrap()).unwrap().embedded.is_empty(), "references by default");
        let packed = d.join("packed.fieldcraft");
        let r = s.execute("file.save_as", &json!({"path": packed, "geometry": "embed"})).unwrap();
        assert!(r["contents"].to_string().contains("a copy of its geometry files (plate.solidcraft): whoever gets the file gets the design"), "{r}");
        let file = fieldcraft_io::read(&std::fs::read(&packed).unwrap()).unwrap();
        assert_eq!(file.geometry.design.as_ref().and_then(|g| g.embedded()), Some("geometry/plate.solidcraft"));

        // With the design next to it, the design file is opened, not the copy.
        let mut t = Session::default();
        let r = t.execute("file.open", &json!({"path": packed})).unwrap();
        assert_eq!(r["warnings"], json!([]), "{r}");
        assert!(t.pack_geometry && t.design.is_some());

        // Alone in another folder: the copy, said so.
        let elsewhere = d.join("elsewhere");
        std::fs::create_dir_all(&elsewhere).unwrap();
        let alone = elsewhere.join("packed.fieldcraft");
        std::fs::copy(&packed, &alone).unwrap();
        std::fs::remove_file(&design).unwrap();
        let mut u = Session::default();
        let r = u.execute("file.open", &json!({"path": alone})).unwrap();
        assert!(r["warnings"][0].as_str().unwrap().contains("the copy embedded in the model was opened instead"), "{r}");
        assert_eq!(u.design.as_ref().map(|d| d.bodies.len()), s.design.as_ref().map(|d| d.bodies.len()));
        assert_eq!(u.model, s.model);
        // Saved again it stays packed, from the copy (the design file is gone).
        u.execute("file.save", &json!({})).unwrap();
        let again = fieldcraft_io::read(&std::fs::read(&alone).unwrap()).unwrap();
        assert_eq!(again.embedded, file.embedded);
        // Saved referring, the copy goes; without the design file it then opens without it.
        u.execute("file.save", &json!({"geometry": "refer"})).unwrap();
        let mut v = Session::default();
        let r = v.execute("file.open", &json!({"path": alone})).unwrap();
        assert!(v.design.is_none() && r["warnings"][0].as_str().unwrap().contains("the model opens without it"), "{r}");
        // Nothing to embed from: the design is gone and the session has no copy.
        let e = v.execute("file.save", &json!({"geometry": "embed"})).unwrap_err().to_string();
        assert!(e.contains("cannot embed"), "{e}");
        assert!(v.execute("file.save", &json!({"geometry": "zip"})).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }
}
