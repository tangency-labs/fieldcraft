//! Filesystem boundary for immutable source bundles. No solver data is read from other files.
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_decks::abaqus::read::bundle::SourceBundle;
use fieldcraft_io::roots::{Access, Roots};
use serde_json::{Value, json};
use std::path::Path;

pub const WRITE_SOURCE: CommandSpec = CommandSpec {
    id: "deck.source.write",
    label: "Write Original Abaqus Source",
    params: r#"{"path":"original.inp"} — replays the original imported Abaqus root file and INCLUDE dependencies, byte for byte, from the current session. This is source recovery, not export of model edits; use deck.write for the edited model. The root filename can change; INCLUDE paths stay relative to its parent. All destinations are checked before writing; includes are written first and root last. Original source is stored by file.save and recovered by file.open without rereading original files. Existing files may be replaced; each file write is atomic, the bundle is not a filesystem transaction."#,
    form: Some(super::super::forms::DECK_SOURCE),
    kind: CommandKind::Session,
    requires: &[],
    run: write,
};

pub(super) fn load(s: &Session, path: &str, text: &str) -> Result<SourceBundle> {
    let absolute = super::super::file::absolute(path);
    let path = Path::new(&absolute);
    let base = path.parent().ok_or_else(|| EngineError::invalid("path", "input needs a parent directory"))?;
    let root = path.file_name().and_then(|s| s.to_str()).ok_or_else(|| EngineError::invalid("path", "input filename is not UTF-8"))?;
    // Always impose the deck-directory boundary, even for unconfined interactive sessions.
    let confined = Roots::new(Some(base), None)?;
    Ok(SourceBundle::load(root, text, |name| {
        let read = || -> Result<String> {
            let path = confined.check(&base.join(name), Access::Read)?;
            let path = s.file_access.check(&path, Access::Read)?;
            let bytes = super::super::file::read_file(&path.to_string_lossy(), fieldcraft_decks::abaqus::read::MAX_INPUT_BYTES as u64)?;
            String::from_utf8(bytes).map_err(|_| EngineError::invalid("path", "INCLUDE file is not UTF-8"))
        };
        read().map_err(|e| fieldcraft_decks::Error::Parse { line: 1, message: e.to_string() })
    })?)
}

fn write(s: &mut Session, p: &Value) -> Result<Value> {
    let path = super::super::file_param(s, p, "path", Access::Write)?;
    let bundle =
        s.deck_source.as_ref().ok_or_else(|| EngineError::invalid("source", "no retained Abaqus source in this session; run deck.read first"))?;
    bundle.validate()?;
    let absolute = super::super::file::absolute(&path);
    let root = Path::new(&absolute);
    let base = root.parent().ok_or_else(|| EngineError::invalid("path", "output needs a parent directory"))?;
    // Existing parent required; no output outside it, including symlinks.
    let confined = Roots::new(None, Some(base))?;
    let mut destinations = Vec::new();
    let mut unique = std::collections::BTreeSet::new();
    let mut bytes = 0usize;
    for (name, doc) in &bundle.files {
        let destination = if name == &bundle.root { root.to_path_buf() } else { base.join(name) };
        let destination = confined.check(&destination, Access::Write)?;
        let destination = s.file_access.check(&destination, Access::Write)?;
        if !unique.insert(destination.clone()) {
            return Err(EngineError::invalid("path", "root output collides with an INCLUDE destination"));
        }
        if destination.exists() && !destination.is_file() {
            return Err(EngineError::invalid("path", "source destination is not a regular file"));
        }
        bytes = bytes.saturating_add(doc.len());
        destinations.push((name == &bundle.root, destination, doc));
    }
    destinations.sort_by_key(|(root, _, _)| *root);
    for (_, path, doc) in &destinations {
        if let Some(parent) = path.parent() {
            fieldcraft_io::fs::create_dir_all(parent)?;
        }
        fieldcraft_io::fs::write_output(path, doc.as_bytes())?;
    }
    Ok(json!({"path":path,"files":destinations.iter().map(|(_,path,_)|path).collect::<Vec<_>>(),
        "bytes":bytes,"original":true,"model_edits_exported":false,"history_interpreted":false}))
}
