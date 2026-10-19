//! FieldCraft's files.
//!
//! A `.fieldcraft` file ([`write`], [`read`]; docs/formats/fieldcraft.md) is a ZIP package of four
//! parts, stored uncompressed so each can be read with any unzip tool:
//!
//! - `fieldcraft.json`: the manifest. The format tag (`fieldcraft/2`), the program that wrote
//!   it, the geometry the model is built on (the SolidCraft design and the imported 3MF meshes,
//!   each a path and the SHA-256 of its content, opened again when the model is), the view state
//!   ([`View`]: what the viewport hides), and the size and SHA-256 of every other part.
//! - `model.json`: the model definition (unit system, sets, materials, sections, reference
//!   points, couplings, boundary conditions, loads, steps) as pretty-printed JSON with a stable
//!   key order, so two versions of a model diff line by line.
//! - `mesh.bin`: nodes, elements and geometry tags in a compact little-endian binary form
//!   ([`mesh`]).
//! - `journal.json`: the commands that built the model, one per line: a script that replays.
//!
//! Files of every earlier format are read: their parts are brought to the current format by the
//! migration chain ([`migrate`]) before they are deserialised; a file of a later format is refused
//! with a message naming the FieldCraft that wrote it. Reading never writes.
//!
//! The same model always gives the same bytes. Reading treats the file as hostile: every size is
//! capped before anything is allocated, every part's checksum is verified, and the model is
//! rebuilt through the same checked methods commands use, so a file cannot hold what commands
//! would refuse. Reading and writing work on bytes (they build for the browser); only [`fs`],
//! [`recovery`] and [`roots`] touch the file system, which in the browser is [`memory`].
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod fs;
pub mod memory;
pub mod mesh;
pub mod migrate;
pub mod recovery;
pub mod roots;
pub mod stl;
pub mod vtk;

use std::collections::{BTreeMap, BTreeSet};
use std::io::{Cursor, Read, Write};

use fieldcraft_model::Model;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};

/// The format version of the files this version writes, and the latest it reads (it reads every
/// earlier one through [`migrate`]).
pub const FORMAT_VERSION: u32 = 2;

/// The format tag of the files this version writes: `fieldcraft/<FORMAT_VERSION>`.
pub const FORMAT: &str = "fieldcraft/2";

/// The file name extension.
pub const EXTENSION: &str = "fieldcraft";

/// Largest file read.
pub const MAX_FILE_BYTES: u64 = 8 << 30;

/// Largest manifest, model definition or journal part.
pub const MAX_TEXT_PART_BYTES: u64 = 256 << 20;

/// Most journal entries read.
pub const MAX_JOURNAL_ENTRIES: usize = 1_000_000;

/// Most imported meshes a model refers to.
pub const MAX_MESH_SOURCES: usize = 10_000;

/// Longest geometry path, command id.
pub const MAX_PATH_LEN: usize = 4096;
const MAX_COMMAND_LEN: usize = 256;

const MANIFEST: &str = "fieldcraft.json";
const MODEL: &str = "model.json";
const MESH: &str = "mesh.bin";
const JOURNAL: &str = "journal.json";

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    /// Not a `.fieldcraft` file, or a damaged or invalid one.
    #[error("{0}")]
    Invalid(String),
    /// A `.fieldcraft` file of a format version this program does not read.
    #[error("{0}")]
    Version(String),
    #[error("limit reached: {0}")]
    Limit(String),
    /// The file system refused.
    #[error("{0}")]
    Io(String),
    /// The path lies outside the folders file access is confined to ([`roots::Roots`]).
    #[error("{0}")]
    Refused(String),
}

pub type Result<T> = std::result::Result<T, Error>;

fn invalid<T>(message: impl Into<String>) -> Result<T> {
    Err(Error::Invalid(message.into()))
}

/// A model error met while reading `part`.
pub(crate) fn model_error(part: &str, e: fieldcraft_model::Error) -> Error {
    match e {
        fieldcraft_model::Error::Limit(m) => Error::Limit(format!("{part}: {m}")),
        e => Error::Invalid(format!("{part}: {e}")),
    }
}

/// The lower-case hexadecimal SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes).iter().map(|b| format!("{b:02x}")).collect()
}

/// A file the model is built on: its path and the SHA-256 of the content that was opened.
///
/// A relative path is relative to the folder of the `.fieldcraft` file that refers to it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "SourceFileRepr", into = "SourceFileRepr")]
pub struct SourceFile {
    path: String,
    sha256: String,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SourceFileRepr {
    path: String,
    sha256: String,
}

impl TryFrom<SourceFileRepr> for SourceFile {
    type Error = Error;
    fn try_from(r: SourceFileRepr) -> Result<SourceFile> {
        SourceFile::new(r.path, r.sha256)
    }
}

impl From<SourceFile> for SourceFileRepr {
    fn from(s: SourceFile) -> Self {
        SourceFileRepr { path: s.path, sha256: s.sha256 }
    }
}

impl SourceFile {
    /// A path (1 to [`MAX_PATH_LEN`] bytes, no control characters) and a lower-case hexadecimal
    /// SHA-256.
    pub fn new(path: impl Into<String>, sha256: impl Into<String>) -> Result<SourceFile> {
        let (path, sha256) = (path.into(), sha256.into());
        if path.is_empty() || path.len() > MAX_PATH_LEN || path.chars().any(char::is_control) {
            return invalid(format!("a geometry path must be 1 to {MAX_PATH_LEN} bytes with no control characters"));
        }
        if sha256.len() != 64 || !sha256.bytes().all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)) {
            return invalid(format!("`{}` is not a SHA-256 (64 lower-case hexadecimal digits)", sha256.chars().take(80).collect::<String>()));
        }
        Ok(SourceFile { path, sha256 })
    }

    /// The file `path`, with the hash of `bytes` (what was read from it).
    pub fn of(path: impl Into<String>, bytes: &[u8]) -> Result<SourceFile> {
        SourceFile::new(path, sha256_hex(bytes))
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    pub fn sha256(&self) -> &str {
        &self.sha256
    }

    /// The same content at another path.
    pub fn with_path(&self, path: impl Into<String>) -> Result<SourceFile> {
        SourceFile::new(path, self.sha256.clone())
    }
}

/// The geometry a model is built on.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Geometry {
    /// The SolidCraft design.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub design: Option<SourceFile>,
    /// Face-tagged 3MF meshes, in the order they were imported (their surfaces follow each other
    /// in that order).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub meshes: Vec<SourceFile>,
}

/// Most items a view hides.
pub const MAX_HIDDEN: usize = 100_000;

/// Longest name of a hidden item (a body's name is a geometry name).
const MAX_HIDDEN_NAME: usize = fieldcraft_model::geometry::MAX_GEOMETRY_NAME_LEN;

/// The kinds of items the viewport draws one by one, and so can hide.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ItemKind {
    /// A body of the open design (or a cell of a partitioned one).
    Body,
    /// An imported surface.
    Surface,
    Datum,
    BoundaryCondition,
    Load,
    ReferencePoint,
    Coupling,
}

impl ItemKind {
    pub const ALL: [ItemKind; 7] = [
        ItemKind::Body,
        ItemKind::Surface,
        ItemKind::Datum,
        ItemKind::BoundaryCondition,
        ItemKind::Load,
        ItemKind::ReferencePoint,
        ItemKind::Coupling,
    ];
}

/// An item of the viewport, by kind and name.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Item {
    pub kind: ItemKind,
    pub name: String,
}

/// View state saved with a model (not part of its definition): the items the viewport hides.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct View {
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub hidden: BTreeSet<Item>,
}

impl View {
    pub fn is_empty(&self) -> bool {
        self.hidden.is_empty()
    }

    pub fn is_hidden(&self, kind: ItemKind, name: &str) -> bool {
        self.hidden.contains(&Item { kind, name: name.to_string() })
    }

    fn check(&self) -> Result<()> {
        if self.hidden.len() > MAX_HIDDEN {
            return Err(Error::Limit(format!("a view hides at most {MAX_HIDDEN} items")));
        }
        if self.hidden.iter().any(|i| i.name.is_empty() || i.name.len() > MAX_HIDDEN_NAME) {
            return invalid(format!("{MANIFEST}: a hidden item's name is 1 to {MAX_HIDDEN_NAME} bytes"));
        }
        Ok(())
    }
}

/// One journal entry: a command id and its parameters.
pub type JournalEntry = (String, Value);

/// What a `.fieldcraft` file holds.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelFile {
    pub model: Model,
    pub geometry: Geometry,
    pub journal: Vec<JournalEntry>,
    pub view: View,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    format: String,
    generator: String,
    #[serde(default)]
    geometry: Geometry,
    #[serde(default, skip_serializing_if = "View::is_empty")]
    view: View,
    parts: BTreeMap<String, Part>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Part {
    bytes: u64,
    sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct JournalLine {
    command: String,
    #[serde(default)]
    params: Value,
}

/// A journal entry as written (borrowed, so saving copies no parameters).
#[derive(Serialize)]
struct JournalLineRef<'a> {
    command: &'a str,
    params: &'a Value,
}

fn json_error(part: &str, e: serde_json::Error) -> Error {
    Error::Invalid(format!("{part}: {e}"))
}

/// The journal as text: one entry per line, so a diff shows the commands added.
fn journal_text(journal: &[JournalEntry]) -> Result<Vec<u8>> {
    let mut out = b"[".to_vec();
    for (i, (command, params)) in journal.iter().enumerate() {
        out.extend_from_slice(if i == 0 { b"\n  " } else { b",\n  " });
        serde_json::to_writer(&mut out, &JournalLineRef { command, params }).map_err(|e| json_error(JOURNAL, e))?;
    }
    out.extend_from_slice(if journal.is_empty() { b"]\n" } else { b"\n]\n" });
    Ok(out)
}

/// The bytes of a `.fieldcraft` file holding `model`, built on `geometry` by the commands in
/// `journal`, shown with `view`. The same arguments always give the same bytes.
pub fn write(model: &Model, geometry: &Geometry, journal: &[JournalEntry], view: &View) -> Result<Vec<u8>> {
    if geometry.meshes.len() > MAX_MESH_SOURCES {
        return Err(Error::Limit(format!("a model refers to at most {MAX_MESH_SOURCES} imported meshes")));
    }
    view.check()?;
    let mut model_json = serde_json::to_vec_pretty(&model.without_mesh()).map_err(|e| json_error(MODEL, e))?;
    model_json.push(b'\n');
    let mesh = mesh::encode(model)?;
    let journal = journal_text(journal)?;
    let parts = [(MODEL, model_json), (MESH, mesh), (JOURNAL, journal)];
    let manifest = Manifest {
        format: FORMAT.into(),
        generator: format!("FieldCraft {}", env!("CARGO_PKG_VERSION")),
        geometry: geometry.clone(),
        view: view.clone(),
        parts: parts.iter().map(|(name, bytes)| (name.to_string(), Part { bytes: bytes.len() as u64, sha256: sha256_hex(bytes) })).collect(),
    };
    let mut manifest_json = serde_json::to_vec_pretty(&manifest).map_err(|e| json_error(MANIFEST, e))?;
    manifest_json.push(b'\n');

    let total: usize = parts.iter().map(|(_, b)| b.len()).sum::<usize>() + manifest_json.len();
    let mut zip = zip::ZipWriter::new(Cursor::new(Vec::with_capacity(total.saturating_add(1024))));
    let zip_error = |e: zip::result::ZipError| Error::Io(format!("cannot build the package: {e}"));
    for (name, bytes) in std::iter::once((MANIFEST, &manifest_json)).chain(parts.iter().map(|(n, b)| (*n, b))) {
        // Fixed time and permissions: the bytes depend on the model only.
        let options = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Stored)
            .last_modified_time(zip::DateTime::default())
            .unix_permissions(0o644)
            .large_file(bytes.len() as u64 >= u64::from(u32::MAX));
        zip.start_file(name, options).map_err(zip_error)?;
        zip.write_all(bytes).map_err(|e| Error::Io(format!("cannot build the package: {e}")))?;
    }
    Ok(zip.finish().map_err(zip_error)?.into_inner())
}

/// Read a part of at most `cap` bytes, checking its size and ZIP checksum.
fn read_part(zip: &mut zip::ZipArchive<Cursor<&[u8]>>, name: &str, cap: u64, file_len: usize) -> Result<Vec<u8>> {
    let mut f = zip.by_name(name).map_err(|e| Error::Invalid(format!("no `{name}` part: {e}")))?;
    let size = f.size();
    if size > cap {
        return Err(Error::Limit(format!("part `{name}` is {size} bytes (at most {cap})")));
    }
    // Declared sizes are not trusted for allocation beyond the input itself.
    let mut buf = Vec::with_capacity(usize::try_from(size).unwrap_or(usize::MAX).min(file_len));
    (&mut f).take(cap.saturating_add(1)).read_to_end(&mut buf).map_err(|e| Error::Invalid(format!("part `{name}` is damaged: {e}")))?;
    if buf.len() as u64 != size {
        return invalid(format!("part `{name}` is damaged: {} bytes where {size} were declared", buf.len()));
    }
    Ok(buf)
}

/// Read a `.fieldcraft` file of any format up to [`FORMAT_VERSION`].
pub fn read(bytes: &[u8]) -> Result<ModelFile> {
    read_versioned(bytes).map(|(file, _)| file)
}

/// Read a `.fieldcraft` file, with the format it was written in: one older than
/// [`FORMAT_VERSION`] was migrated, and saving it again writes the current format.
pub fn read_versioned(bytes: &[u8]) -> Result<(ModelFile, u32)> {
    if bytes.len() as u64 > MAX_FILE_BYTES {
        return Err(Error::Limit(format!("the file is {} MB (at most {} MB)", bytes.len() >> 20, MAX_FILE_BYTES >> 20)));
    }
    if !bytes.starts_with(b"PK\x03\x04") {
        return invalid("not a FieldCraft file (not a ZIP package)");
    }
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).map_err(|e| Error::Invalid(format!("not a FieldCraft file: {e}")))?;

    // The version first: a later format may have parts and fields this version does not know,
    // and its reader is told to update FieldCraft, not that a field is unknown.
    let manifest_bytes = read_part(&mut zip, MANIFEST, MAX_TEXT_PART_BYTES, bytes.len())?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes).map_err(|e| json_error(MANIFEST, e))?;
    let format = migrate::version(&manifest)?;

    let names = [MANIFEST, MODEL, MESH, JOURNAL];
    let mut seen = Vec::with_capacity(zip.len());
    for name in zip.file_names() {
        let name = name.map_err(|e| Error::Invalid(format!("not a FieldCraft file: {e}")))?;
        if !names.contains(&name.as_ref()) {
            return invalid(format!("unexpected part `{}`", name.chars().take(80).collect::<String>()));
        }
        if seen.contains(&name) {
            return invalid(format!("part `{name}` appears twice"));
        }
        seen.push(name);
    }
    // The package (its parts' sizes and hashes) is checked as it is on disk, before migration.
    let listed: BTreeMap<String, Part> =
        serde_json::from_value(manifest.get("parts").cloned().unwrap_or(Value::Null)).map_err(|e| json_error(MANIFEST, e))?;
    if listed.len() != 3 {
        return invalid(format!("`{MANIFEST}` must list exactly the parts {MODEL}, {MESH} and {JOURNAL}"));
    }
    let mut part = |name: &str, cap: u64| -> Result<Vec<u8>> {
        let entry = listed.get(name).ok_or_else(|| Error::Invalid(format!("`{MANIFEST}` does not list `{name}`")))?;
        let data = read_part(&mut zip, name, cap, bytes.len())?;
        if data.len() as u64 != entry.bytes || sha256_hex(&data) != entry.sha256 {
            return invalid(format!("part `{name}` is damaged: its size or SHA-256 differs from the manifest"));
        }
        Ok(data)
    };
    let model_json = part(MODEL, MAX_TEXT_PART_BYTES)?;
    let mesh = part(MESH, mesh::MAX_BYTES)?;
    let journal_json = part(JOURNAL, MAX_TEXT_PART_BYTES)?;
    // A current file is deserialised straight from its bytes (a duplicate key is an error there,
    // and no copy of the parts is made); an older one goes through plain JSON for the steps.
    let (manifest, mut model, journal, mesh_bin): (Manifest, Model, Vec<JournalLine>, Vec<u8>) = if format == FORMAT_VERSION {
        drop(manifest);
        (
            serde_json::from_slice(&manifest_bytes).map_err(|e| json_error(MANIFEST, e))?,
            serde_json::from_slice(&model_json).map_err(|e| json_error(MODEL, e))?,
            serde_json::from_slice(&journal_json).map_err(|e| json_error(JOURNAL, e))?,
            mesh,
        )
    } else {
        let mut parts = migrate::Parts {
            manifest,
            model: serde_json::from_slice(&model_json).map_err(|e| json_error(MODEL, e))?,
            mesh,
            journal: serde_json::from_slice(&journal_json).map_err(|e| json_error(JOURNAL, e))?,
        };
        migrate::migrate(&mut parts)?;
        (
            serde_json::from_value(parts.manifest).map_err(|e| json_error(MANIFEST, e))?,
            serde_json::from_value(parts.model).map_err(|e| json_error(MODEL, e))?,
            serde_json::from_value(parts.journal).map_err(|e| json_error(JOURNAL, e))?,
            parts.mesh,
        )
    };
    drop((manifest_bytes, model_json, journal_json));
    if manifest.geometry.meshes.len() > MAX_MESH_SOURCES {
        return Err(Error::Limit(format!("a model refers to at most {MAX_MESH_SOURCES} imported meshes")));
    }
    manifest.view.check()?;
    if !model.nodes().is_empty() || !model.elements().is_empty() || !model.classification().entities().is_empty() {
        return invalid(format!("`{MODEL}` holds mesh data; the mesh belongs in `{MESH}`"));
    }
    mesh::decode_into(&mut model, &mesh_bin)?;
    drop(mesh_bin);
    let journal = read_journal(journal)?;
    Ok((ModelFile { model, geometry: manifest.geometry, journal, view: manifest.view }, format))
}

/// The format version of the `.fieldcraft` file `reader` holds, from its manifest alone (nothing
/// else is read or checked); `None` when it is not a FieldCraft file with a format tag. A later
/// format than this version reads gives its number.
pub fn format_of<R: Read + std::io::Seek>(reader: R) -> Option<u32> {
    let mut zip = zip::ZipArchive::new(reader).ok()?;
    let f = zip.by_name(MANIFEST).ok()?;
    let mut buf = Vec::new();
    f.take(MAX_TEXT_PART_BYTES).read_to_end(&mut buf).ok()?;
    let manifest: Value = serde_json::from_slice(&buf).ok()?;
    match migrate::version(&manifest) {
        Ok(n) => Some(n),
        Err(Error::Version(_)) => manifest.get("format")?.as_str()?.strip_prefix("fieldcraft/")?.parse().ok().or(Some(u32::MAX)),
        Err(_) => None,
    }
}

fn read_journal(lines: Vec<JournalLine>) -> Result<Vec<JournalEntry>> {
    if lines.len() > MAX_JOURNAL_ENTRIES {
        return Err(Error::Limit(format!("a journal holds at most {MAX_JOURNAL_ENTRIES} commands")));
    }
    lines
        .into_iter()
        .enumerate()
        .map(|(i, line)| {
            if line.command.is_empty() || line.command.len() > MAX_COMMAND_LEN || line.command.chars().any(char::is_control) {
                return invalid(format!("{JOURNAL}: entry {}: not a command id", i + 1));
            }
            Ok((line.command, line.params))
        })
        .collect()
}
