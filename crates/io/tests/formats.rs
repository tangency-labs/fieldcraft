//! Files of every released format open, migrate to what their own writer's reader gave, and
//! save again at the current format (crates/io/tests/formats/v<N>/README.md).

use std::io::Cursor;
use std::path::{Path, PathBuf};

use fieldcraft_io::migrate::{self, Parts};
use fieldcraft_io::{FORMAT, FORMAT_VERSION, ModelFile, format_of, read_versioned, write};
use serde_json::{Value, json};

/// A fixture: the program that wrote it (its folder), the file, and what that program's own
/// reader gave for it (`<name>.read.json`).
struct Fixture {
    writer: String,
    path: PathBuf,
    recorded: Value,
}

impl Fixture {
    fn name(&self) -> String {
        format!("{}/{}", self.writer, self.path.file_name().unwrap().to_string_lossy())
    }
}

/// Every fixture of format `n`, from every writer's folder.
fn fixtures(n: u32) -> Vec<Fixture> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats").join(format!("v{n}"));
    let mut out = Vec::new();
    for writer in std::fs::read_dir(&root).unwrap().map(|e| e.unwrap().path()).filter(|p| p.is_dir()) {
        for path in std::fs::read_dir(&writer).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|x| x == "fieldcraft")) {
            let recorded = serde_json::from_slice(&std::fs::read(path.with_extension("read.json")).unwrap()).unwrap();
            out.push(Fixture { writer: writer.file_name().unwrap().to_string_lossy().into_owned(), path, recorded });
        }
    }
    out.sort_by(|a, b| a.path.cmp(&b.path));
    out
}

/// What the current reader must give for a file of format `from` whose writer's reader gave
/// `recorded`: the steps from `from` to the current format, applied to the recorded JSON.
/// Format 1 → 2 changes nothing. A record from a reader that had no view state (0.4.0, 0.5.0)
/// holds `null`: the model has no hidden items.
fn expected(mut recorded: Value, from: u32) -> Value {
    assert!((1..=FORMAT_VERSION).contains(&from));
    if recorded["view"].is_null() {
        recorded["view"] = json!({});
    }
    recorded
}

/// The current reader's result in the record's shape, as raw JSON.
fn as_recorded(file: &ModelFile) -> Value {
    json!({
        "model": serde_json::to_value(&file.model).unwrap(),
        "geometry": serde_json::to_value(&file.geometry).unwrap(),
        "journal": file.journal.iter().map(|(c, p)| json!({"command": c, "params": p})).collect::<Vec<_>>(),
        "view": serde_json::to_value(&file.view).unwrap(),
    })
}

#[test]
fn every_format_1_fixture_opens_as_its_writer_read_it() {
    let all = fixtures(1);
    let writers: std::collections::BTreeSet<&str> = all.iter().map(|f| f.writer.as_str()).collect();
    assert_eq!(writers.into_iter().collect::<Vec<_>>(), ["0.4.0", "0.5.0", "0.6.0", "0.7.0", "develop-2eb3515c"]);
    for f in &all {
        let name = f.name();
        let bytes = std::fs::read(&f.path).unwrap();
        assert_eq!(format_of(Cursor::new(&bytes)), Some(1), "{name}");
        let (file, format) = read_versioned(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(format, 1, "{name}");
        let (got, want) = (as_recorded(&file), expected(f.recorded.clone(), 1));
        for key in ["model", "geometry", "journal", "view"] {
            assert!(got[key] == want[key], "{name}: {key} differs from what its writer read");
        }
        assert_eq!(got, want, "{name}");
    }
}

/// The step 1 → 2 on every real format 1 file: only the tag changes; every part's content is
/// the same.
#[test]
fn the_step_from_1_to_2_changes_only_the_tag() {
    for f in fixtures(1) {
        let name = f.name();
        let ps = parts(&std::fs::read(&f.path).unwrap());
        let json = |part: &str| -> Value { serde_json::from_slice(&ps[part]).unwrap() };
        let before =
            Parts { manifest: json("fieldcraft.json"), model: json("model.json"), mesh: ps["mesh.bin"].clone(), journal: json("journal.json") };
        let mut after = before.clone();
        migrate::step(&mut after, 1).unwrap();
        assert_eq!(after.manifest["format"], "fieldcraft/2", "{name}");
        let mut manifest = after.manifest.clone();
        manifest["format"] = before.manifest["format"].clone();
        assert_eq!(
            (&manifest, &after.model, &after.mesh, &after.journal),
            (&before.manifest, &before.model, &before.mesh, &before.journal),
            "{name}"
        );
    }
}

/// Saving a migrated model writes the current format; reading that gives the same model, and
/// saving it again the same bytes.
#[test]
fn migrated_fixtures_save_and_reopen_bit_identical() {
    for f in fixtures(1) {
        let name = f.name();
        let (file, _) = read_versioned(&std::fs::read(&f.path).unwrap()).unwrap();
        let saved = write(&file.model, &file.geometry, &file.journal, &file.view).unwrap();
        assert_eq!(format_of(Cursor::new(&saved)), Some(FORMAT_VERSION), "{name}");
        let (again, format) = read_versioned(&saved).unwrap();
        assert_eq!(format, FORMAT_VERSION, "{name}");
        assert_eq!(again, file, "{name}");
        assert_eq!(write(&again.model, &again.geometry, &again.journal, &again.view).unwrap(), saved, "{name}: bytes");
        let manifest: Value = serde_json::from_slice(&parts(&saved)["fieldcraft.json"]).unwrap();
        assert_eq!(manifest["format"], FORMAT);
    }
}

fn parts(bytes: &[u8]) -> std::collections::BTreeMap<String, Vec<u8>> {
    use std::io::Read;
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).unwrap();
    (0..zip.len())
        .map(|i| {
            let mut f = zip.by_index(i).unwrap();
            let mut b = Vec::new();
            f.read_to_end(&mut b).unwrap();
            (f.name().unwrap().to_string(), b)
        })
        .collect()
}

/// What `format_of` makes of things that are not a FieldCraft file, or are a later one.
#[test]
fn format_of_reads_only_the_tag() {
    assert_eq!(format_of(Cursor::new(b"not a zip")), None);
    let (file, _) = read_versioned(&std::fs::read(&fixtures(1)[0].path).unwrap()).unwrap();
    let saved = write(&file.model, &file.geometry, &file.journal, &file.view).unwrap();
    let mut ps = parts(&saved);
    let mut manifest: Value = serde_json::from_slice(&ps["fieldcraft.json"]).unwrap();
    for (tag, expect) in [("fieldcraft/9", Some(9)), ("fieldcraft/99999999999", Some(u32::MAX)), ("other/1", None)] {
        manifest["format"] = Value::from(tag);
        ps.insert("fieldcraft.json".into(), serde_json::to_vec(&manifest).unwrap());
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for (name, b) in &ps {
            zip.start_file(name.as_str(), zip::write::SimpleFileOptions::default()).unwrap();
            std::io::Write::write_all(&mut zip, b).unwrap();
        }
        let bytes = zip.finish().unwrap().into_inner();
        assert_eq!(format_of(Cursor::new(&bytes)), expect, "{tag}");
    }
}
