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

/// A fresh folder holding a copy of the format 1 fixture `name` as `bracket.fieldcraft`.
fn folder_with_v1(tag: &str, name: &str) -> (PathBuf, PathBuf) {
    let dir = std::env::temp_dir().join(format!("fieldcraft-formats-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let model = dir.join("bracket.fieldcraft");
    std::fs::copy(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats/v1/0.7.0").join(name), &model).unwrap();
    (dir, model)
}

fn listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().into_owned()).collect();
    names.sort();
    names
}

#[cfg(not(target_arch = "wasm32"))]
#[test]
fn saving_over_an_older_format_keeps_a_copy_once() {
    use fieldcraft_io::fs::{backup_path, format_of_file, save_model};
    let (dir, model) = folder_with_v1("backup", "hex-static.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    let (file, _) = read_versioned(&original).unwrap();
    let bytes = write(&file.model, &file.geometry, &file.journal, &file.view).unwrap();
    assert_eq!(format_of_file(&model), Some(1));

    let kept = save_model(&model, &bytes).unwrap().expect("a copy of the format 1 file");
    assert_eq!(kept, dir.join("bracket.v1.fieldcraft"));
    assert_eq!(std::fs::read(&kept).unwrap(), original, "the copy is the original, byte for byte");
    assert_eq!(std::fs::read(&model).unwrap(), bytes);
    assert_eq!(format_of_file(&model), Some(FORMAT_VERSION));

    // The next save replaces `.bak` as always, but makes no copy and leaves the first one.
    assert_eq!(save_model(&model, &bytes).unwrap(), None);
    assert_eq!(std::fs::read(&kept).unwrap(), original);
    assert_eq!(std::fs::read(backup_path(&model)).unwrap(), bytes);
    assert_eq!(listing(&dir), ["bracket.fieldcraft", "bracket.fieldcraft.bak", "bracket.v1.fieldcraft"]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// A name already taken (an earlier copy, or a file of the user's) is never replaced: the copy
/// takes the next free name.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn a_format_backup_never_replaces_a_file() {
    use fieldcraft_io::fs::save_model;
    let (dir, model) = folder_with_v1("collision", "every-kind.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    std::fs::write(dir.join("bracket.v1.fieldcraft"), b"the user's own file").unwrap();
    std::fs::write(dir.join("bracket.v1-2.fieldcraft"), b"another").unwrap();
    let (file, _) = read_versioned(&original).unwrap();
    let bytes = write(&file.model, &file.geometry, &file.journal, &file.view).unwrap();
    let kept = save_model(&model, &bytes).unwrap().unwrap();
    assert_eq!(kept, dir.join("bracket.v1-3.fieldcraft"));
    assert_eq!(std::fs::read(&kept).unwrap(), original);
    assert_eq!(std::fs::read(dir.join("bracket.v1.fieldcraft")).unwrap(), b"the user's own file");
    assert_eq!(std::fs::read(dir.join("bracket.v1-2.fieldcraft")).unwrap(), b"another");

    // A later format's file is kept the same way before this version writes over it; a file
    // that is not a model, or a current one, needs no copy.
    let later = dir.join("later.fieldcraft");
    let mut ps = parts(&bytes);
    let mut manifest: Value = serde_json::from_slice(&ps["fieldcraft.json"]).unwrap();
    manifest["format"] = Value::from(format!("fieldcraft/{}", FORMAT_VERSION + 1));
    ps.insert("fieldcraft.json".into(), serde_json::to_vec(&manifest).unwrap());
    let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for (name, b) in &ps {
        zip.start_file(name.as_str(), zip::write::SimpleFileOptions::default()).unwrap();
        std::io::Write::write_all(&mut zip, b).unwrap();
    }
    std::fs::write(&later, zip.finish().unwrap().into_inner()).unwrap();
    assert_eq!(save_model(&later, &bytes).unwrap(), Some(dir.join(format!("later.v{}.fieldcraft", FORMAT_VERSION + 1))));
    // Anything that is not positively a current model is kept: a file that is not a model at
    // all (a Save As over it), a package whose manifest cannot be read, a malformed tag.
    let plain = dir.join("notes.txt");
    std::fs::write(&plain, b"the user's notes").unwrap();
    let kept = save_model(&plain, &bytes).unwrap().unwrap();
    assert_eq!((kept.file_name().unwrap().to_str().unwrap(), std::fs::read(&kept).unwrap()), ("notes.backup.txt", b"the user's notes".to_vec()));
    for (name, manifest) in
        [("unreadable", b"{not json".to_vec()), ("malformed", br#"{"format": "fieldcraft/01"}"#.to_vec()), ("untagged", b"{}".to_vec())]
    {
        let mut ps = parts(&bytes);
        ps.insert("fieldcraft.json".into(), manifest);
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for (part, b) in &ps {
            zip.start_file(part.as_str(), zip::write::SimpleFileOptions::default()).unwrap();
            std::io::Write::write_all(&mut zip, b).unwrap();
        }
        let odd = zip.finish().unwrap().into_inner();
        let path = dir.join(format!("{name}.fieldcraft"));
        std::fs::write(&path, &odd).unwrap();
        let kept = save_model(&path, &bytes).unwrap().unwrap_or_else(|| panic!("{name}: no copy"));
        assert_eq!(kept, dir.join(format!("{name}.backup.fieldcraft")));
        assert_eq!(std::fs::read(&kept).unwrap(), odd, "{name}");
    }
    // A current model, or no file at all, needs no copy.
    assert_eq!(save_model(&dir.join("new.fieldcraft"), &bytes).unwrap(), None);
    assert_eq!(save_model(&dir.join("new.fieldcraft"), &bytes).unwrap(), None);
    let _ = std::fs::remove_dir_all(&dir);
}

/// When the copy cannot be made, nothing is saved: the original stays as it was.
#[cfg(unix)]
#[test]
fn no_copy_no_save() {
    use std::os::unix::fs::PermissionsExt;
    let (dir, model) = folder_with_v1("readonly", "hex-static.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o555)).unwrap();
    let r = fieldcraft_io::fs::save_model(&model, b"new bytes");
    std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o755)).unwrap();
    // Running as root writes anyway; then the copy was made and the save went through.
    if let Err(e) = r {
        assert!(e.to_string().contains("could not be kept"), "{e}");
        assert_eq!(std::fs::read(&model).unwrap(), original);
        assert_eq!(listing(&dir), ["bracket.fieldcraft"]);
    }
    let _ = std::fs::remove_dir_all(&dir);
}
