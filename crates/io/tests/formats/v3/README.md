# Format 3 fixtures

Real `.fieldcraft` files of format `fieldcraft/3`, one folder per release that wrote them. Format 3
was written by FieldCraft 0.8.0 and 0.8.1 (format 4, #1049, follows):

| Folder | Written by | Commit |
|---|---|---|
| `0.8.0/` | FieldCraft 0.8.0 (release tag `v0.8.0`, its `fieldcraft-cli-0.8.0-macos-aarch64` release binary) | the tag |
| `0.8.1/` | FieldCraft 0.8.1 (`v0.8.1`, `fieldcraft-cli-0.8.1-macos-aarch64`) | the tag |

`crates/io/tests/formats.rs` opens every one with the current reader and checks that it gives what
its writer's own reader gave (step 3 → 4 changes only the tag); the engine's tests open them
through `file.open`, which reports `migrated_from: 3`.

**Never edit, regenerate or replace these files.** They stand for the files users have on disk.

| Fixture | What it holds |
|---|---|
| `design-mesh-results` | a tet4 mesh of the design `hollow-box.solidcraft` (a copy of `crates/geometry/tests/designs/hollow-box.solidcraft`) beside an imported untagged 3MF surface `probe.3mf`, with a set, a material, a section and a frequency step, solved: the file keeps the step's results (`results.bin`) and embeds both geometry files (`geometry/hollow-box.solidcraft`, `geometry/probe.3mf`), the parts only format 3 and later list |

## How they were made

`design-mesh-results.recipe.json` is the script that wrote `design-mesh-results.fieldcraft`, run in
each folder with that release's `fieldcraft-cli run design-mesh-results.recipe.json` (the SHA-256
of each release binary checked against the release's `SHA256SUMS.txt`). `probe.3mf` is a unit
tetrahedron, 10 mm a side at x = 200 mm, written for these fixtures (a stored ZIP holding
`3D/3dmodel.model` alone).

`design-mesh-results.read.json` is what that release's reader (`fieldcraft_io::read` at the tag)
gave for the file, as one line of JSON: the model (mesh included), the geometry, the journal with
its header, the view, the number of results and the embedded parts. It was written by this test,
added to the tag's `crates/io/tests/` and run once with `FX` set to the folder:

```rust
//! Records what this release's reader gives for every `.fieldcraft` file in `$FX`.
#[test]
fn dump() {
    let dir = std::path::PathBuf::from(std::env::var("FX").unwrap());
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        if p.extension().and_then(|x| x.to_str()) != Some("fieldcraft") { continue; }
        let f = fieldcraft_io::read(&std::fs::read(&p).unwrap()).unwrap();
        let v = serde_json::json!({
            "model": serde_json::to_value(&f.model).unwrap(),
            "geometry": serde_json::to_value(&f.geometry).unwrap(),
            "journal": {
                "units": serde_json::to_value(&f.journal_header.units).unwrap(),
                "model_units": serde_json::to_value(&f.journal_header.model_units).unwrap(),
                "commands": f.journal.iter().map(|(c, p)| serde_json::json!({"command": c, "params": p})).collect::<Vec<_>>(),
            },
            "view": serde_json::to_value(&f.view).unwrap(),
            "results": f.results.len(),
            "embedded": f.embedded.keys().collect::<Vec<_>>(),
        });
        let mut text = serde_json::to_vec(&v).unwrap();
        text.push(b'\n');
        std::fs::write(p.with_extension("read.json"), text).unwrap();
    }
}
```
