//! `.fieldcraft` files: round trips, stable bytes and hostile files.

use std::io::{Cursor, Read, Write};

use fieldcraft_io::{Error, Geometry, Item, ItemKind, ModelFile, SourceFile, View, mesh, read, write};
use fieldcraft_model::{DisplaySystem, ElementKind, Entity, GeometryRef, Model};
use serde_json::{Value, json};

/// A hex and a tet sharing a face, tagged on geometry, with a few definitions.
fn sample() -> ModelFile {
    let mut m = Model::default();
    let corners = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.], [0.5, 0.5, -1.0]];
    for (i, c) in corners.into_iter().enumerate() {
        m.set_node(i as u32 + 1, c).unwrap();
    }
    m.set_node(42, [-0.0, 1e-300, 123456.789]).unwrap();
    m.set_element(1, ElementKind::Hex8, (1..=8).collect()).unwrap();
    m.set_element(7, ElementKind::Tet4, vec![1, 3, 2, 9]).unwrap();
    let body = m.geometry_tag(&GeometryRef::body("Body1").unwrap()).unwrap();
    let top = m.geometry_tag(&GeometryRef::face("Body1", "F2:top").unwrap()).unwrap();
    let edge = m.geometry_tag(&GeometryRef::new("Body1", Entity::Edge("E1".into())).unwrap()).unwrap();
    let vertex = m.geometry_tag(&GeometryRef::new("Body1", Entity::Vertex("V1".into())).unwrap()).unwrap();
    m.tag_element(1, body).unwrap();
    m.tag_element_face(1, 2, top).unwrap();
    m.tag_element_edge(1, 1, edge).unwrap();
    m.tag_element(7, body).unwrap();
    m.tag_node(1, vertex).unwrap();
    m.set_display_units(DisplaySystem::from_label("mm-t-s-N-MPa"));
    let set = serde_json::from_value(json!({"geometry": [{"body": "Body1"}]})).unwrap();
    m.add_set("all", set).unwrap();
    let steel = serde_json::from_value(json!({"elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}, "density": 7.85e-9})).unwrap();
    m.add_material("Steel", steel).unwrap();
    let geometry = Geometry {
        design: Some(SourceFile::of("parts/plate.solidcraft", b"design bytes").unwrap()),
        meshes: vec![SourceFile::of("/abs/tet.3mf", b"3mf bytes").unwrap()],
    };
    let journal =
        vec![("mesh.node.add".to_string(), json!({"label": 1, "xyz": [0, 0, 0]})), ("model.units".to_string(), json!({"system": "mm-t-s-N-MPa"}))];
    let view = View {
        hidden: [Item { kind: ItemKind::Datum, name: "Sym".into() }, Item { kind: ItemKind::Body, name: "Body1".into() }].into(),
        ..View::default()
    };
    ModelFile { model: m, geometry, journal, journal_header: Default::default(), view, ..ModelFile::default() }
}

fn bytes_of(f: &ModelFile) -> Vec<u8> {
    write(&f.model, &f.geometry, &f.journal, &f.journal_header, &f.view).unwrap()
}

/// The package's parts, by name.
fn parts(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
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

/// A package of these parts, written as given (no checksums fixed up).
fn package(parts: &[(String, Vec<u8>)]) -> Vec<u8> {
    let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for (name, bytes) in parts {
        zip.start_file(name, zip::write::SimpleFileOptions::default().compression_method(zip::CompressionMethod::Stored)).unwrap();
        zip.write_all(bytes).unwrap();
    }
    zip.finish().unwrap().into_inner()
}

/// Replace one part and fix the manifest's size and hash for it, so the change reaches the
/// part's own checks.
fn with_part(bytes: &[u8], name: &str, new: Vec<u8>) -> Vec<u8> {
    let mut ps = parts(bytes);
    let mut manifest: Value = serde_json::from_slice(&ps[0].1).unwrap();
    manifest["parts"][name] = json!({"bytes": new.len(), "sha256": fieldcraft_io::sha256_hex(&new)});
    ps[0].1 = serde_json::to_vec(&manifest).unwrap();
    for p in &mut ps {
        if p.0 == name {
            p.1 = new.clone();
        }
    }
    package(&ps)
}

fn with_manifest(bytes: &[u8], edit: impl FnOnce(&mut Value)) -> Vec<u8> {
    let mut ps = parts(bytes);
    let mut manifest: Value = serde_json::from_slice(&ps[0].1).unwrap();
    edit(&mut manifest);
    ps[0].1 = serde_json::to_vec(&manifest).unwrap();
    package(&ps)
}

#[test]
fn round_trip_keeps_everything_and_bytes_are_stable() {
    let f = sample();
    let bytes = bytes_of(&f);
    let back = read(&bytes).unwrap();
    assert_eq!(back, f);
    assert_eq!(back.model.node(42).unwrap().xyz[0].to_bits(), (-0.0f64).to_bits(), "signed zero survives");
    assert_eq!(bytes_of(&back), bytes, "the same model gives the same bytes");

    // The parts, in order, all stored as plain bytes.
    let names: Vec<String> = parts(&bytes).into_iter().map(|(n, _)| n).collect();
    assert_eq!(names, ["fieldcraft.json", "model.json", "mesh.bin", "journal.json"]);
    let ps = parts(&bytes);
    let manifest: Value = serde_json::from_slice(&ps[0].1).unwrap();
    assert_eq!(manifest["format"], fieldcraft_io::FORMAT);
    assert_eq!(manifest["geometry"]["design"]["path"], "parts/plate.solidcraft");
    let model_text = String::from_utf8(ps[1].1.clone()).unwrap();
    assert!(model_text.contains("\"Steel\"") && !model_text.contains("\"nodes\""), "{model_text}");
    let journal_text = String::from_utf8(ps[3].1.clone()).unwrap();
    assert_eq!(journal_text.lines().count(), 4, "one command per line: {journal_text}");

    // An empty model too.
    let empty = ModelFile::default();
    assert_eq!(read(&bytes_of(&empty)).unwrap(), empty);
}

#[test]
fn a_definition_change_changes_only_its_lines() {
    let f = sample();
    let mut g = f.clone();
    let alu = serde_json::from_value(json!({"elastic": {"type": "isotropic", "young": 70000.0, "poisson": 0.33}})).unwrap();
    g.model.add_material("Aluminium", alu).unwrap();
    let (a, b) = (parts(&bytes_of(&f)), parts(&bytes_of(&g)));
    assert_eq!(a[2].1, b[2].1, "the mesh is unchanged");
    let (ta, tb) = (String::from_utf8(a[1].1.clone()).unwrap(), String::from_utf8(b[1].1.clone()).unwrap());
    let removed = ta.lines().filter(|l| !tb.lines().any(|m| m == *l)).count();
    assert_eq!(removed, 0, "adding a material only adds lines:\n{ta}\n---\n{tb}");
}

#[test]
fn hostile_files_are_refused_cleanly() {
    let good = bytes_of(&sample());
    let mesh_bin = parts(&good)[2].1.clone();

    // Truncated anywhere.
    for cut in (0..good.len()).step_by(7).chain([good.len() - 1]) {
        assert!(read(&good[..cut]).is_err(), "cut at {cut}");
    }
    // A flipped byte anywhere is caught by a checksum or a check.
    for at in (0..good.len()).step_by(11) {
        let mut b = good.clone();
        b[at] ^= 0x5a;
        if let Ok(f) = read(&b) {
            // Only bytes outside every part's content (ZIP metadata the reader ignores) may flip.
            assert_eq!(f, sample(), "byte {at}");
        }
    }

    // Not a model file at all.
    for junk in [&b""[..], b"PK", b"PK\x03\x04garbage", b"{\"format\": \"fieldcraft/1\"}", &[0u8; 4096]] {
        assert!(matches!(read(junk), Err(Error::Invalid(_))), "{junk:?}");
    }

    // Wrong and later versions. A later format is refused for its version even when it holds
    // fields and parts this version does not know (which would otherwise be the first error).
    let later = format!("fieldcraft/{}", fieldcraft_io::FORMAT_VERSION + 1);
    let newer = with_manifest(&good, |m| m["format"] = json!(later));
    assert!(matches!(read(&newer), Err(Error::Version(_))));
    let newer_fields = with_manifest(&good, |m| {
        m["format"] = json!(later);
        m["generator"] = json!("FieldCraft 9.0.0");
        m["colour"] = json!("red");
        m["parts"]["results.bin"] = json!({"bytes": 1, "sha256": "0".repeat(64)});
    });
    let mut ps = parts(&newer_fields);
    ps.push(("results.bin".into(), b"r".to_vec()));
    ps[1].1 = br#"{"units": "m-kg-s-N-Pa", "temperature_scale": "kelvin"}"#.to_vec();
    match read(&package(&ps)) {
        Err(Error::Version(m)) => {
            assert!(m.contains("made by FieldCraft 9.0.0") && m.contains(&format!("format {}", fieldcraft_io::FORMAT_VERSION + 1)), "{m}");
            assert!(m.contains("Update FieldCraft"), "{m}");
        }
        other => panic!("{other:?}"),
    }
    for f in [json!("solidcraft/16"), json!(null), json!(1), json!("fieldcraft/"), json!("fieldcraft/0"), json!("fieldcraft/two")] {
        let b = with_manifest(&good, |m| m["format"] = f.clone());
        assert!(matches!(read(&b), Err(Error::Invalid(_))), "{f}");
    }
    let extra = with_manifest(&good, |m| m["colour"] = json!("red"));
    assert!(read(&extra).is_err());

    // Bad checksums: a part that differs from its manifest entry.
    let mut ps = parts(&good);
    ps[2].1[60] ^= 1;
    assert!(read(&package(&ps)).unwrap_err().to_string().contains("damaged"));
    let bad_hash = with_manifest(&good, |m| m["parts"]["model.json"]["sha256"] = json!("0".repeat(64)));
    assert!(read(&bad_hash).unwrap_err().to_string().contains("damaged"));
    let bad_geometry_hash = with_manifest(&good, |m| m["geometry"]["design"]["sha256"] = json!("XYZ"));
    assert!(read(&bad_geometry_hash).is_err());
    let bad_path = with_manifest(&good, |m| m["geometry"]["design"]["path"] = json!("a\u{0}b"));
    assert!(read(&bad_path).is_err());

    // Missing, extra and repeated parts.
    let mut ps = parts(&good);
    ps.remove(3);
    assert!(read(&package(&ps)).is_err());
    let mut ps = parts(&good);
    ps.push(("../evil".into(), b"x".to_vec()));
    assert!(read(&package(&ps)).is_err());

    // Mesh parts: huge counts, NaN, out of order, unknown kinds, trailing bytes.
    let mut huge = mesh_bin.clone();
    huge[8..16].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(matches!(read(&with_part(&good, "mesh.bin", huge)), Err(Error::Limit(_))));
    let mut many = mesh_bin.clone();
    many[8..16].copy_from_slice(&40_000_000u64.to_le_bytes());
    assert!(read(&with_part(&good, "mesh.bin", many)).unwrap_err().to_string().contains("too short"));
    let first_x = 48 + 4;
    for bad in [f64::NAN, f64::INFINITY, -1e300] {
        let mut b = mesh_bin.clone();
        b[first_x..first_x + 8].copy_from_slice(&bad.to_le_bytes());
        assert!(read(&with_part(&good, "mesh.bin", b)).is_err(), "{bad}");
    }
    let mut swapped = mesh_bin.clone();
    swapped[48..52].copy_from_slice(&5u32.to_le_bytes());
    assert!(read(&with_part(&good, "mesh.bin", swapped)).unwrap_err().to_string().contains("order"));
    let first_kind = 48 + 10 * 28 + 4;
    let mut kind = mesh_bin.clone();
    kind[first_kind] = 99;
    assert!(read(&with_part(&good, "mesh.bin", kind)).unwrap_err().to_string().contains("unknown kind"));
    let mut trailing = mesh_bin.clone();
    trailing.push(0);
    assert!(read(&with_part(&good, "mesh.bin", trailing)).is_err());
    let mut magic = mesh_bin.clone();
    magic[7] = b'9';
    let e = read(&with_part(&good, "mesh.bin", magic)).unwrap_err();
    assert!(matches!(&e, Error::Version(m) if m.contains("mesh encoding `09`")), "{e}");
    let mut magic = mesh_bin.clone();
    magic[0] = b'X';
    assert!(matches!(read(&with_part(&good, "mesh.bin", magic)), Err(Error::Invalid(_))));

    // Model definitions that commands would refuse, and mesh data in the text part.
    for model in [
        json!({"materials": {"1bad": {}}}),
        json!({"sections": {"S": {"set": "nope", "material": "nope", "kind": {"type": "solid"}}}}),
        json!({"nodes": {"1": {"xyz": [0, 0, 0]}}}),
        json!({"display_units": "furlongs"}),
        json!({"units": "mm-t-s-N-MPa"}),
        json!({"display_units": {"preset": "mm-t-s-N-MPa", "overrides": {"length": "MPa"}}}),
    ] {
        let b = with_part(&good, "model.json", serde_json::to_vec(&model).unwrap());
        assert!(read(&b).is_err(), "{model}");
    }

    // A key given twice is an error in a current file, not the last one silently winning.
    let twice = with_part(&good, "model.json", br#"{"display_units": "mm-t-s-N-MPa", "display_units": "m-kg-s-N-Pa"}"#.to_vec());
    let e = read(&twice).unwrap_err();
    assert!(e.to_string().contains("duplicate field"), "{e}");
    let mut ps = parts(&good);
    let manifest = String::from_utf8(ps[0].1.clone()).unwrap();
    ps[0].1 = manifest.replacen("\"generator\"", "\"generator\": \"FieldCraft 0.0.1\",\n  \"generator\"", 1).into_bytes();
    assert!(read(&package(&ps)).unwrap_err().to_string().contains("duplicate field"));

    // Journals: not a list, entries without a command, unknown fields.
    for journal in [
        json!({}),
        json!([]),
        json!({"units": "mm-t-s-N-MPa", "commands": [{"params": {}}]}),
        json!({"units": "mm-t-s-N-MPa", "commands": [{"command": ""}]}),
        json!({"units": "mm-t-s-N-MPa", "commands": [{"command": "a.b", "when": 1}]}),
        json!({"units": "mm-t-s-N-MPa", "commands": [7]}),
        json!({"units": "furlongs", "commands": []}),
        json!({"commands": []}),
        json!({"units": "mm-t-s-N-MPa", "commands": [], "extra": 1}),
    ] {
        let b = with_part(&good, "journal.json", serde_json::to_vec(&journal).unwrap());
        assert!(read(&b).is_err(), "{journal}");
    }
    // A bare list (format 2's journal) in a format 3 file says what it lacks.
    let bare = with_part(&good, "journal.json", br#"[{"command": "model.check"}]"#.to_vec());
    let e = read(&bare).unwrap_err().to_string();
    assert!(e.contains("naming the units") && e.contains("bare list"), "{e}");
    let ok = with_part(&good, "journal.json", br#"{"units": "mm-t-s-N-MPa", "commands": [{"command": "model.check"}]}"#.to_vec());
    assert_eq!(read(&ok).unwrap().journal, vec![("model.check".to_string(), Value::Null)]);
}

#[test]
fn mesh_part_rejects_what_the_model_rejects() {
    let f = sample();
    let good = mesh::encode(&f.model).unwrap();
    // Decoding into a model that has a mesh already is refused.
    let mut m = f.model.clone();
    assert!(mesh::decode_into(&mut m, &good).is_err());
    // An element naming a missing node.
    let mut b = good.clone();
    let first_conn = 48 + 10 * 28 + 5;
    b[first_conn..first_conn + 4].copy_from_slice(&1000u32.to_le_bytes());
    assert!(mesh::decode_into(&mut Model::default(), &b).is_err());
    // Every prefix fails without panicking.
    for cut in 0..good.len() {
        assert!(mesh::decode_into(&mut Model::default(), &good[..cut]).is_err(), "{cut}");
    }
    let mut back = f.model.without_mesh();
    mesh::decode_into(&mut back, &good).unwrap();
    assert_eq!(back, f.model);
}

/// The view's hidden items come back as written; a hostile hidden list is refused.
#[test]
fn the_view_state_round_trips_and_is_checked() {
    let f = sample();
    let back = read(&bytes_of(&f)).unwrap();
    assert_eq!(back.view, f.view);
    assert!(back.view.is_hidden(ItemKind::Datum, "Sym") && !back.view.is_hidden(ItemKind::Load, "Sym"));
    let bytes = bytes_of(&f);
    for bad in [json!([{"kind": "body", "name": ""}]), json!([{"kind": "nope", "name": "A"}]), json!("x")] {
        let mut ps = parts(&bytes);
        let mut manifest: Value = serde_json::from_slice(&ps[0].1).unwrap();
        manifest["view"]["hidden"] = bad.clone();
        ps[0].1 = serde_json::to_vec(&manifest).unwrap();
        assert!(read(&package(&ps)).is_err(), "{bad}");
    }
    let empty = View::default();
    let plain = write(&f.model, &f.geometry, &f.journal, &f.journal_header, &empty).unwrap();
    let manifest: Value = serde_json::from_slice(&parts(&plain)[0].1).unwrap();
    assert!(manifest.get("view").is_none(), "an empty view writes nothing: {manifest}");
}

/// A small heat transfer result of step `Heat` on nodes 1 to 3.
fn heat_result() -> fieldcraft_solver::HeatResult {
    fieldcraft_solver::HeatResult {
        step: "Heat".into(),
        nodes: vec![1, 2, 3],
        times: vec![0.0, 0.5],
        frames: vec![vec![293.15; 3], vec![300.0, 310.0, f64::NAN]],
        flux: vec![[1.0, 0.0, -0.0]; 3],
        elements: Vec::new(),
        reaction: vec![0.0, -2.5, 0.0],
        loads: vec![fieldcraft_solver::heat::LoadFlow { name: "Film".into(), kind: "convection".into(), heat: 1.5 }],
        applied_total: 1.5,
        reaction_total: -2.5,
        stats: Default::default(),
        warnings: vec!["a note".into()],
    }
}

fn with_results(f: &ModelFile, results: &[(&str, fieldcraft_io::results::StoredRef<'_>)]) -> Vec<u8> {
    let extras = fieldcraft_io::Extras { results: results.to_vec(), ..Default::default() };
    fieldcraft_io::write_with(&f.model, &f.geometry, &f.journal, &f.journal_header, &f.view, &extras).unwrap()
}

/// The results part: kept bit for bit with its definition hash, listed in the manifest, checked
/// as hostile input; one of another results encoding leaves the file openable without it.
#[test]
fn results_round_trip_and_hostile_results_are_refused() {
    use fieldcraft_io::results::{self, Solved, SolvedRef, StoredRef};
    let f = sample();
    let heat = fieldcraft_solver::StepResult::Heat(heat_result());
    let hash = [7u8; 32];
    let bytes = with_results(&f, &[("Heat", StoredRef { definition: hash, result: SolvedRef::Implicit(&heat) })]);
    let back = read(&bytes).unwrap();
    assert_eq!(back.results.len(), 1);
    assert_eq!(back.results[0].definition, hash);
    let Solved::Implicit(fieldcraft_solver::StepResult::Heat(h)) = &back.results[0].result else { panic!() };
    assert!(h.frames[1][2].is_nan() && h.flux[0][2].to_bits() == (-0.0f64).to_bits(), "every bit kept");
    assert_eq!(h.loads, heat_result().loads);
    assert_eq!((h.nodes.clone(), h.times.clone(), h.warnings.clone()), (vec![1, 2, 3], vec![0.0, 0.5], vec!["a note".to_string()]));
    let names: Vec<String> = parts(&bytes).into_iter().map(|(n, _)| n).collect();
    assert_eq!(names, ["fieldcraft.json", "model.json", "mesh.bin", "journal.json", "results.bin"]);
    let manifest: Value = serde_json::from_slice(&parts(&bytes)[0].1).unwrap();
    let rb = parts(&bytes)[4].1.clone();
    assert_eq!(manifest["parts"]["results.bin"], json!({"bytes": rb.len(), "sha256": fieldcraft_io::sha256_hex(&rb)}));
    assert_eq!(&rb[..8], results::MAGIC);
    // The same results, the same bytes.
    let again: Vec<(&str, StoredRef<'_>)> =
        back.results.iter().map(|s| (s.result.step(), StoredRef { definition: s.definition, result: s.result.as_ref() })).collect();
    assert_eq!(with_results(&back, &again), bytes);
    // Without results the part is not written.
    assert_eq!(parts(&bytes_of(&f)).len(), 4);

    // Hostile results parts (their manifest entries fixed, so the part's own checks run).
    let refused = |b: Vec<u8>, what: &str| {
        let e = read(&with_part(&bytes, "results.bin", b)).unwrap_err();
        assert!(matches!(e, Error::Invalid(_) | Error::Limit(_)), "{what}: {e}");
        e.to_string()
    };
    for cut in (0..rb.len()).step_by(5) {
        refused(rb[..cut].to_vec(), &format!("cut at {cut}"));
    }
    let mut trailing = rb.clone();
    trailing.push(0);
    assert!(refused(trailing, "trailing").contains("trailing"));
    let mut many = rb.clone();
    many[8..12].copy_from_slice(&(results::MAX_STEPS as u32 + 1).to_le_bytes());
    assert!(matches!(read(&with_part(&bytes, "results.bin", many)), Err(Error::Limit(_))));
    let mut long_name = rb.clone();
    long_name[12..20].copy_from_slice(&u64::MAX.to_le_bytes());
    refused(long_name, "a name longer than the part");
    let kind_at = 12 + 8 + 4 + 32;
    let mut kind = rb.clone();
    kind[kind_at] = 9;
    assert!(refused(kind, "kind").contains("unknown kind"));
    let mut wrong_kind = rb.clone();
    wrong_kind[kind_at] = 4;
    assert!(refused(wrong_kind, "a heat result read as a fluid one").contains("damaged"));
    let mut payload_len = rb.clone();
    payload_len[kind_at + 1..kind_at + 9].copy_from_slice(&(1u64 << 40).to_le_bytes());
    refused(payload_len, "a payload longer than the part");
    let mut name = rb.clone();
    name[20] = b'1';
    refused(name, "a step name commands refuse");
    // A result that names another step, and one whose nodes do not ascend.
    let other = fieldcraft_solver::StepResult::Heat(fieldcraft_solver::HeatResult { step: "Other".into(), ..heat_result() });
    let wrong = results::encode(&[("Other", StoredRef { definition: hash, result: SolvedRef::Implicit(&other) })]).unwrap();
    let mut renamed = wrong.clone();
    renamed[20..25].copy_from_slice(b"Otxer");
    assert!(refused(renamed, "renamed").contains("step `Other`'s"));
    let unsorted = fieldcraft_solver::StepResult::Heat(fieldcraft_solver::HeatResult { nodes: vec![3, 2, 1], ..heat_result() });
    let b = results::encode(&[("Heat", StoredRef { definition: hash, result: SolvedRef::Implicit(&unsorted) })]).unwrap();
    assert!(refused(b, "unsorted").contains("do not ascend"));
    let short = fieldcraft_solver::StepResult::Heat(fieldcraft_solver::HeatResult { frames: vec![vec![1.0; 3], vec![1.0; 2]], ..heat_result() });
    let b = results::encode(&[("Heat", StoredRef { definition: hash, result: SolvedRef::Implicit(&short) })]).unwrap();
    assert!(refused(b, "a frame too short").contains("2 values for 3"));
    assert!(results::encode(&[("Heat", StoredRef { definition: hash, result: SolvedRef::Implicit(&other) })]).is_err(), "listed under another name");
    let twice = [("Heat", StoredRef { definition: hash, result: SolvedRef::Implicit(&heat) }); 2];
    assert!(results::encode(&twice).is_err());

    // Another results encoding: the model opens without them, saying why.
    let mut later = rb.clone();
    later[6..8].copy_from_slice(b"99");
    let f2 = read(&with_part(&bytes, "results.bin", later)).unwrap();
    assert!(f2.results.is_empty());
    assert!(f2.results_note.as_deref().is_some_and(|n| n.contains("MCRSLT99") && n.contains("solve again")), "{:?}", f2.results_note);
    let mut magic = rb.clone();
    magic[0] = b'X';
    refused(magic, "not a results part");

    // The part must be listed, and listed parts must be there; older formats hold none.
    let mut ps = parts(&bytes_of(&f));
    ps.push(("results.bin".into(), rb.clone()));
    assert!(read(&package(&ps)).unwrap_err().to_string().contains("unexpected part"));
    let mut ps = parts(&bytes);
    ps.pop();
    assert!(read(&package(&ps)).is_err());
    let unknown = with_manifest(&bytes, |m| m["parts"]["notes.txt"] = json!({"bytes": 1, "sha256": "0".repeat(64)}));
    assert!(read(&unknown).unwrap_err().to_string().contains("unexpected part"));
}

/// Saved views, the current one and renders round trip in the manifest's `view`, in the commands'
/// parameter shapes, outside the model part; references and names are checked.
#[test]
fn saved_views_round_trip_and_are_checked() {
    use fieldcraft_io::views::{ResultsSetup, SavedView};
    let mut f = sample();
    let camera = json!({"yaw": -70, "pitch": 18, "fit": true}).as_object().cloned();
    let contour = json!({"field": "temperature", "range": {"min": 25, "max": 45}}).as_object().cloned();
    let slice = json!({"axis": "x", "position": 14.9}).as_object().cloned();
    f.view.views = vec![SavedView {
        name: "Cutaway".into(),
        units: DisplaySystem::default(),
        camera,
        results: Some(ResultsSetup { step: "Flow".into(), contour, slice, frame: Some(10), ..ResultsSetup::default() }),
    }];
    f.view.current = Some("Cutaway".into());
    f.view.renders = vec![fieldcraft_io::views::Render {
        name: "Airflow clip".into(),
        view: "Cutaway".into(),
        caption: Some("AIO radiator".into()),
        ..Default::default()
    }];
    let bytes = bytes_of(&f);
    let back = read(&bytes).unwrap();
    assert_eq!(back.view, f.view);
    let manifest: Value = serde_json::from_slice(&parts(&bytes)[0].1).unwrap();
    assert_eq!(manifest["view"]["renders"][0]["fps"], 12);
    assert_eq!(manifest["view"]["renders"][0]["frames"], "all");
    let manifest: Value = serde_json::from_slice(&parts(&bytes)[0].1).unwrap();
    assert_eq!(manifest["view"]["current"], "Cutaway");
    assert_eq!(manifest["view"]["views"][0]["units"], "mm-t-s-N-MPa");
    assert_eq!(manifest["view"]["views"][0]["results"]["contour"]["range"]["max"], 45);
    // The model part does not change with the view: it is outside the definition.
    assert_eq!(parts(&bytes)[1].1, parts(&bytes_of(&sample()))[1].1);
    for (what, key, value) in [
        ("a current view not saved", "current", json!("Nope")),
        ("an unknown key in a view", "views", json!([{"name": "V", "units": "mm-t-s-N-MPa", "zoom": 2}])),
        ("a view with a bad step", "views", json!([{"name": "V", "units": "mm-t-s-N-MPa", "results": {"step": ""}}])),
        ("a view without units", "views", json!([{"name": "V"}])),
        ("two views of one name", "views", json!([{"name": "V", "units": "mm-t-s-N-MPa"}, {"name": "V", "units": "mm-t-s-N-MPa"}])),
        ("a render of a view not saved", "renders", json!([{"name": "R", "view": "Nope"}])),
        ("a render at 0 fps", "renders", json!([{"name": "R", "view": "Cutaway", "fps": 0}])),
        ("a render as mp4", "renders", json!([{"name": "R", "view": "Cutaway", "format": "mp4"}])),
    ] {
        let b = with_manifest(&bytes, |m| m["view"][key] = value.clone());
        assert!(read(&b).is_err(), "{what}");
    }
}

/// Geometry files embedded as parts: listed with their hashes, referenced from `geometry`, each
/// once; a copy no file names, a name the manifest does not list or bytes that are not the file
/// are refused.
#[test]
fn embedded_geometry_round_trips_and_is_checked() {
    let mut f = sample();
    let design = b"design bytes".to_vec();
    let mesh_file = b"3mf bytes".to_vec();
    // The same mesh imported twice is embedded once.
    f.geometry.meshes.push(SourceFile::of("/other/tet.3mf", &mesh_file).unwrap());
    let embed = [(fieldcraft_io::sha256_hex(&design), design.as_slice()), (fieldcraft_io::sha256_hex(&mesh_file), mesh_file.as_slice())]
        .into_iter()
        .collect();
    let extras = fieldcraft_io::Extras { embed, ..Default::default() };
    let bytes = fieldcraft_io::write_with(&f.model, &f.geometry, &f.journal, &f.journal_header, &f.view, &extras).unwrap();
    let names: Vec<String> = parts(&bytes).into_iter().map(|(n, _)| n).collect();
    assert_eq!(names, ["fieldcraft.json", "model.json", "mesh.bin", "journal.json", "geometry/plate.solidcraft", "geometry/tet.3mf"]);
    let back = read(&bytes).unwrap();
    assert_eq!(back.geometry.design.as_ref().unwrap().embedded(), Some("geometry/plate.solidcraft"));
    assert_eq!(back.geometry.design.as_ref().unwrap().path(), "parts/plate.solidcraft", "the reference stays");
    assert!(back.geometry.meshes.iter().all(|m| m.embedded() == Some("geometry/tet.3mf")));
    assert_eq!(back.embedded.get(&fieldcraft_io::sha256_hex(&design)), Some(&design));
    assert_eq!(back.embedded.len(), 2);
    // Without embedding, the references only.
    assert!(read(&bytes_of(&back)).unwrap().geometry.design.unwrap().embedded().is_none());

    // Two files of one name get two parts.
    let mut g = sample();
    g.geometry.meshes = vec![SourceFile::of("a/plate.solidcraft", b"other").unwrap()];
    let embed = [(fieldcraft_io::sha256_hex(&design), design.as_slice()), (fieldcraft_io::sha256_hex(b"other"), &b"other"[..])].into_iter().collect();
    let b = fieldcraft_io::write_with(
        &g.model,
        &g.geometry,
        &g.journal,
        &g.journal_header,
        &g.view,
        &fieldcraft_io::Extras { embed, ..Default::default() },
    )
    .unwrap();
    let names: Vec<String> = parts(&b).into_iter().map(|(n, _)| n).collect();
    assert!(names.contains(&"geometry/plate.solidcraft".to_string()) && names.contains(&"geometry/plate-2.solidcraft".to_string()), "{names:?}");
    // Bytes that are not the file's content are not embedded.
    let wrong = [(fieldcraft_io::sha256_hex(&design), &b"tampered"[..])].into_iter().collect();
    assert!(
        fieldcraft_io::write_with(
            &g.model,
            &g.geometry,
            &g.journal,
            &g.journal_header,
            &g.view,
            &fieldcraft_io::Extras { embed: wrong, ..Default::default() }
        )
        .is_err()
    );

    let refused = |b: Vec<u8>, what: &str| assert!(matches!(read(&b), Err(Error::Invalid(_))), "{what}");
    // A copy no file names.
    refused(with_manifest(&bytes, |m| m["geometry"]["design"].as_object_mut().unwrap().remove("embedded").map(|_| ()).unwrap()), "orphan");
    // A name that is not listed, or not a geometry part.
    refused(with_manifest(&bytes, |m| m["geometry"]["design"]["embedded"] = json!("geometry/other.solidcraft")), "unlisted");
    for bad in ["geometry/../x", "geometry/", "geometry/.hidden", "model.json", "geometry/a/b", "geometry/sp ace"] {
        refused(with_manifest(&bytes, |m| m["geometry"]["design"]["embedded"] = json!(bad)), bad);
    }
    // A part whose bytes are not the file named.
    refused(with_part(&bytes, "geometry/plate.solidcraft", b"not the design".to_vec()), "other bytes");
    // An embedded part in a file of an older format.
    let old = with_manifest(&bytes, |m| m["format"] = json!("fieldcraft/2"));
    assert!(read(&old).is_err());
}

/// Preferred run settings round trip in the manifest, outside the model part, and are checked.
#[test]
fn run_settings_round_trip_and_are_checked() {
    use fieldcraft_io::{Backend, RunSettings};
    let f = sample();
    let run: fieldcraft_io::Run = [("Flow".to_string(), RunSettings { backend: Backend::Gpu, fallback: Some(Backend::Cpu) })].into_iter().collect();
    let extras = fieldcraft_io::Extras { run: run.clone(), ..Default::default() };
    let bytes = fieldcraft_io::write_with(&f.model, &f.geometry, &f.journal, &f.journal_header, &f.view, &extras).unwrap();
    assert_eq!(read(&bytes).unwrap().run, run);
    let manifest: Value = serde_json::from_slice(&parts(&bytes)[0].1).unwrap();
    assert_eq!(manifest["run"], json!({"Flow": {"backend": "gpu", "fallback": "cpu"}}));
    assert_eq!(parts(&bytes)[1].1, parts(&bytes_of(&f))[1].1, "outside the definition");
    for (what, run) in [
        ("falling back on itself", json!({"Flow": {"backend": "gpu", "fallback": "gpu"}})),
        ("an unknown backend", json!({"Flow": {"backend": "tpu"}})),
        ("a bad step name", json!({"1 x": {"backend": "cpu"}})),
        ("an unknown key", json!({"Flow": {"backend": "cpu", "threads": 4}})),
    ] {
        assert!(read(&with_manifest(&bytes, |m| m["run"] = run.clone())).is_err(), "{what}");
    }
    let bad = fieldcraft_io::Extras {
        run: [("Flow".to_string(), RunSettings { backend: Backend::Cpu, fallback: Some(Backend::Cpu) })].into_iter().collect(),
        ..Default::default()
    };
    assert!(fieldcraft_io::write_with(&f.model, &f.geometry, &f.journal, &f.journal_header, &f.view, &bad).is_err());
}
