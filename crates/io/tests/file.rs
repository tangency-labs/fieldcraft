//! `.fieldcraft` files: round trips, stable bytes and hostile files.

use std::io::{Cursor, Read, Write};

use fieldcraft_io::{Error, Geometry, Item, ItemKind, ModelFile, SourceFile, View, mesh, read, write};
use fieldcraft_model::{ElementKind, Entity, GeometryRef, Model, UnitSystem};
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
    m.set_units(UnitSystem::from_label("mm-t-s-N-MPa").unwrap());
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
    let view = View { hidden: [Item { kind: ItemKind::Datum, name: "Sym".into() }, Item { kind: ItemKind::Body, name: "Body1".into() }].into() };
    ModelFile { model: m, geometry, journal, view }
}

fn bytes_of(f: &ModelFile) -> Vec<u8> {
    write(&f.model, &f.geometry, &f.journal, &f.view).unwrap()
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
        json!({"units": "furlongs"}),
    ] {
        let b = with_part(&good, "model.json", serde_json::to_vec(&model).unwrap());
        assert!(read(&b).is_err(), "{model}");
    }

    // A key given twice is an error in a current file, not the last one silently winning.
    let twice = with_part(&good, "model.json", br#"{"units": "mm-t-s-N-MPa", "units": "m-kg-s-N-Pa"}"#.to_vec());
    let e = read(&twice).unwrap_err();
    assert!(e.to_string().contains("duplicate field"), "{e}");
    let mut ps = parts(&good);
    let manifest = String::from_utf8(ps[0].1.clone()).unwrap();
    ps[0].1 = manifest.replacen("\"generator\"", "\"generator\": \"FieldCraft 0.0.1\",\n  \"generator\"", 1).into_bytes();
    assert!(read(&package(&ps)).unwrap_err().to_string().contains("duplicate field"));

    // Journals: not a list, entries without a command, unknown fields.
    for journal in [json!({}), json!([{"params": {}}]), json!([{"command": ""}]), json!([{"command": "a.b", "when": 1}]), json!([7])] {
        let b = with_part(&good, "journal.json", serde_json::to_vec(&journal).unwrap());
        assert!(read(&b).is_err(), "{journal}");
    }
    let ok = with_part(&good, "journal.json", br#"[{"command": "model.check"}]"#.to_vec());
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
    let plain = write(&f.model, &f.geometry, &f.journal, &empty).unwrap();
    let manifest: Value = serde_json::from_slice(&parts(&plain)[0].1).unwrap();
    assert!(manifest.get("view").is_none(), "an empty view writes nothing: {manifest}");
}
