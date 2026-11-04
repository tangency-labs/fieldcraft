//! Every design SolidCraft saved for these tests opens (`tests/designs`): the format 21 designs
//! and those of later SolidCraft releases, which SolidCraft's reader migrates. A format 21 design
//! re-saved by a later release (`resaved/`) opens with every face name on the same face, as a
//! person's design does after they save it again; a recipe built afresh by a later release opens
//! as the same bodies with the same face names (#605: FieldCraft 0.7.0 refused SolidCraft 0.16.2's
//! designs; #655: the 0.16.2 pin refused 0.17's and 0.18's).

use std::path::{Path, PathBuf};

use fieldcraft_geometry::{DESIGN_FORMAT, Design, open_design};

fn designs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/designs")
}

/// Every `.solidcraft` file under `dir`, sorted.
fn saved(dir: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            files.extend(saved(&path));
        } else if path.extension().is_some_and(|e| e == "solidcraft") {
            files.push(path);
        }
    }
    files.sort();
    files
}

fn open(path: &Path) -> Design {
    let design = open_design(&std::fs::read(path).unwrap()).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert!(design.errors.is_empty(), "{}: {:?}", path.display(), design.errors);
    design
}

/// Bodies by name with their volume and face names, in the design's order.
fn summary(design: &Design) -> Vec<(String, f64, Vec<String>)> {
    design.bodies.iter().map(|b| (b.name.clone(), b.measure().unwrap().volume, b.faces().iter().map(|f| f.name.clone()).collect())).collect()
}

#[test]
fn every_saved_design_opens() {
    let files = saved(&designs());
    assert!(files.len() >= 4, "{files:?}");
    for path in &files {
        let design = open(path);
        assert!(!design.bodies.is_empty(), "{}", path.display());
        for body in &design.bodies {
            let mut names: Vec<&str> = body.faces().iter().map(|f| f.name.as_str()).collect();
            let n = names.len();
            names.sort_unstable();
            names.dedup();
            assert!(n > 0 && names.len() == n, "{}: every face of {} has its own name", path.display(), body.name);
        }
    }
}

/// Face names sorted, for a recipe built afresh: SolidCraft 0.30 numbers the pieces of a split
/// face in its kernel's order for features made in it (a migrated feature keeps the old order,
/// `kernel_piece_order`), so a fresh build may give the same names to the pieces in another order.
fn sorted(mut names: Vec<String>) -> Vec<String> {
    names.sort_unstable();
    names
}

#[test]
fn a_later_release_saves_the_same_design() {
    let (mut fresh, mut resaved) = (0, 0);
    for path in saved(&designs()).iter().filter(|p| p.parent() != Some(designs().as_path())) {
        let reference = designs().join(path.file_name().unwrap());
        if !reference.exists() {
            continue;
        }
        let again = path.parent().and_then(Path::file_name).is_some_and(|d| d == "resaved");
        let (a, b) = (summary(&open(&reference)), summary(&open(path)));
        assert_eq!(a.len(), b.len(), "{}", path.display());
        for ((name_a, volume_a, faces_a), (name_b, volume_b, faces_b)) in a.into_iter().zip(b) {
            assert_eq!(name_a, name_b, "{}", path.display());
            assert!((volume_a - volume_b).abs() <= 1e-9 * volume_a, "{}: {name_a} {volume_a} vs {volume_b}", path.display());
            if again {
                // A re-saved design keeps each name on its face.
                assert_eq!(faces_a, faces_b, "{}: {name_a}", path.display());
            } else {
                assert_eq!(sorted(faces_a), sorted(faces_b), "{}: {name_a}", path.display());
            }
        }
        if again {
            resaved += 1;
        } else {
            fresh += 1;
        }
    }
    assert!(fresh >= 1 && resaved >= 1, "saved copies: {fresh} fresh, {resaved} re-saved");
}

/// SolidCraft 0.30 saves design format 48; this FieldCraft reads it.
const _: () = assert!(DESIGN_FORMAT >= 48);

/// The hollow box with a post as each release saved it: one body, `Case`, of 15 faces.
#[test]
fn reads_every_release_s_designs() {
    for release in ["0.16.2", "0.17.0", "0.18.1", "0.18.2", "0.30.1"] {
        let path = designs().join(format!("solidcraft-{release}/hollow-box-post.solidcraft"));
        let design = open(&path);
        let [body] = design.bodies.as_slice() else { panic!("{release}: {:?}", summary(&design)) };
        assert_eq!(body.name, "Case", "{release}");
        assert_eq!(body.faces().len(), 15, "{release}");
    }
}
