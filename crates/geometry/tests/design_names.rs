//! The persistent names SolidCraft gives a split body's faces, edges and vertices, frozen
//! (`tests/names/*.txt`). Models name geometry by these names, so a SolidCraft that names the
//! same entity otherwise leaves a saved model's sets, loads and seeds on names the design no
//! longer has (#605: SolidCraft 0.16 named `F9:face:4` `F9:section`). A pin bump that renames
//! anything fails here and shows the renames; bless them (`FIELDCRAFT_BLESS=1`) only once
//! opening older models finds them again (`fieldcraft_engine`'s `names` module) or they are
//! listed for the release notes.
#![cfg(feature = "scripts")]

use std::path::PathBuf;

use fieldcraft_geometry::{open_design, script_design};

/// Every body's face, edge and vertex names, a line each, in the design's order.
fn names(recipe: &str) -> String {
    let design = open_design(&script_design(recipe).unwrap()).unwrap();
    assert!(design.errors.is_empty(), "{:?}", design.errors);
    let mut out = String::new();
    for b in &design.bodies {
        out += &format!("body {}\n", b.name);
        for f in b.faces() {
            out += &format!("  face {}\n", f.name);
        }
        for e in b.edges().unwrap() {
            out += &format!("  edge {}\n", e.name);
        }
        for v in b.vertices().unwrap() {
            out += &format!("  vertex {}\n", v.name);
        }
    }
    out
}

fn golden(name: &str, now: &str) {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/names").join(format!("{name}.txt"));
    if std::env::var_os("FIELDCRAFT_BLESS").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, now).unwrap();
    }
    let frozen = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e} (FIELDCRAFT_BLESS=1 writes it)", path.display()));
    if frozen != now {
        let (a, b): (Vec<&str>, Vec<&str>) = (frozen.lines().collect(), now.lines().collect());
        let gone: Vec<&&str> = a.iter().filter(|l| !b.contains(l)).collect();
        let new: Vec<&&str> = b.iter().filter(|l| !a.contains(l)).collect();
        panic!("{name}: SolidCraft names the geometry otherwise ({} names gone, {} new)\ngone: {gone:#?}\nnew: {new:#?}", gone.len(), new.len());
    }
}

#[test]
fn the_split_bracket_keeps_its_names() {
    golden("filleted-bracket-blocks", &names(include_str!("../../../tests/corpus/filleted-bracket-blocks/recipe.json")));
}
