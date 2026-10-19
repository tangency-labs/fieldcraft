//! A solid with a void shell (a box shelled, then capped) opens with every shell and everything
//! joined to them (#364); a feature that fails to rebuild is reported naming the body it changes.

use fieldcraft_geometry::{Design, open_design};

const PI: f64 = std::f64::consts::PI;
/// The outer box, 40 × 30 × 20 mm, and its cavity with 2 mm walls and floor.
const BOX: f64 = 40.0 * 30.0 * 20.0;
const CAVITY: f64 = 36.0 * 26.0 * 18.0;
const LID: f64 = 40.0 * 30.0 * 2.0;
/// A post ⌀6 × 5 mm standing on the floor.
const POST: f64 = PI * 9.0 * 5.0;

fn open(bytes: &[u8]) -> Design {
    let design = open_design(bytes).unwrap();
    assert!(design.errors.is_empty(), "{:?}", design.errors);
    assert_eq!(design.bodies.len(), 1);
    design
}

fn check(design: &Design, volume: f64, faces: usize) {
    let body = design.body("Case").unwrap();
    let m = body.measure().unwrap();
    assert!((m.volume - volume).abs() < 1e-5 * volume, "{} vs {volume}", m.volume);
    assert_eq!(m.faces, faces);
    let mut names: Vec<&str> = body.faces().iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names.len(), faces);
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), faces, "every face has its own name");
    // The mesh sees the same solid: every face tessellated, the same enclosed volume.
    let mesh = body.tessellate(0.01).unwrap();
    let mut tags = mesh.tri_face.clone();
    tags.sort_unstable();
    tags.dedup();
    assert_eq!(tags.len(), faces);
    assert!((mesh.enclosed_volume() - volume).abs() < 1e-3 * volume, "{} vs {volume}", mesh.enclosed_volume());
}

#[test]
fn a_closed_hollow_box_keeps_its_void_shell() {
    // Six outside faces and six cavity faces.
    check(&open(include_bytes!("designs/hollow-box.solidcraft")), BOX - CAVITY + LID, 12);
}

#[test]
fn a_post_joined_to_the_void_floor_comes_through() {
    // The hollow box's twelve faces, the post's top and its side (two halves).
    check(&open(include_bytes!("designs/hollow-box-post.solidcraft")), BOX - CAVITY + LID + POST, 15);
}

#[test]
fn a_post_in_an_open_tray_comes_through() {
    // Five outside faces, the rim, five cavity faces and the post's three.
    check(&open(include_bytes!("designs/open-tray-post.solidcraft")), BOX - CAVITY + POST, 14);
}

#[test]
fn a_feature_that_fails_names_the_body_it_changes() {
    let text = std::str::from_utf8(include_bytes!("designs/hollow-box-post.solidcraft")).unwrap();
    assert!(text.contains(r#""height": "5""#));
    let broken = text.replace(r#""height": "5""#, r#""height": "0""#);
    let design = open_design(broken.as_bytes()).unwrap();
    let [error] = design.errors.as_slice() else { panic!("{:?}", design.errors) };
    assert_eq!(error.feature, "Cylinder1");
    assert_eq!(error.bodies, ["Case"]);
    let warning = error.warning();
    assert!(warning.starts_with(r#"Cylinder1 failed to build, so body "Case" may differ from the design (SolidCraft: "#), "{warning}");
}
