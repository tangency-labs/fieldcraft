//! Contact definitions in every deck format: Abaqus surface interactions, contact pairs and
//! general contact; LS-DYNA automatic contacts; Nastran refuses them by name.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::Error;
use fieldcraft_model::{
    BoundaryCondition, Contact, ContactSurfaces, Elastic, ElementFace, ElementKind, Friction, InitialPenetration, Material, Model, Procedure,
    Section, SectionKind, Set, Step, Target, Tie,
};

/// Two unit hex8 cubes, the second on top of the first (separate nodes), with face sets for the
/// touching faces and a static step.
fn model() -> Model {
    let mut m = Model::default();
    let corners = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    for cube in 0..2u32 {
        for (i, c) in corners.iter().enumerate() {
            m.set_node(8 * cube + i as u32 + 1, [c[0], c[1], c[2] + cube as f64]).unwrap();
        }
        m.set_element(cube + 1, ElementKind::Hex8, (8 * cube + 1..=8 * cube + 8).collect()).unwrap();
    }
    m.add_set("all", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    m.add_set("lower", Set { faces: vec![ElementFace { element: 1, face: 2 }], ..Set::default() }).unwrap();
    m.add_set("upper", Set { faces: vec![ElementFace { element: 2, face: 1 }], ..Set::default() }).unwrap();
    m.add_set("base", Set { nodes: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 210_000.0, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() },
    )
    .unwrap();
    m.add_section("all", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Default::default() })
        .unwrap();
    m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    let mut step = Step::new("Load", Procedure::Static {});
    step.boundary_conditions.insert("fix".into());
    m.add_step(step, None).unwrap();
    let pair = Contact {
        friction: Friction { static_coefficient: 0.3, kinetic: Some(0.2), decay: 5.0 },
        initial_penetration: InitialPenetration::Adjust,
        ..Contact::new(ContactSurfaces::Pair { main: "lower".into(), secondary: "upper".into() })
    };
    m.add_contact("stack", pair).unwrap();
    m.add_contact("everything", Contact { friction: Friction::coulomb(0.1), ..Contact::new(ContactSurfaces::General { surface: None }) }).unwrap();
    m.add_contact("fold", Contact { penalty_scale: 0.5, ..Contact::new(ContactSurfaces::General { surface: Some("upper".into()) }) }).unwrap();
    m.add_tie("glue", Tie { position_tolerance: Some(0.05), ..Tie::new("lower", "upper") }).unwrap();
    m
}

#[test]
fn abaqus_writes_interactions_pairs_and_general_contact() {
    let deck = abaqus::write(&model()).unwrap();
    for expected in [
        "*SURFACE, TYPE=ELEMENT, NAME=lower\n1, S2\n",
        "*SURFACE INTERACTION, NAME=stack-interaction\n*FRICTION, EXPONENTIAL DECAY\n0.3, 0.2, 5.0\n",
        "*CONTACT PAIR, INTERACTION=stack-interaction, TYPE=SURFACE TO SURFACE, ADJUST=0.0\nupper, lower\n",
        "*SURFACE INTERACTION, NAME=everything-interaction\n*FRICTION\n0.1\n",
        "*SURFACE INTERACTION, NAME=fold-interaction\n*FRICTION\n0.0\n",
        "*CONTACT\n*CONTACT INCLUSIONS, ALL EXTERIOR\n*CONTACT INCLUSIONS\nupper,\n*CONTACT PROPERTY ASSIGNMENT\n, , everything-interaction\nupper, upper, fold-interaction\n",
        "** contact `fold`: penalty scale 0.5",
        "*TIE, NAME=glue, ADJUST=NO, POSITION TOLERANCE=0.05\nupper, lower\n",
    ] {
        assert!(deck.contains(expected), "{expected}\n{deck}");
    }
    // Model data: before the first step.
    assert!(deck.find("*CONTACT PAIR") < deck.find("*STEP"));
}

#[test]
fn lsdyna_writes_automatic_contacts_on_segment_sets() {
    let deck = lsdyna::write(&model()).unwrap().text;
    let after = |keyword: &str| -> Vec<String> {
        let at = deck.find(keyword).unwrap_or_else(|| panic!("{keyword}\n{deck}"));
        deck[at..].lines().skip(1).filter(|l| !l.starts_with('$')).take(3).map(|l| l.split_whitespace().collect::<Vec<_>>().join(" ")).collect()
    };
    // Sets are numbered in name order: all 1, base 2, lower 3, upper 4.
    assert_eq!(after("*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE"), ["4 3 0 0", "0.3 0.2 5.0", "1.0 1.0"], "{deck}");
    let single: Vec<&str> = deck.match_indices("*CONTACT_AUTOMATIC_SINGLE_SURFACE").map(|(i, _)| &deck[i..]).collect();
    assert_eq!(single.len(), 2, "{deck}");
    let cards = |text: &str| -> Vec<String> {
        text.lines().skip(1).filter(|l| !l.starts_with('$')).take(3).map(|l| l.split_whitespace().collect::<Vec<_>>().join(" ")).collect()
    };
    // `everything` (every part), then `fold` (the upper face's segment set).
    assert_eq!(cards(single[0]), ["0 5", "0.1 0.1 0.0", "1.0 1.0"], "{deck}");
    assert_eq!(cards(single[1]), ["4 0", "0.0 0.0 0.0", "0.5 0.5"], "{deck}");
    for set in ["lower", "upper"] {
        assert!(deck.contains(&format!("*SET_SEGMENT_TITLE\n{set}\n")), "{set}\n{deck}");
    }
    assert_eq!(after("*CONTACT_TIED_SURFACE_TO_SURFACE"), ["4 3 0 0", "0.0 0.0 0.0", "1.0 1.0"], "{deck}");
}

#[test]
fn nastran_refuses_contact_by_name() {
    match nastran::write(&model()) {
        Err(Error::Unsupported { what, .. }) => {
            for name in ["everything", "fold", "stack"] {
                assert!(what.contains(&format!("contact `{name}`")), "{what}");
            }
            assert!(what.contains("tie `glue`"), "{what}");
        }
        other => panic!("{other:?}"),
    }
}
