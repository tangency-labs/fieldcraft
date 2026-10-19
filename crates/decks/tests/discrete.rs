//! Trusses (rod3), point masses, rotary inertia, springs and dashpots in every deck format:
//! what each writes, and what each reports as unsupported.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::Error;
use fieldcraft_model::{
    BoundaryCondition, Direction, Elastic, ElementKind, Formulation, Mass, Material, Model, PointMass, Procedure, ReferencePoint, Section,
    SectionKind, Set, Spring, Step, Target,
};

/// A rod2 (or rod3) bar from node 1 to node 2 and a rod2 on to node 4, a reference point beyond
/// it; a motor mass at the tip, a hub with inertia on the point, a spring to the ground across
/// the tip, a spring along the line from the tip to the hub, a dashpot to the ground.
fn model(rod3: bool) -> Model {
    let mut m = Model::default();
    for (l, x) in [(1, 0.0), (2, 1000.0), (3, 500.0), (4, 2000.0)] {
        m.set_node(l, [x, 0.0, 0.0]).unwrap();
    }
    if rod3 {
        m.set_element(1, ElementKind::Rod3, vec![1, 2, 3]).unwrap();
    } else {
        m.set_element(1, ElementKind::Rod2, vec![1, 2]).unwrap();
    }
    m.set_element(2, ElementKind::Rod2, vec![2, 4]).unwrap();
    m.add_set("bars", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    m.add_set("base", Set { nodes: vec![1], ..Set::default() }).unwrap();
    m.add_set("tip", Set { nodes: vec![4], ..Set::default() }).unwrap();
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 210_000.0, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() },
    )
    .unwrap();
    m.add_section(
        "bars",
        Section { set: "bars".into(), material: "Steel".into(), kind: SectionKind::Rod { area: 20.0 }, formulation: Formulation::default() },
    )
    .unwrap();
    m.add_reference_point("Hub", ReferencePoint { xyz: [2500.0, 0.0, 0.0], label: 0 }).unwrap();
    let tip = Target::Set("tip".into());
    let hub = Target::ReferencePoint("Hub".into());
    m.add_mass("motor", PointMass { target: tip.clone(), mass: Some(Mass::Isotropic(0.002)), inertia: None }).unwrap();
    m.add_mass("hub", PointMass { target: hub.clone(), mass: Some(Mass::Isotropic(0.01)), inertia: Some([1.0, 2.0, 3.0, 0.1, 0.0, 0.0]) }).unwrap();
    m.add_spring("mount", Spring::Spring { first: tip.clone(), second: None, direction: Direction::Dof(2), stiffness: 50.0 }).unwrap();
    m.add_spring("link", Spring::Spring { first: tip.clone(), second: Some(hub), direction: Direction::Line, stiffness: 1000.0 }).unwrap();
    m.add_spring("damper", Spring::Dashpot { first: tip, second: None, direction: Direction::Vector([0.0, 0.0, -2.0]), coefficient: 0.5 }).unwrap();
    m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    let mut step = Step::new("Modes", Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
    step.boundary_conditions.insert("fix".into());
    m.add_step(step, None).unwrap();
    m
}

fn unsupported(e: std::result::Result<impl std::fmt::Debug, Error>) -> String {
    match e {
        Err(Error::Unsupported { what, .. }) => what,
        other => panic!("expected unsupported, got {other:?}"),
    }
}

#[test]
fn abaqus_writes_every_discrete_element_and_rod3() {
    let deck = abaqus::write(&model(true)).unwrap();
    for expected in [
        "*ELEMENT, TYPE=T3D3\n1, 1, 2, 3\n",
        "*ELEMENT, TYPE=MASS, ELSET=hub-mass\n3, 99999999\n",
        "*ELEMENT, TYPE=MASS, ELSET=motor-mass\n4, 4\n",
        "*ELEMENT, TYPE=ROTARYI, ELSET=hub-inertia\n5, 99999999\n",
        "*ELEMENT, TYPE=DASHPOT1, ELSET=damper-dashpot\n6, 4\n",
        "*ELEMENT, TYPE=SPRINGA, ELSET=link-spring\n7, 4, 99999999\n",
        "*ELEMENT, TYPE=SPRING1, ELSET=mount-spring\n8, 4\n",
        "*MASS, ELSET=hub-mass\n0.01\n",
        "*ROTARY INERTIA, ELSET=hub-inertia\n1.0, 2.0, 3.0, 0.1, 0.0, 0.0\n",
        "*DASHPOT, ELSET=damper-dashpot\n3\n0.5\n",
        "*SPRING, ELSET=link-spring\n,\n1000.0\n",
        "*SPRING, ELSET=mount-spring\n2\n50.0\n",
    ] {
        assert!(deck.contains(expected), "{expected}\n{deck}");
    }
    // An element set name that a model set already has is made unique.
    let mut m = model(false);
    m.add_set("motor-mass", Set { nodes: vec![4], ..Set::default() }).unwrap();
    assert!(abaqus::write(&m).unwrap().contains("ELSET=motor-mass-2\n"));
    // Skew directions and masses per direction need orientations.
    let mut m = model(false);
    m.replace_spring(
        "mount",
        Spring::Spring { first: Target::Set("tip".into()), second: None, direction: Direction::Vector([0.0, 1.0, 1.0]), stiffness: 1.0 },
    )
    .unwrap();
    m.replace_mass("motor", PointMass { target: Target::Set("tip".into()), mass: Some(Mass::Directional([1.0, 2.0, 3.0])), inertia: None }).unwrap();
    let what = unsupported(abaqus::write(&m));
    assert!(what.contains("spring `mount`") && what.contains("point mass `motor`"), "{what}");
}

#[test]
fn nastran_writes_conm1_celas2_cdamp2() {
    let deck = nastran::write(&model(false)).unwrap();
    eprintln!("{deck}");
    for expected in ["CONM1*", "CELAS2*", "CDAMP2*", "$ point mass hub", "$ spring link", "$ dashpot damper"] {
        assert!(deck.contains(expected), "{expected}\n{deck}");
    }
    // The hub's 6 × 6 mass matrix: mass on the translations, the inertia tensor on the rotations.
    let conm1: Vec<&str> = deck.lines().skip_while(|l| !l.starts_with("CONM1*")).take(6).collect();
    let fields: Vec<f64> = conm1
        .iter()
        .flat_map(|l| l.get(8..).unwrap_or("").as_bytes().chunks(16).map(|c| std::str::from_utf8(c).unwrap().trim().to_string()))
        .filter(|t| !t.is_empty())
        .filter_map(|t| t.parse().ok())
        .collect();
    assert_eq!(fields.first().copied(), Some(3.0), "element id after the mesh's: {conm1:?}");
    assert!(
        fields.ends_with(&[0.01, 0.0, 0.01, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.1, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0]),
        "{fields:?}"
    );
    // The line spring runs along X: component 1 at both grids; the dashpot along −Z: component 3.
    let link = deck.lines().skip_while(|l| !l.contains("spring link")).nth(1).unwrap().to_string();
    assert!(link.starts_with("CELAS2*") && link.contains("1000.") && link.ends_with("4               1"), "{link}");
    let rest = deck.lines().skip_while(|l| !l.contains("spring link")).nth(2).unwrap().to_string();
    assert!(rest.contains("99999999") && rest.trim_end().ends_with('1'), "{rest}");
    let what = unsupported(nastran::write(&model(true)));
    assert!(what.contains("rod3"), "{what}");
    let mut m = model(false);
    m.replace_spring(
        "link",
        Spring::Spring { first: Target::Set("tip".into()), second: None, direction: Direction::Vector([1.0, 1.0, 0.0]), stiffness: 1.0 },
    )
    .unwrap();
    assert!(unsupported(nastran::write(&m)).contains("spring `link`"));
}

#[test]
fn lsdyna_writes_masses_inertia_and_discrete_springs() {
    let deck = lsdyna::write(&model(false)).unwrap().text;
    eprintln!("{deck}");
    for expected in [
        "*ELEMENT_MASS",
        "*ELEMENT_INERTIA",
        "*MAT_SPRING_ELASTIC_TITLE\nlink\n",
        "*MAT_DAMPER_VISCOUS_TITLE\ndamper\n",
        "*SECTION_DISCRETE",
        "*DEFINE_SD_ORIENTATION",
        "*ELEMENT_DISCRETE",
        "*BOUNDARY_SPC_NODE",
    ] {
        assert!(deck.contains(expected), "{expected}\n{deck}");
    }
    // Two grounded links (damper, mount) get two ground nodes after the mesh's (5, 6); the
    // line spring has no orientation vector (vid 0, along its nodes).
    let discrete: Vec<Vec<i64>> = deck
        .lines()
        .skip_while(|l| *l != "*ELEMENT_DISCRETE")
        .skip(1)
        .take_while(|l| !l.starts_with('*'))
        .filter(|l| !l.starts_with('$'))
        .map(|l| l.split_whitespace().map(|f| f.parse().unwrap()).collect())
        .collect();
    assert_eq!(discrete.len(), 3);
    assert_eq!(discrete[0][2..4], [4, 5], "damper to its ground node");
    assert_eq!(discrete[1][2..4], [4, 99_999_999], "link to the hub");
    assert_eq!(discrete[1].get(4).copied().unwrap_or(0), 0, "along the line");
    assert_eq!(discrete[2][2..4], [4, 6], "mount to its ground node");
    assert!(unsupported(lsdyna::write(&model(true))).contains("rod3"));
    let mut m = model(false);
    m.replace_mass("motor", PointMass { target: Target::Set("tip".into()), mass: Some(Mass::Directional([1.0, 2.0, 3.0])), inertia: None }).unwrap();
    assert!(unsupported(lsdyna::write(&m)).contains("point mass `motor`"));
}
