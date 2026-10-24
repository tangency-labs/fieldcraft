//! 2D models (plane stress, plane strain, axisymmetric) in every deck format: element names,
//! section data, edge surfaces and pressures, centrifugal loads, and what a format cannot write.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::{Error, edge_pressure_loads};
use fieldcraft_model::{
    BoundaryCondition, Elastic, ElementEdge, ElementKind, Formulation, Integration, Load, Material, Model, Procedure, Section, SectionKind, Set,
    Step, Target,
};

/// Two quad4 over [1, 3] × [0, 1] (or the given kind on the same corners), the left edge held,
/// a pressure on the right edge and a centrifugal load.
fn plate(kind: SectionKind, element: ElementKind, formulation: Formulation) -> Model {
    let mut m = Model::default();
    for (i, x) in [1.0, 2.0, 3.0].iter().enumerate() {
        m.set_node(i as u32 + 1, [*x, 0.0, 0.0]).unwrap();
        m.set_node(i as u32 + 4, [*x, 1.0, 0.0]).unwrap();
    }
    match element {
        ElementKind::Tri3 => {
            for (l, n) in [(1, [1, 2, 5]), (2, [1, 5, 4]), (3, [2, 3, 6]), (4, [2, 6, 5])] {
                m.set_element(l, ElementKind::Tri3, n.to_vec()).unwrap();
            }
        }
        _ => {
            m.set_element(1, ElementKind::Quad4, vec![1, 2, 5, 4]).unwrap();
            m.set_element(2, ElementKind::Quad4, vec![2, 3, 6, 5]).unwrap();
        }
    }
    let right = if element == ElementKind::Tri3 { ElementEdge { element: 3, edge: 2 } } else { ElementEdge { element: 2, edge: 2 } };
    let all: Vec<u32> = m.elements().keys().copied().collect();
    m.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    m.add_set("left", Set { nodes: vec![1, 4], ..Set::default() }).unwrap();
    m.add_set("right", Set { edges: vec![right], ..Set::default() }).unwrap();
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 210_000.0, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() },
    )
    .unwrap();
    m.add_section("plate", Section { set: "all".into(), material: "Steel".into(), kind, formulation }).unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("left".into()) }).unwrap();
    m.add_load("press", Load::Pressure { set: "right".into(), magnitude: 5.0 }).unwrap();
    let axis = if kind == (SectionKind::Axisymmetric {}) { [0.0, 1.0, 0.0] } else { [0.0, 0.0, 1.0] };
    m.add_load("spin", Load::Centrifugal { angular_velocity: 100.0, point: [0.0; 3], axis, set: None }).unwrap();
    let mut s = Step::new("Load", Procedure::Static {});
    s.boundary_conditions.insert("hold".into());
    s.loads.extend(["press".to_string(), "spin".to_string()]);
    m.add_step(s, None).unwrap();
    m
}

fn unsupported(r: Result<String, Error>) -> String {
    match r {
        Err(Error::Unsupported { what, .. }) => what,
        other => panic!("expected unsupported, got {other:?}"),
    }
}

#[test]
fn production_reader_retains_planar_sections_and_edge_surfaces() {
    use fieldcraft_model::units::Preset;
    for kind in [SectionKind::PlaneStress { thickness: 2.0 }, SectionKind::PlaneStrain { thickness: 2.0 }, SectionKind::Axisymmetric {}] {
        for integration in [Integration::Full, Integration::Reduced] {
            let m = plate(kind, ElementKind::Quad4, Formulation { integration, ..Formulation::default() });
            let written = abaqus::write(&m).unwrap();
            let read = fieldcraft_decks::abaqus::read_in(&written, Preset::MmTonneS).unwrap();
            assert_eq!(read.model.nodes(), common::si(&m).nodes());
            assert_eq!(read.model.set("right").unwrap().edges, m.set("right").unwrap().edges);
            assert!(read.model.sections().values().any(|s| common::si(&m).sections().values().any(|other| other == s)));
        }
    }
}

#[test]
fn abaqus_writes_plane_and_axisymmetric_elements() {
    let full = Formulation::default();
    let im = Formulation { incompatible_modes: true, ..Formulation::default() };
    let red = Formulation { integration: Integration::Reduced, ..Formulation::default() };
    for (kind, f, want) in [
        (SectionKind::PlaneStress { thickness: 2.0 }, full, "CPS4"),
        (SectionKind::PlaneStress { thickness: 2.0 }, im, "CPS4I"),
        (SectionKind::PlaneStrain { thickness: 1.0 }, red, "CPE4R"),
        (SectionKind::Axisymmetric {}, im, "CAX4I"),
    ] {
        let deck = abaqus::write(&plate(kind, ElementKind::Quad4, f)).unwrap();
        assert!(deck.contains(&format!("*ELEMENT, TYPE={want}\n")), "{want}\n{deck}");
        // The right edge is surface S2 of element 2; the pressure is a surface load on it.
        assert!(deck.contains("*SURFACE, TYPE=ELEMENT, NAME=right\n2, S2\n"), "{deck}");
        assert!(deck.contains("right, P, 5.0\n"), "{deck}");
        assert!(deck.contains("all, CENTRIF, 10000.0, 0.0, 0.0, 0.0, "), "{deck}");
        let section = match kind {
            SectionKind::Axisymmetric {} => "*SOLID SECTION, ELSET=all, MATERIAL=Steel\n,\n",
            SectionKind::PlaneStress { .. } => "*SOLID SECTION, ELSET=all, MATERIAL=Steel\n2.0\n",
            _ => "*SOLID SECTION, ELSET=all, MATERIAL=Steel\n1.0\n",
        };
        assert!(deck.contains(section), "{deck}");
    }
    let deck = abaqus::write(&plate(SectionKind::PlaneStrain { thickness: 1.0 }, ElementKind::Tri3, full)).unwrap();
    assert!(deck.contains("*ELEMENT, TYPE=CPE3\n") && deck.contains("3, S2\n"), "{deck}");
    // Axisymmetric body forces are radial and axial.
    let mut m = plate(SectionKind::Axisymmetric {}, ElementKind::Quad4, full);
    m.add_load("swell", Load::BodyForce { set: "all".into(), vector: [1.0, 2.0, 0.0] }).unwrap();
    m.set_load_active("Load", "swell", true).unwrap();
    let deck = abaqus::write(&m).unwrap();
    assert!(deck.contains("all, BR, 1.0\nall, BZ, 2.0\n"), "{deck}");
    // Abaqus has no reduced six-node triangle.
    let mut m = plate(SectionKind::PlaneStress { thickness: 1.0 }, ElementKind::Quad4, red);
    for (l, n) in [(10, 10), (11, 11), (12, 12), (13, 13), (14, 14), (15, 15)] {
        m.set_node(l, [1.0 + f64::from(n - 10) * 0.1, 2.0 + f64::from(n % 2), 0.0]).unwrap();
    }
    m.set_element(9, ElementKind::Tri6, vec![10, 12, 14, 11, 13, 15]).unwrap();
    m.replace_set("all", Set { elements: vec![1, 2, 9], ..Set::default() }).unwrap();
    assert!(abaqus::write(&m).is_err());
}

#[test]
fn nastran_writes_plane_stress_membranes_and_refuses_the_rest() {
    let m = plate(SectionKind::PlaneStress { thickness: 2.0 }, ElementKind::Quad4, Formulation::default());
    let deck = nastran::write(&m).unwrap();
    // A membrane: the material on MID1 only.
    assert!(deck.contains("PSHELL*                1               1              2.\n$"), "{deck}");
    assert!(deck.contains("CQUAD4         1       1       1       2       5       4"), "{deck}");
    // Every grid holds its out-of-plane degrees of freedom (PS on the continuation).
    let grids = deck.lines().filter(|l| l.starts_with("GRID")).count();
    assert_eq!(grids, 6);
    assert_eq!(deck.lines().filter(|l| l.starts_with('*') && l.trim_end().ends_with("3456")).count(), grids, "{deck}");
    // The edge pressure as consistent forces, pushing in (−X on the right edge).
    let forces: Vec<&str> = deck.lines().filter(|l| l.starts_with("FORCE")).collect();
    assert_eq!(forces.len(), 2, "{deck}");
    assert!(deck.contains("RFORCE"), "{deck}");
    for kind in [SectionKind::PlaneStrain { thickness: 1.0 }, SectionKind::Axisymmetric {}] {
        let why = unsupported(nastran::write(&plate(kind, ElementKind::Quad4, Formulation::default())));
        assert!(why.contains("differ between dialects"), "{why}");
    }
    let mut m = plate(SectionKind::PlaneStress { thickness: 2.0 }, ElementKind::Quad4, Formulation::default());
    m.replace_load("spin", Load::Centrifugal { angular_velocity: 1.0, point: [5.0, 0.0, 0.0], axis: [0.0, 0.0, 1.0], set: None }).unwrap();
    assert!(unsupported(nastran::write(&m)).contains("misses the origin"));
}

#[test]
fn lsdyna_writes_one_point_plane_forms_and_refuses_the_rest() {
    let red = Formulation { integration: Integration::Reduced, ..Formulation::default() };
    let mut m = plate(SectionKind::PlaneStress { thickness: 2.0 }, ElementKind::Quad4, red);
    m.remove_step("Load").unwrap();
    let mut s = Step::new("Load", Procedure::Static {});
    s.boundary_conditions.insert("hold".into());
    s.loads.insert("press".into());
    m.add_step(s, None).unwrap();
    let deck = lsdyna::write(&m).unwrap().text;
    assert!(deck.contains("*SECTION_SHELL") && deck.contains("*HOURGLASS"), "{deck}");
    assert!(deck.contains("*LOAD_NODE_POINT") && !deck.contains("*LOAD_SEGMENT"), "{deck}");
    let why = |m: &Model| match lsdyna::write(m) {
        Err(Error::Unsupported { what, .. }) => what,
        other => panic!("{other:?}"),
    };
    assert!(why(&plate(SectionKind::PlaneStress { thickness: 2.0 }, ElementKind::Quad4, Formulation::default())).contains("one point"));
    assert!(why(&plate(SectionKind::Axisymmetric {}, ElementKind::Tri3, red)).contains("axisymmetric"));
    assert!(why(&plate(SectionKind::PlaneStrain { thickness: 2.0 }, ElementKind::Tri3, red)).contains("per unit thickness"));
    assert!(why(&plate(SectionKind::PlaneStrain { thickness: 1.0 }, ElementKind::Tri3, red)).contains("centrifugal"));
}

#[test]
fn edge_pressures_are_the_consistent_nodal_forces() {
    // Pressure 5 on the right edge (x = 3, length 1): plane stress t = 2 gives −5 per node in X;
    // axisymmetric, the ring 2π · 3 · 1 shares 5 · 6π equally on a straight edge at fixed r.
    let m = plate(SectionKind::PlaneStress { thickness: 2.0 }, ElementKind::Quad4, Formulation::default());
    let r = m.resolve().unwrap();
    let f = edge_pressure_loads(&m, &r, ["press"]).unwrap();
    assert_eq!(f.len(), 2);
    assert!(f.values().all(|v| (v + 5.0).abs() < 1e-12), "{f:?}");
    let m = plate(SectionKind::Axisymmetric {}, ElementKind::Quad4, Formulation::default());
    let r = m.resolve().unwrap();
    let total: f64 = edge_pressure_loads(&m, &r, ["press"]).unwrap().values().sum();
    assert!((total + 5.0 * 6.0 * std::f64::consts::PI).abs() < 1e-10, "{total}");
}
