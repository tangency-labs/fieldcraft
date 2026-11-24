//! 2D models (plane stress, plane strain, axisymmetric) in every deck format: element names,
//! section data, edge surfaces and pressures, centrifugal loads, and what a format cannot write.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::{Error, edge_pressure_loads};
use fieldcraft_model::{
    BeamProfile, BeamTheory, BoundaryCondition, Elastic, ElementEdge, ElementKind, Formulation, Integration, Load, Material, Model, Procedure,
    Section, SectionKind, Set, Step, Symmetry, Target,
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

/// The plane stress plate with a beam along its bottom edge (nodes 1, 2, 3), the left edge held
/// and the right edge pressed.
fn stiffened(formulation: Formulation) -> Model {
    let mut m = plate(SectionKind::PlaneStress { thickness: 2.0 }, ElementKind::Quad4, formulation);
    m.set_element(10, ElementKind::Beam2, vec![1, 2]).unwrap();
    m.set_element(11, ElementKind::Beam2, vec![2, 3]).unwrap();
    m.add_set("rib", Set { elements: vec![10, 11], ..Set::default() }).unwrap();
    let kind = SectionKind::Beam { profile: BeamProfile::Rectangle { width: 2.0, height: 4.0 }, orientation: [0.0, 0.0, 1.0] };
    m.add_section("rib", Section { set: "rib".into(), material: "Steel".into(), kind, formulation: Formulation::default() }).unwrap();
    m.remove_step("Load").unwrap();
    let mut s = Step::new("Load", Procedure::Static {});
    s.boundary_conditions.insert("hold".into());
    s.loads.insert("press".into());
    m.add_step(s, None).unwrap();
    m
}

#[test]
fn beams_of_a_2d_model_are_planar_beams_held_in_plane_or_refused_by_name() {
    use fieldcraft_model::units::Preset;
    let m = stiffened(Formulation::default());
    // Abaqus: the planar beam, its first section axis out of the plane.
    let deck = abaqus::write(&m).unwrap();
    assert!(deck.contains("*ELEMENT, TYPE=B21\n10, 1, 2\n11, 2, 3\n"), "{deck}");
    assert!(deck.contains("*BEAM SECTION, SECTION=RECT, ELSET=rib, MATERIAL=Steel\n2.0, 4.0\n0.0, 0.0, -1.0\n"), "{deck}");
    assert!(!deck.contains("B31"), "{deck}");
    // And back: a beam section along Z on the same elements, in a model that is still 2D.
    let read = fieldcraft_decks::abaqus::read_in(&deck, Preset::MmTonneS).unwrap();
    let beams: Vec<_> = read.model.sections().values().filter(|s| matches!(s.kind, SectionKind::Beam { .. })).collect();
    assert_eq!(beams.len(), 1);
    assert!(matches!(beams[0].kind, SectionKind::Beam { orientation, .. } if orientation == [0.0, 0.0, -1.0]), "{:?}", beams[0]);
    assert_eq!(read.model.element(10).map(|e| e.kind), Some(ElementKind::Beam2));
    // A planar beam without a 2D continuum, and a spatial beam with one, are refused on import.
    let frame = "*NODE\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=B21,ELSET=Bar\n1,1,2\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*BEAM SECTION,ELSET=Bar,MATERIAL=Steel,SECTION=CIRC\n0.1\n0,0,-1\n";
    let e = fieldcraft_decks::abaqus::read_in(frame, Preset::MmTonneS).unwrap_err().to_string();
    assert!(e.contains("B21 beams without a planar continuum"), "{e}");
    let e = fieldcraft_decks::abaqus::read_in(&deck.replace("TYPE=B21", "TYPE=B31"), Preset::MmTonneS).unwrap_err().to_string();
    assert!(e.contains("B31 and B33 beams with a planar continuum"), "{e}");

    // Nastran: CBEAM, its three grids holding 3, 4, 5 and the plate's other three 3–6.
    let deck = nastran::write(&m).unwrap();
    assert_eq!(deck.lines().filter(|l| l.starts_with("CBEAM")).count(), 2, "{deck}");
    let held = |ps: &str| deck.lines().filter(|l| l.starts_with('*') && l.trim_end().ends_with(ps)).count();
    assert_eq!((held(" 345"), held("3456")), (3, 3), "{deck}");
    // The beam is along X with its orientation along Z: the orientation vector v = axis × orientation
    // = (0, -1, 0), and PBEAML's BAR has the profile's width and height (2 and 4).
    let card = |name: &str| -> Vec<String> {
        let mut fields = Vec::new();
        let mut lines = deck.lines().skip_while(|l| !l.starts_with(name));
        if let Some(first) = lines.next() {
            for line in std::iter::once(first).chain(lines.take_while(|l| l.starts_with('*'))) {
                fields.extend(line.as_bytes().get(8..).unwrap_or_default().chunks(16).map(|c| String::from_utf8_lossy(c).trim().to_string()));
            }
        }
        fields
    };
    let real = |s: &str| -> f64 { s.parse().unwrap_or_else(|_| panic!("`{s}` is not a real\n{deck}")) };
    let cbeam = card("CBEAM");
    assert_eq!(cbeam.get(4..7).map(|v| v.iter().map(|s| real(s)).collect::<Vec<_>>()), Some(vec![0.0, -1.0, 0.0]), "{cbeam:?}\n{deck}");
    let pbeaml = card("PBEAML");
    assert_eq!(pbeaml.get(3).map(String::as_str), Some("BAR"), "{pbeaml:?}");
    assert_eq!(pbeaml.get(4..6).map(|v| v.iter().map(|s| real(s)).collect::<Vec<_>>()), Some(vec![2.0, 4.0]), "{pbeaml:?}\n{deck}");

    // An Euler–Bernoulli beam and a pipe pressure in a 2D model are refused before any deck is
    // written: the planar beam is a Timoshenko beam with no pipe wall.
    let shear_rigid = Formulation { beam_theory: BeamTheory::EulerBernoulli, ..Formulation::default() };
    let mut eb = stiffened(Formulation::default());
    let mut rib = eb.section("rib").unwrap().clone();
    rib.formulation = shear_rigid;
    eb.replace_section("rib", rib).unwrap();
    let mut piped = stiffened(Formulation::default());
    piped.add_load("steam", Load::PipePressure { set: "rib".into(), internal: 1.0, external: 0.0 }).unwrap();
    for (model, want) in
        [(&eb, "section `rib`: an Euler-Bernoulli beam is not supported in a 2D model"), (&piped, "a pipe pressure is not supported in a 2D model")]
    {
        for e in [abaqus::write(model).unwrap_err().to_string(), nastran::write(model).unwrap_err().to_string()] {
            assert!(e.contains(want), "{want}: {e}");
        }
    }

    // LS-DYNA: its 2D analyses take no beam; the refusal names the section.
    let red = Formulation { integration: Integration::Reduced, ..Formulation::default() };
    match lsdyna::write(&stiffened(red)) {
        Err(Error::Unsupported { what, .. }) => assert!(what.contains("section `rib`: beams in a 2D model"), "{what}"),
        other => panic!("{other:?}"),
    }
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

/// An inclined symmetry plane through a planar beam's node: the equations hold the translation
/// along the normal and the rotation about Z (an in-plane axis of the plane), and have no term on
/// the rotations about X and Y, which a 2D model's node lacks (Abaqus B21 has no dofs 4 and 5;
/// Nastran holds them in PS = 345).
#[test]
fn inclined_symmetry_at_a_planar_beam_node_has_no_terms_on_rotations_about_x_and_y() {
    let mut m = stiffened(Formulation::default());
    m.add_set("mid", Set { nodes: vec![2], ..Set::default() }).unwrap();
    let plane = Symmetry::from_normal([1.0, -1.0, 0.0]).unwrap();
    m.add_boundary_condition("sym", BoundaryCondition::Symmetry { target: Target::Set("mid".into()), plane }).unwrap();
    m.set_boundary_condition_active("Load", "sym", true).unwrap();

    let deck = abaqus::write(&m).unwrap();
    let lines: Vec<&str> = deck.lines().collect();
    let mut dofs = Vec::new();
    for (i, _) in lines.iter().enumerate().filter(|(_, l)| **l == "*EQUATION") {
        let f: Vec<&str> = lines[i + 2].split(',').map(str::trim).collect();
        dofs.push(f.chunks(3).map(|t| t[1].parse::<u8>().unwrap()).collect::<Vec<_>>());
    }
    assert!(!dofs.is_empty() && dofs.iter().flatten().all(|d| !matches!(d, 4 | 5)), "{dofs:?}\n{deck}");
    assert!(dofs.iter().flatten().any(|d| *d == 6), "{dofs:?}");

    // Nastran: the MPC entries' components, the second field of each (grid, component, value)
    // triple of the large-field lines (16 columns per field after the 8-column name).
    let deck = nastran::write(&m).unwrap();
    let mut components = Vec::new();
    let mut in_mpc = false;
    for line in deck.lines() {
        let (head, rest) = line.split_at(8.min(line.len()));
        match head.trim() {
            "MPC*" => in_mpc = true,
            "*" if in_mpc => {}
            _ => in_mpc = false,
        }
        if in_mpc {
            let fields: Vec<String> = rest.as_bytes().chunks(16).map(|c| String::from_utf8_lossy(c).trim().to_string()).collect();
            // Fields: [SID or blank, grid, component, value] then [blank-or-grid, ...]; take the
            // small integers that follow a grid label.
            components.extend(fields.iter().filter(|f| matches!(f.as_str(), "4" | "5" | "6")).cloned());
        }
    }
    assert!(!components.is_empty(), "{deck}");
    assert!(components.iter().all(|c| c == "6"), "{components:?}\n{deck}");
}
