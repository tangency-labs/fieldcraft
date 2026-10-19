//! LS-DYNA deck writer: golden files of our own small models, node-order mappings checked by the
//! sign of the element Jacobian, a re-read of our own output, and what the format cannot take.
//!
//! Golden files live in `tests/lsdyna/`; `FIELDCRAFT_BLESS=1 cargo test -p fieldcraft-decks`
//! rewrites them after an intended change (review the diff).

use std::collections::{BTreeMap, BTreeSet};

mod common;

use common::lsdyna;
use fieldcraft_decks::Error;
use fieldcraft_model::*;

// ---------------------------------------------------------------------------------------------
// Our own models.

fn steel() -> Material {
    Material { elastic: Some(Elastic::Isotropic { young: 210000.0, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() }
}

fn step(name: &str, procedure: Procedure, bcs: &[&str], loads: &[&str]) -> Step {
    let mut s = Step::new(name, procedure);
    s.boundary_conditions = bcs.iter().map(|b| b.to_string()).collect();
    s.loads = loads.iter().map(|l| l.to_string()).collect();
    s
}

/// A 4 × 1 × 1 cantilever of four hex8 elements along x, fixed at x = 0. Nodes are numbered
/// x-fastest: label = 1 + i + 5·(j + 2·k).
fn cantilever_mesh(m: &mut Model) {
    let label = |i: u32, j: u32, k: u32| 1 + i + 5 * (j + 2 * k);
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..5 {
                m.set_node(label(i, j, k), [f64::from(i), f64::from(j), f64::from(k)]).unwrap();
            }
        }
    }
    for i in 0..4 {
        let nodes = vec![
            label(i, 0, 0),
            label(i + 1, 0, 0),
            label(i + 1, 1, 0),
            label(i, 1, 0),
            label(i, 0, 1),
            label(i + 1, 0, 1),
            label(i + 1, 1, 1),
            label(i, 1, 1),
        ];
        m.set_element(i + 1, ElementKind::Hex8, nodes).unwrap();
    }
    m.add_set("all", Set { elements: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    m.add_set("root", Set { nodes: vec![label(0, 0, 0), label(0, 1, 0), label(0, 0, 1), label(0, 1, 1)], ..Set::default() }).unwrap();
    // Face 4 of a hex8 is the side along edge 2-3: x = 4 for the last element.
    m.add_set("tip", Set { faces: vec![ElementFace { element: 4, face: 4 }], ..Set::default() }).unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
}

fn cantilever() -> Model {
    let mut m = Model::default();
    cantilever_mesh(&mut m);
    m.add_material("Steel", steel()).unwrap();
    m.add_section("beam", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.add_load("tip_load", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, -100.0], distribution: Distribution::Total }).unwrap();
    m.add_step(step("Bend", Procedure::Static {}, &["clamp"], &["tip_load"]), None).unwrap();
    m
}

/// The cantilever in an orthotropic material with incompatible modes, loaded through a coupled
/// reference point, then its natural frequencies: two steps, so two cases.
fn cantilever_modes() -> Model {
    let mut m = Model::default();
    cantilever_mesh(&mut m);
    let ortho = Elastic::EngineeringConstants { young: [140000.0, 10000.0, 10000.0], poisson: [0.3, 0.3, 0.45], shear: [5000.0, 5000.0, 3500.0] };
    m.add_material("Ply", Material { elastic: Some(ortho), density: Some(1.6e-9), ..Material::default() }).unwrap();
    let f = Formulation { incompatible_modes: true, ..Formulation::default() };
    m.add_section("beam", Section { set: "all".into(), material: "Ply".into(), kind: SectionKind::Solid {}, formulation: f }).unwrap();
    m.add_reference_point("RP", ReferencePoint { xyz: [5.0, 0.5, 0.5], label: 0 }).unwrap();
    m.add_coupling("pad", Coupling { reference_point: "RP".into(), set: "tip".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
    m.add_load("twist", Load::Moment { target: Target::ReferencePoint("RP".into()), vector: [250.0, 0.0, 0.0], distribution: Distribution::Total })
        .unwrap();
    m.add_load("push", Load::Force { target: Target::ReferencePoint("RP".into()), vector: [0.0, -10.0, 0.0], distribution: Distribution::Total })
        .unwrap();
    m.add_step(step("Bend", Procedure::Static {}, &["clamp"], &["push", "twist"]), None).unwrap();
    m.add_step(step("Modes", Procedure::Frequency { modes: 6, min_hz: Some(1.0), max_hz: Some(5000.0) }, &["clamp"], &[]), None).unwrap();
    m
}

/// A 2 × 1 plate at z = 0: two quad4 and two tri3 shells, held on the edge x = 0, pressed and
/// under gravity.
fn plate() -> Model {
    let mut m = Model::default();
    for (l, x, y) in
        [(1, 0.0, 0.0), (2, 10.0, 0.0), (3, 20.0, 0.0), (4, 0.0, 10.0), (5, 10.0, 10.0), (6, 20.0, 10.0), (7, 30.0, 0.0), (8, 30.0, 10.0)]
    {
        m.set_node(l, [x, y, 0.0]).unwrap();
    }
    m.set_element(1, ElementKind::Quad4, vec![1, 2, 5, 4]).unwrap();
    m.set_element(2, ElementKind::Quad4, vec![2, 3, 6, 5]).unwrap();
    m.set_element(3, ElementKind::Tri3, vec![3, 7, 8]).unwrap();
    m.set_element(4, ElementKind::Tri3, vec![3, 8, 6]).unwrap();
    m.add_set("plate", Set { elements: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    m.add_set("edge", Set { nodes: vec![1, 4], ..Set::default() }).unwrap();
    let top: Vec<ElementFace> = (1..=4).map(|element| ElementFace { element, face: 1 }).collect();
    m.add_set("top", Set { faces: top, ..Set::default() }).unwrap();
    m.add_material("Steel", steel()).unwrap();
    m.add_section(
        "skin",
        Section {
            set: "plate".into(),
            material: "Steel".into(),
            kind: SectionKind::Shell { thickness: 0.5, offset: 0.25, points: 5 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("edge".into()) }).unwrap();
    m.add_load("press", Load::Pressure { set: "top".into(), magnitude: 0.01 }).unwrap();
    m.add_load("weight", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: None }).unwrap();
    m.add_step(step("Sag", Procedure::Static {}, &["hold"], &["press", "weight"]), None).unwrap();
    m
}

/// One element of each solid kind besides hex8, side by side (each its own section), with a
/// distributing coupling, a prescribed displacement and a traction.
fn solids() -> Model {
    let mut m = Model::default();
    let mut next = 0u32;
    let mut node = |m: &mut Model, xyz: [f64; 3]| {
        next += 1;
        m.set_node(next, xyz).unwrap();
        next
    };
    // tet4
    let t: Vec<u32> = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]].iter().map(|p| node(&mut m, *p)).collect();
    m.set_element(1, ElementKind::Tet4, t.clone()).unwrap();
    // tet10: corners then mid-sides on 1-2, 2-3, 3-1, 1-4, 2-4, 3-4.
    let c = [[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0]];
    let mid = |a: usize, b: usize| [(c[a][0] + c[b][0]) / 2.0, (c[a][1] + c[b][1]) / 2.0, (c[a][2] + c[b][2]) / 2.0];
    let mut t10: Vec<u32> = c.iter().map(|p| node(&mut m, *p)).collect();
    for (a, b) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
        t10.push(node(&mut m, mid(a, b)));
    }
    m.set_element(2, ElementKind::Tet10, t10).unwrap();
    // wedge6
    let w: Vec<u32> = [[4.0, 0.0, 0.0], [5.0, 0.0, 0.0], [4.0, 1.0, 0.0], [4.0, 0.0, 1.0], [5.0, 0.0, 1.0], [4.0, 1.0, 1.0]]
        .iter()
        .map(|p| node(&mut m, *p))
        .collect();
    m.set_element(3, ElementKind::Wedge6, w).unwrap();
    // hex20
    let h = hex20_nodes([6.0, 0.0, 0.0]);
    let h: Vec<u32> = h.iter().map(|p| node(&mut m, *p)).collect();
    m.set_element(4, ElementKind::Hex20, h).unwrap();

    m.add_material(
        "Aluminium",
        Material { elastic: Some(Elastic::Isotropic { young: 70e9, poisson: 0.33 }), density: Some(2700.0), ..Material::default() },
    )
    .unwrap();
    for (e, name) in [(1, "tet4"), (2, "tet10"), (3, "wedge"), (4, "hex20")] {
        m.add_set(name, Set { elements: vec![e], ..Set::default() }).unwrap();
        let f = Formulation { hybrid: e == 1, ..Formulation::default() };
        m.add_section(name, Section { set: name.into(), material: "Aluminium".into(), kind: SectionKind::Solid {}, formulation: f }).unwrap();
    }
    m.add_set("base", Set { nodes: t.clone(), ..Set::default() }).unwrap();
    // Face 2 of the hex20 is its top (nodes 5-8 and their mid-sides).
    m.add_set("lid", Set { faces: vec![ElementFace { element: 4, face: 2 }], ..Set::default() }).unwrap();
    m.add_reference_point("RP", ReferencePoint { xyz: [6.5, 0.5, 2.0], label: 0 }).unwrap();
    m.add_coupling("spider", Coupling { reference_point: "RP".into(), set: "lid".into(), kind: CouplingKind::Distributing, dofs: Dofs::all() })
        .unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    m.add_boundary_condition(
        "lift",
        BoundaryCondition::Displacement { target: Target::Set("tet10".into()), values: [None, None, Some(0.001), None, None, None] },
    )
    .unwrap();
    m.add_load("drag", Load::Traction { set: "lid".into(), vector: [1000.0, 0.0, 0.0] }).unwrap();
    m.add_load("pull", Load::Force { target: Target::ReferencePoint("RP".into()), vector: [0.0, 0.0, 50.0], distribution: Distribution::Total })
        .unwrap();
    m.add_step(step("Load", Procedure::Static {}, &["hold", "lift"], &["drag", "pull"]), None).unwrap();
    m
}

/// A unit 20-node hexahedron at `o` in model order.
fn hex20_nodes(o: [f64; 3]) -> Vec<[f64; 3]> {
    let corners = CUBE.map(|c| c.map(|v| (v + 1.0) / 2.0));
    let mut out: Vec<[f64; 3]> = corners.to_vec();
    out.extend(HEX_EDGES.map(|(a, b)| [0, 1, 2].map(|k| (corners[a][k] + corners[b][k]) / 2.0)));
    out.iter().map(|p| [p[0] + o[0], p[1] + o[1], p[2] + o[2]]).collect()
}

// ---------------------------------------------------------------------------------------------
// Golden files.

fn golden(name: &str, model: &Model) -> lsdyna::Written {
    // The solids are in SI numbers; the others in mm-t-s.
    let write = |m: &Model| if name == "solids" { lsdyna::write_si(m) } else { lsdyna::write(m) };
    let written = write(model).unwrap();
    let again = write(model).unwrap();
    assert_eq!(written, again, "deterministic");
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/lsdyna").join(format!("{name}.k"));
    if std::env::var_os("FIELDCRAFT_BLESS").is_some() {
        std::fs::write(&path, &written.text).unwrap();
    }
    let expected = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e} (run with FIELDCRAFT_BLESS=1)", path.display()));
    assert!(written.text == expected, "{name}.k differs from the golden file (FIELDCRAFT_BLESS=1 rewrites it)");
    assert!(written.text.is_ascii());
    for line in written.text.lines() {
        assert!(line.len() <= 200 && line == line.trim_end(), "{line:?}");
        assert!(!line.contains('\t'));
    }
    assert_eq!(written.rounded_values, 0);
    written
}

#[test]
fn golden_hex8_cantilever() {
    let w = golden("cantilever", &cantilever());
    assert_eq!(w.cases, 0);
    let deck = read(&w.text);
    assert_eq!(deck.nodes.len(), 20);
    // 100 N on the tip face: 25 N on each of its four nodes, written once on the node set.
    assert!(w.text.contains("*LOAD_NODE_SET"));
    assert_eq!(deck.node_sets["tip"], vec![5, 10, 15, 20]);
    assert_eq!(deck.node_sets["root"], vec![1, 6, 11, 16]);
}

#[test]
fn golden_cantilever_modes_in_cases() {
    let w = golden("cantilever_modes", &cantilever_modes());
    assert_eq!(w.cases, 2);
    for k in ["*CASE_BEGIN_1", "*CASE_END_1", "*CASE_BEGIN_2", "*CASE_END_2", "*CONTROL_IMPLICIT_EIGENVALUE", "*MAT_ORTHOTROPIC_ELASTIC_TITLE"] {
        assert!(w.text.contains(k), "{k}");
    }
    // The rigid body's set leads with the reference point.
    let deck = read(&w.text);
    let rp = cantilever_modes().reference_point("RP").unwrap().label;
    assert_eq!(deck.node_sets["pad (rigid body)"].first(), Some(&rp));
    assert_eq!(deck.nodes[&rp], [5.0, 0.5, 0.5]);
}

#[test]
fn golden_plate_with_pressure_and_gravity() {
    let w = golden("plate", &plate());
    assert!(w.text.contains("*LOAD_SEGMENT_SET") && w.text.contains("*LOAD_BODY_Z"));
    let deck = read(&w.text);
    // Triangles are degenerate quadrilaterals.
    assert_eq!(deck.elements[&3], vec![3, 7, 8, 8]);
    assert_eq!(deck.parts.len(), 2, "quad4 and tri3 take different formulations");
}

#[test]
fn golden_solids_and_distributing_coupling() {
    let w = golden("solids", &solids());
    for k in ["*ELEMENT_SOLID_H20", "*CONSTRAINED_INTERPOLATION", "*BOUNDARY_PRESCRIBED_MOTION_SET_ID", "*LOAD_NODE_POINT"] {
        assert!(w.text.contains(k), "{k}");
    }
}

// ---------------------------------------------------------------------------------------------
// A minimal reader of our own output: nodes, elements and node sets.

#[derive(Default)]
struct Read {
    nodes: BTreeMap<u32, [f64; 3]>,
    /// Element → its nodes as written.
    elements: BTreeMap<u32, Vec<u32>>,
    /// The part ids elements are in.
    parts: BTreeSet<u32>,
    node_sets: BTreeMap<String, Vec<u32>>,
}

fn fields(line: &str) -> Vec<&str> {
    (0..line.len()).step_by(20).map(|i| line[i..(i + 20).min(line.len())].trim()).collect()
}

fn read(text: &str) -> Read {
    let mut out = Read::default();
    let mut keyword = String::new();
    let mut cards: Vec<&str> = Vec::new();
    let flush = |keyword: &str, cards: &[&str], out: &mut Read| match keyword {
        "*NODE" => {
            for c in cards {
                let f = fields(c);
                out.nodes.insert(f[0].parse().unwrap(), [f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap()]);
            }
        }
        "*ELEMENT_SOLID" | "*ELEMENT_SOLID_H20" => {
            let per = if keyword == "*ELEMENT_SOLID" { 2 } else { 3 };
            for e in cards.chunks(per) {
                let head = fields(e[0]);
                let nodes: Vec<u32> = e[1..].iter().flat_map(|c| fields(c)).map(|v| v.parse().unwrap()).collect();
                out.elements.insert(head[0].parse().unwrap(), nodes);
                out.parts.insert(head[1].parse().unwrap());
            }
        }
        "*ELEMENT_SHELL" | "*ELEMENT_BEAM" => {
            for c in cards {
                let f: Vec<u32> = fields(c).iter().map(|v| v.parse().unwrap()).collect();
                out.elements.insert(f[0], f[2..].to_vec());
                out.parts.insert(f[1]);
            }
        }
        "*SET_NODE_LIST_TITLE" => {
            let nodes = cards[2..].iter().flat_map(|c| fields(c)).map(|v| v.parse().unwrap()).collect();
            out.node_sets.insert(cards[0].to_string(), nodes);
        }
        _ => {}
    };
    for line in text.lines() {
        if line.starts_with('$') {
            continue;
        }
        if line.starts_with('*') {
            flush(&keyword, &cards, &mut out);
            keyword = line.to_string();
            cards.clear();
        } else {
            cards.push(line);
        }
    }
    flush(&keyword, &cards, &mut out);
    out
}

#[test]
fn re_read_gives_back_nodes_elements_and_sets() {
    for model in [cantilever(), cantilever_modes(), plate(), solids()] {
        let w = lsdyna::write(&model).unwrap();
        let deck = read(&w.text);
        for (label, n) in model.nodes() {
            assert_eq!(deck.nodes[label], n.xyz, "node {label} round-trips exactly");
        }
        for (label, e) in model.elements() {
            let written = &deck.elements[label];
            let unique: BTreeSet<u32> = written.iter().copied().collect();
            assert_eq!(unique, e.nodes.iter().copied().collect::<BTreeSet<u32>>(), "element {label}");
        }
        let r = model.resolve().unwrap();
        for (name, nodes) in &deck.node_sets {
            if let Some(s) = r.sets.get(name) {
                assert_eq!(nodes, &s.nodes, "set {name}");
            }
        }
    }
}

#[test]
fn awkward_coordinates_round_trip_or_are_counted() {
    let mut m = cantilever();
    // Full 17-digit mantissas at every magnitude a model allows.
    let values = [0.1 + 0.2, -123_456_789.123_456_78, 1e-300, -2.2250738585072014e-308, 0.012_345_678_901_234_567, 999_999_999.999_999_9];
    for (i, v) in values.iter().enumerate() {
        m.set_node(100 + i as u32, [*v, -*v, *v / 3.0]).unwrap();
    }
    let w = lsdyna::write_si(&m).unwrap();
    let deck = read(&w.text);
    let mut inexact = 0;
    for (i, _) in values.iter().enumerate() {
        let label = 100 + i as u32;
        let want = m.node(label).unwrap().xyz;
        let got = deck.nodes[&label];
        for (a, b) in want.iter().zip(got) {
            if *a != b {
                inexact += 1;
                assert!(((a - b) / a).abs() < 5e-13, "{a:e} {b:e}");
            }
        }
    }
    assert_eq!(inexact, w.rounded_values);
}

// ---------------------------------------------------------------------------------------------
// Node order: the sign of the Jacobian in LS-DYNA's element conventions.

/// Hexahedron corners in the reference cube, in LS-DYNA (and model) order.
const CUBE: [[f64; 3]; 8] =
    [[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.], [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]];

/// Hexahedron edges as corner pairs, in mid-side node order.
const HEX_EDGES: [(usize, usize); 12] = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];

/// det J of the map with shape-function gradients `dn[i]` (reference axes) at nodes `x[i]`.
fn jacobian_det(x: &[[f64; 3]], dn: impl Iterator<Item = [f64; 3]>) -> f64 {
    let mut j = [[0.0; 3]; 3];
    for (p, d) in x.iter().zip(dn) {
        for a in 0..3 {
            for b in 0..3 {
                j[a][b] += d[a] * p[b];
            }
        }
    }
    det(j)
}

/// Trilinear hexahedron in LS-DYNA order (1-4 one face, 5-8 over them): det J at a point.
fn hex8_det(x: &[[f64; 3]], (r, s, t): (f64, f64, f64)) -> f64 {
    jacobian_det(
        x,
        CUBE.iter().map(|c| {
            [
                c[0] * (1.0 + c[1] * s) * (1.0 + c[2] * t) / 8.0,
                c[1] * (1.0 + c[0] * r) * (1.0 + c[2] * t) / 8.0,
                c[2] * (1.0 + c[0] * r) * (1.0 + c[1] * s) / 8.0,
            ]
        }),
    )
}

fn det(j: [[f64; 3]; 3]) -> f64 {
    j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
        + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0])
}

const GAUSS: f64 = 0.577_350_269_189_625_8;

/// Write a one-element model of `kind` with these model-order coordinates and read back the
/// written coordinates in LS-DYNA order.
fn written_coordinates(kind: ElementKind, xyz: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let mut m = Model::default();
    for (i, p) in xyz.iter().enumerate() {
        m.set_node(i as u32 + 1, *p).unwrap();
    }
    m.set_element(1, kind, (1..=xyz.len() as u32).collect()).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_set("fix", Set { nodes: vec![1], ..Set::default() }).unwrap();
    m.add_material("M", steel()).unwrap();
    let kind_of = if kind.dimension() == 3 { SectionKind::Solid {} } else { SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 3 } };
    m.add_section("s", Section { set: "all".into(), material: "M".into(), kind: kind_of, formulation: Formulation::default() }).unwrap();
    m.add_boundary_condition("b", BoundaryCondition::Fixed { target: Target::Set("fix".into()) }).unwrap();
    m.add_step(step("S", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None }, &["b"], &[]), None).unwrap();
    let deck = read(&lsdyna::write_si(&m).unwrap().text);
    deck.elements[&1].iter().map(|n| deck.nodes[n]).collect()
}

/// Slightly distorted coordinates, so a wrong mapping cannot pass by symmetry (short decimals,
/// so every coordinate is written exactly).
fn jiggle(p: [f64; 3], i: usize) -> [f64; 3] {
    let k = ((i as f64 * 0.37).sin() * 0.03 * 1e4).round() / 1e4;
    [p[0] + k, p[1] - k * 0.5, p[2] + k * 0.25].map(|c| (c * 1e6).round() / 1e6)
}

#[test]
fn degenerate_solids_have_positive_jacobians_in_lsdyna_order() {
    let tet = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let wedge = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]];
    let hex =
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]];
    for (kind, xyz) in [(ElementKind::Tet4, &tet[..]), (ElementKind::Wedge6, &wedge[..]), (ElementKind::Hex8, &hex[..])] {
        let xyz: Vec<[f64; 3]> = xyz.iter().enumerate().map(|(i, p)| jiggle(*p, i)).collect();
        let x = written_coordinates(kind, &xyz);
        assert_eq!(x.len(), 8, "{}", kind.name());
        for r in [-GAUSS, GAUSS] {
            for s in [-GAUSS, GAUSS] {
                for t in [-GAUSS, GAUSS] {
                    assert!(hex8_det(&x, (r, s, t)) > 0.0, "{} at ({r}, {s}, {t})", kind.name());
                }
            }
        }
        // The model's own order mirrored is negative (the test can fail).
        let mut mirrored = xyz.clone();
        for p in &mut mirrored {
            p[2] = -p[2];
        }
        let x = written_coordinates_unchecked(kind, &mirrored);
        assert!(hex8_det(&x, (0.1, -0.1, 0.1)) < 0.0, "{}", kind.name());
    }
}

/// [`written_coordinates`] without writing: the writer's node order applied to coordinates the
/// model would reject as inverted.
fn written_coordinates_unchecked(kind: ElementKind, xyz: &[[f64; 3]]) -> Vec<[f64; 3]> {
    lsdyna::node_order(kind).iter().map(|&i| xyz[i]).collect()
}

/// Quadratic tetrahedron in LS-DYNA order (mid-sides 5, 6, 7 on 1-2, 2-3, 3-1 and 8, 9, 10 on
/// 2-4, 3-4, 1-4): det J at volume coordinates (l1, l2, l3).
fn tet10_det(x: &[[f64; 3]], l: [f64; 3]) -> f64 {
    let ls = [1.0 - l[0] - l[1] - l[2], l[0], l[1], l[2]];
    // dL/d(ξ, η, ζ) with L1 = 1 − ξ − η − ζ, L2 = ξ, L3 = η, L4 = ζ.
    let dl = [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let edges = [(0, 1), (1, 2), (2, 0), (1, 3), (2, 3), (0, 3)];
    let corners = (0..4).map(|i| dl[i].map(|d| (4.0 * ls[i] - 1.0) * d));
    let mids = edges.iter().map(|&(a, b)| [0, 1, 2].map(|c| 4.0 * (dl[a][c] * ls[b] + ls[a] * dl[b][c])));
    jacobian_det(x, corners.chain(mids))
}

#[test]
fn tet10_mid_side_nodes_sit_on_lsdyna_edges_with_positive_jacobian() {
    let c = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut xyz: Vec<[f64; 3]> = c.iter().enumerate().map(|(i, p)| jiggle(*p, i)).collect();
    for (a, b) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
        xyz.push([0, 1, 2].map(|k| (xyz[a][k] + xyz[b][k]) / 2.0));
    }
    let x = written_coordinates(ElementKind::Tet10, &xyz);
    assert_eq!(x.len(), 10);
    for (m, (a, b)) in [(0, 1), (1, 2), (2, 0), (1, 3), (2, 3), (0, 3)].iter().enumerate() {
        let want = [0, 1, 2].map(|k| (x[*a][k] + x[*b][k]) / 2.0);
        assert_eq!(x[4 + m], want, "node {} is on LS-DYNA edge {}-{}", 5 + m, a + 1, b + 1);
    }
    for l in [[0.25, 0.25, 0.25], [0.1, 0.1, 0.1], [0.6, 0.1, 0.1], [0.1, 0.6, 0.1], [0.1, 0.1, 0.6]] {
        assert!(tet10_det(&x, l) > 0.0, "{l:?}");
    }
}

/// 20-node serendipity hexahedron in LS-DYNA order: det J at a point.
fn hex20_det(x: &[[f64; 3]], p: [f64; 3]) -> f64 {
    let mut nodes: Vec<[f64; 3]> = CUBE.to_vec();
    nodes.extend(HEX_EDGES.map(|(a, b)| [0, 1, 2].map(|k| (CUBE[a][k] + CUBE[b][k]) / 2.0)));
    let gradient = |n: &[f64; 3]| -> [f64; 3] {
        let q = [0, 1, 2].map(|k| 1.0 + n[k] * p[k]);
        match (0..3).find(|k| n[*k] == 0.0) {
            // Corner: N = q0 q1 q2 (n·p − 2) / 8.
            None => {
                let s = n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - 2.0;
                [0, 1, 2].map(|k| n[k] * q[(k + 1) % 3] * q[(k + 2) % 3] * (s + q[k]) / 8.0)
            }
            // Mid-side on an edge along axis z: N = (1 − p_z²) q_a q_b / 4.
            Some(z) => [0, 1, 2].map(|k| {
                let others = |skip: usize| (0..3).filter(|m| *m != skip && *m != z).map(|m| q[m]).product::<f64>();
                if k == z { -2.0 * p[z] * others(z) / 4.0 } else { n[k] * (1.0 - p[z] * p[z]) * others(k) / 4.0 }
            }),
        }
    };
    jacobian_det(x, nodes.iter().map(gradient))
}

#[test]
fn hex20_keeps_the_model_order_with_positive_jacobian() {
    let xyz: Vec<[f64; 3]> = hex20_nodes([0.0; 3]).into_iter().enumerate().map(|(i, p)| jiggle(p, i)).collect();
    let x = written_coordinates(ElementKind::Hex20, &xyz);
    assert_eq!(x, xyz);
    for r in [-GAUSS, GAUSS] {
        for s in [-GAUSS, GAUSS] {
            for t in [-GAUSS, GAUSS] {
                assert!(hex20_det(&x, [r, s, t]) > 0.0);
            }
        }
    }
    // The unit cube's Jacobian is 1/8 everywhere: the shape functions are right.
    assert!((hex20_det(&hex20_nodes([0.0; 3]), [0.3, -0.2, 0.7]) - 0.125).abs() < 1e-12);
}

#[test]
fn shells_keep_their_normal() {
    let normal = |x: &[[f64; 3]]| {
        // Diagonals of the (possibly degenerate) quadrilateral 1-2-3-4.
        let d1 = [0, 1, 2].map(|k| x[2][k] - x[0][k]);
        let d2 = [0, 1, 2].map(|k| x[3][k] - x[1][k]);
        [d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0]]
    };
    let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let quad = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
    let mid = |p: &[[f64; 3]], a: usize, b: usize| [0, 1, 2].map(|k| (p[a][k] + p[b][k]) / 2.0);
    let tri6: Vec<[f64; 3]> = tri.iter().copied().chain([mid(&tri, 0, 1), mid(&tri, 1, 2), mid(&tri, 2, 0)]).collect();
    let quad8: Vec<[f64; 3]> = quad.iter().copied().chain([mid(&quad, 0, 1), mid(&quad, 1, 2), mid(&quad, 2, 3), mid(&quad, 3, 0)]).collect();
    for (kind, xyz) in
        [(ElementKind::Tri3, tri.to_vec()), (ElementKind::Quad4, quad.to_vec()), (ElementKind::Tri6, tri6), (ElementKind::Quad8, quad8)]
    {
        let x = written_coordinates(kind, &xyz);
        assert!(normal(&x)[2] > 0.0, "{} keeps +z as its top", kind.name());
        if kind == ElementKind::Tri6 {
            assert_eq!(&x[4..7], &xyz[3..6]);
        }
        if kind == ElementKind::Quad8 {
            assert_eq!(&x[4..8], &xyz[4..8]);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// What the format cannot take, and hostile models.

fn unsupported(m: &Model, needle: &str) {
    match lsdyna::write(m) {
        Err(Error::Unsupported { format, what }) => {
            assert_eq!(format, "LS-DYNA");
            assert!(what.contains(needle), "{what}");
        }
        other => panic!("expected unsupported `{needle}`, got {other:?}"),
    }
}

#[test]
fn what_lsdyna_cannot_take_is_named() {
    let mut m = cantilever();
    m.add_load("heat", Load::BodyForce { set: "all".into(), vector: [0.0, 0.0, 1.0] }).unwrap();
    m.replace_step(step("Bend", Procedure::Static {}, &["clamp"], &["tip_load", "heat"])).unwrap();
    unsupported(&m, "load `heat`");

    let mut m = cantilever();
    m.replace_section(
        "beam",
        Section {
            set: "all".into(),
            material: "Steel".into(),
            kind: SectionKind::Solid {},
            formulation: Formulation { hybrid: true, ..Formulation::default() },
        },
    )
    .unwrap();
    unsupported(&m, "hybrid");

    let mut m = cantilever();
    m.replace_material("Steel", Material { density: None, ..steel() }).unwrap();
    unsupported(&m, "density");

    let mut m = cantilever_modes();
    m.replace_coupling(
        "pad",
        Coupling { reference_point: "RP".into(), set: "tip".into(), kind: CouplingKind::Kinematic, dofs: Dofs::try_from(vec![1, 2, 3]).unwrap() },
    )
    .unwrap();
    unsupported(&m, "coupling `pad`");

    let mut m = cantilever_modes();
    m.add_boundary_condition("hold_rp", BoundaryCondition::Fixed { target: Target::ReferencePoint("RP".into()) }).unwrap();
    m.replace_step(step("Bend", Procedure::Static {}, &["clamp", "hold_rp"], &["push"])).unwrap();
    unsupported(&m, "boundary condition `hold_rp`");

    let mut m = plate();
    m.replace_section(
        "skin",
        Section {
            set: "plate".into(),
            material: "Steel".into(),
            kind: SectionKind::Shell { thickness: 0.5, offset: 0.0, points: 11 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    unsupported(&m, "11 integration points");

    let mut m = plate();
    m.add_set("half", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.replace_load("weight", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: Some("half".into()) }).unwrap();
    unsupported(&m, "gravity on set `half`");
}

#[test]
fn gravity_on_whole_sections_uses_a_part_set() {
    let mut m = plate();
    m.replace_load("weight", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: Some("plate".into()) }).unwrap();
    let w = lsdyna::write(&m).unwrap();
    assert!(w.text.contains("*SET_PART_LIST") && w.text.contains("*LOAD_BODY_PARTS"));
}

#[test]
fn models_that_are_not_ready_are_refused_not_written() {
    assert!(matches!(lsdyna::write(&Model::default()), Err(Error::Model(_))));
    let mut m = cantilever();
    m.remove_step("Bend").unwrap();
    assert!(matches!(lsdyna::write(&m), Err(Error::Model(_))));
}

#[test]
fn hostile_numbers_stay_within_their_fields() {
    let mut m = cantilever();
    m.set_node(1, [1e9, -1e9, 1e-300]).unwrap();
    m.set_node(2, [f64::MIN_POSITIVE, -f64::MIN_POSITIVE, 5e-324]).unwrap();
    m.replace_material(
        "Steel",
        Material {
            elastic: Some(Elastic::Isotropic { young: 1e300, poisson: 0.499_999_999_999_999_9 }),
            density: Some(1e-300),
            ..Material::default()
        },
    )
    .unwrap();
    m.replace_load(
        "tip_load",
        Load::Force { target: Target::Set("tip".into()), vector: [1e-300, -1e300, 1.0 / 3.0], distribution: Distribution::PerNode },
    )
    .unwrap();
    let w = lsdyna::write(&m).unwrap();
    for line in w.text.lines().filter(|l| !l.starts_with('$') && !l.starts_with('*')) {
        for f in fields(line) {
            assert!(f.len() <= 20, "{line}");
        }
    }
    let deck = read(&w.text);
    assert_eq!(deck.nodes[&1], [1e9, -1e9, 1e-300]);
}

/// Explicit dynamic steps are not written by any deck format yet: each says so by name.
#[test]
fn explicit_steps_are_refused_by_every_format() {
    let mut m = cantilever();
    m.add_step(step("Impact", Procedure::ExplicitDynamic(ExplicitDynamic::new(1e-3)), &["clamp"], &[]), None).unwrap();
    for f in fieldcraft_decks::FORMATS {
        let e = f.write(&common::si(&m), fieldcraft_model::units::Preset::MmTonneS).unwrap_err();
        assert!(matches!(e, Error::Unsupported { .. }), "{}: {e}", f.name);
        assert!(e.to_string().contains("step `Impact`: explicit dynamic"), "{}: {e}", f.name);
    }
}

/// Plasticity: `*MAT_PLASTIC_KINEMATIC` (bilinear) and `*MAT_PIECEWISE_LINEAR_PLASTICITY` with a
/// curve numbered after the steps' ramps; the erosion strain is the failure strain; Johnson–Cook
/// and damage are reported.
#[test]
fn plasticity_is_written_as_plastic_materials() {
    let with = |p: Plastic, erosion: Option<f64>| {
        let mut m = cantilever();
        m.replace_material("Steel", Material { plastic: Some(p), erosion: erosion.map(|e| Erosion { plastic_strain: e }), ..steel() }).unwrap();
        m
    };
    let w = lsdyna::write(&with(Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 1000.0, hardening: Hardening::Isotropic }, Some(0.5)))
        .unwrap();
    let mat = w.text.split("*MAT_PLASTIC_KINEMATIC_TITLE\nSteel\n").nth(1).unwrap();
    let cards: Vec<Vec<&str>> = mat.lines().filter(|l| !l.starts_with('$')).take(2).map(|l| l.split_whitespace().collect()).collect();
    // E_t = E H / (E + H) = 210000 · 1000 / 211000.
    let etan: f64 = cards[0][5].parse().unwrap();
    assert!((etan - 210000.0 * 1000.0 / 211000.0).abs() < 1e-9, "{cards:?}");
    assert_eq!((cards[0][4], cards[0][6], cards[1][0]), ("250.0", "1.0", "0.5"), "{cards:?}");
    assert!(!w.text.contains("*MAT_ELASTIC"));
    let table = Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: 200.0 }, HardeningPoint { strain: 0.1, stress: 300.0 }] };
    let w = lsdyna::write(&with(table, None)).unwrap();
    let mat = w.text.split("*MAT_PIECEWISE_LINEAR_PLASTICITY_TITLE\nSteel\n").nth(1).unwrap();
    let cards: Vec<&str> = mat.lines().filter(|l| !l.starts_with('$')).take(2).collect();
    assert_eq!(cards[1].split_whitespace().collect::<Vec<_>>(), ["2"], "the curve after the one step's ramp: {cards:?}");
    let curve = w.text.split("*DEFINE_CURVE_TITLE\nSteel hardening\n").nth(1).unwrap();
    let rows: Vec<Vec<&str>> = curve.lines().filter(|l| !l.starts_with('$')).take(3).map(|l| l.split_whitespace().collect()).collect();
    assert_eq!(rows, [vec!["2"], vec!["0.0", "200.0"], vec!["0.1", "300.0"]]);
    let jc = Plastic::JohnsonCook {
        a: 90.0,
        b: 292.0,
        n: 0.31,
        c: 0.0,
        m: 1.0,
        reference_strain_rate: 1.0,
        room_temperature: 294.0,
        melt_temperature: 1356.0,
        specific_heat: None,
        taylor_quinney: 0.9,
    };
    let Err(Error::Unsupported { what, .. }) = lsdyna::write(&with(jc, None)) else { panic!() };
    assert!(what.contains("Johnson–Cook"), "{what}");
}

/// A nonlinear static step: the ramped implicit run with full Newton and automatic increments;
/// what LS-DYNA cannot take is named.
#[test]
fn nonlinear_static_steps() {
    let mut m = cantilever();
    let controls =
        NonlinearStatic { initial_increment: 0.05, min_increment: 1e-4, max_increment: 0.5, max_iterations: 12, ..NonlinearStatic::default() };
    m.replace_step(step("Bend", Procedure::NonlinearStatic(controls.clone()), &["clamp"], &["tip_load"])).unwrap();
    let deck = lsdyna::write(&m).unwrap().text;
    let after = |keyword: &str| -> Vec<String> {
        let mut lines = deck.lines().skip_while(|l| *l != keyword).skip(1).filter(|l| !l.starts_with('$'));
        lines.next().unwrap().split_whitespace().map(str::to_string).collect()
    };
    let num = |s: &String| s.parse::<f64>().unwrap();
    let general = after("*CONTROL_IMPLICIT_GENERAL");
    assert_eq!((num(&general[0]), num(&general[1])), (1.0, 0.05), "{general:?}");
    let solution = after("*CONTROL_IMPLICIT_SOLUTION");
    assert_eq!(solution.iter().map(num).collect::<Vec<_>>(), [2.0, 1.0, 12.0], "{solution:?}");
    let auto = after("*CONTROL_IMPLICIT_AUTO");
    assert_eq!(auto.iter().map(num).collect::<Vec<_>>(), [1.0, 1e-4, 0.5], "blank fields are skipped: {auto:?}");
    assert!(deck.contains("*DEFINE_CURVE_TITLE\nBend ramp\n"), "the loads follow the ramp");
    assert!(deck.contains("*LOAD_NODE"), "{deck}");
    // A load path is the curve.
    let path = NonlinearStatic { load_path: Some(vec![[0.0, 0.0], [0.5, 1.0], [1.0, -1.0]]), ..controls.clone() };
    m.replace_step(step("Bend", Procedure::NonlinearStatic(path), &["clamp"], &["tip_load"])).unwrap();
    let deck = lsdyna::write(&m).unwrap().text;
    let curve: Vec<Vec<f64>> = deck
        .lines()
        .skip_while(|l| *l != "Bend ramp")
        .filter(|l| !l.starts_with('$'))
        .skip(2)
        .take(3)
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect();
    assert_eq!(curve, [vec![0.0, 0.0], vec![0.5, 1.0], vec![1.0, -1.0]], "{deck}");

    m.replace_step(step("Bend", Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: false, ..controls }), &["clamp"], &["tip_load"]))
        .unwrap();
    unsupported(&m, "a nonlinear step without nonlinear geometry");
    m.replace_step(step("Bend", Procedure::NonlinearStatic(NonlinearStatic::default()), &["clamp"], &["tip_load"])).unwrap();
    m.replace_material("Steel", Material { hyperelastic: Some(Hyperelastic::NeoHookean {}), ..steel() }).unwrap();
    unsupported(&m, "hyperelastic material `Steel` in a nonlinear step");
    // In a linear step the elastic constants are the whole of a hyperelastic material.
    m.replace_step(step("Bend", Procedure::Static {}, &["clamp"], &["tip_load"])).unwrap();
    assert!(lsdyna::write(&m).unwrap().text.contains("*MAT_ELASTIC_TITLE"));
}
