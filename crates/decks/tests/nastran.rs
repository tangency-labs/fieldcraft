//! The Nastran writer on our own small models: golden files, a minimal re-reader of our own
//! output (grids, elements, constraint sets, case control), node order checked against Nastran's
//! own element conventions by Jacobian sign and mid-side positions, and what is unsupported.
//!
//! Golden files are rewritten with `FIELDCRAFT_BLESS=1 cargo test -p fieldcraft-decks`.

use std::collections::BTreeMap;

mod common;

use common::nastran;
use fieldcraft_decks::{Error, with_steps};
use fieldcraft_model::{
    BeamProfile, BoundaryCondition, Coupling, CouplingKind, Distribution, Dofs, Elastic, ElementFace, ElementKind, Formulation, Integration, Load,
    Material, Model, Procedure, ReferencePoint, Section, SectionKind, Set, Step, Target,
};

// ---------------------------------------------------------------------------------------------
// A minimal reader of our own output.

#[derive(Debug, Default)]
struct Read {
    sol: u32,
    /// Per subcase: its case control lines as `KEY = value`.
    subcases: Vec<BTreeMap<String, String>>,
    /// Bulk entries: name and data fields (trimmed, blanks as "").
    cards: Vec<(String, Vec<String>)>,
}

fn read(deck: &str) -> Read {
    let mut r = Read::default();
    let mut lines = deck.lines();
    for line in lines.by_ref() {
        let l = line.trim();
        if l.starts_with('$') {
            continue;
        }
        if let Some(s) = l.strip_prefix("SOL ") {
            r.sol = s.parse().unwrap();
        } else if l.starts_with("SUBCASE") {
            r.subcases.push(BTreeMap::new());
        } else if l == "BEGIN BULK" {
            break;
        } else if let (Some((k, v)), Some(sub)) = (l.split_once('='), r.subcases.last_mut()) {
            sub.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    for line in lines {
        assert!(line.len() <= 80, "{line}");
        assert!(!line.contains('\t'));
        if line.starts_with('$') || line == "ENDDATA" {
            continue;
        }
        let (first, rest) = line.split_at(8.min(line.len()));
        let large = first.trim_end().ends_with('*');
        let width = if large { 16 } else { 8 };
        let per_line = if large { 4 } else { 8 };
        // Every physical line holds its full count of fields; trailing blanks are trimmed away.
        let mut fields: Vec<String> = rest.as_bytes().chunks(width).map(|c| String::from_utf8_lossy(c).trim().to_string()).collect();
        fields.resize(per_line, String::new());
        let continuation = first.trim().is_empty() || first.trim() == "*";
        if continuation {
            r.cards.last_mut().expect("a continuation follows an entry").1.extend(fields);
        } else {
            r.cards.push((first.trim().trim_end_matches('*').to_string(), fields));
        }
    }
    r
}

/// A Nastran real: an exponent may come without `E`.
fn real(s: &str) -> f64 {
    let split = s.char_indices().skip(1).find(|(i, c)| (*c == '+' || *c == '-') && !s[..*i].ends_with(['E', 'e']));
    match split {
        Some((i, _)) => format!("{}E{}", &s[..i], &s[i..]).parse().unwrap(),
        None => s.parse().unwrap(),
    }
}

fn int(s: &str) -> u32 {
    s.parse().unwrap_or_else(|_| panic!("`{s}` is not an integer"))
}

impl Read {
    fn cards<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Vec<String>> + 'a {
        self.cards.iter().filter(move |(n, _)| n == name).map(|(_, f)| f)
    }

    fn grids(&self) -> BTreeMap<u32, [f64; 3]> {
        self.cards("GRID").map(|f| (int(&f[0]), [real(&f[2]), real(&f[3]), real(&f[4])])).collect()
    }

    /// Element connectivity in Nastran order, by label (solids and shells).
    fn elements(&self) -> BTreeMap<u32, (String, Vec<u32>)> {
        let mut out = BTreeMap::new();
        for (name, f) in &self.cards {
            let n = match name.as_str() {
                "CTETRA" | "CPENTA" | "CHEXA" => f.len() - 2,
                "CTRIA3" => 3,
                "CQUAD4" => 4,
                "CTRIA6" => 6,
                "CQUAD8" => 8,
                "CBEAM" | "CROD" => 2,
                _ => continue,
            };
            let nodes = f[2..2 + n].iter().filter(|s| !s.is_empty()).map(|s| int(s)).collect();
            out.insert(int(&f[0]), (name.clone(), nodes));
        }
        out
    }

    /// Constraint set `sid`: node → (degree of freedom → value).
    fn constraints(&self, sid: u32) -> BTreeMap<u32, BTreeMap<u8, f64>> {
        let mut out: BTreeMap<u32, BTreeMap<u8, f64>> = BTreeMap::new();
        let dofs = |c: &str| c.bytes().map(|b| b - b'0').collect::<Vec<u8>>();
        for f in self.cards("SPC1").filter(|f| int(&f[0]) == sid) {
            for g in f[2..].iter().filter(|s| !s.is_empty()) {
                for d in dofs(&f[1]) {
                    assert!(out.entry(int(g)).or_default().insert(d, 0.0).is_none(), "constrained twice");
                }
            }
        }
        for f in self.cards("SPC").filter(|f| int(&f[0]) == sid) {
            for d in dofs(&f[2]) {
                assert!(out.entry(int(&f[1])).or_default().insert(d, real(&f[3])).is_none(), "constrained twice");
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------------------------
// Our own models.

fn steel() -> Material {
    Material { elastic: Some(Elastic::Isotropic { young: 210_000.0, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() }
}

/// A bar of `n` unit hex8 along X; node (i, y, z) = 1 + 4i + 2y + z.
fn bar(n: u32) -> Model {
    let mut m = Model::default();
    for i in 0..=n {
        for (j, (y, z)) in [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)].into_iter().enumerate() {
            m.set_node(1 + 4 * i + j as u32, [f64::from(i), y, z]).unwrap();
        }
    }
    let g = |i: u32, y: u32, z: u32| 1 + 4 * i + 2 * y + z;
    for e in 0..n {
        let nodes = vec![g(e, 0, 0), g(e + 1, 0, 0), g(e + 1, 1, 0), g(e, 1, 0), g(e, 0, 1), g(e + 1, 0, 1), g(e + 1, 1, 1), g(e, 1, 1)];
        m.set_element(e + 1, ElementKind::Hex8, nodes).unwrap();
    }
    m.add_set("all", Set { elements: (1..=n).collect(), ..Set::default() }).unwrap();
    m.add_set("root", Set { faces: vec![ElementFace { element: 1, face: 6 }], ..Set::default() }).unwrap();
    m.add_set("end", Set { faces: vec![ElementFace { element: n, face: 4 }], ..Set::default() }).unwrap();
    m.add_material("Steel", steel()).unwrap();
    m.add_section("bar", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    m
}

fn static_step(name: &str, bcs: &[&str], loads: &[&str]) -> Step {
    let mut s = Step::new(name, Procedure::Static {});
    s.boundary_conditions = bcs.iter().map(|s| s.to_string()).collect();
    s.loads = loads.iter().map(|s| s.to_string()).collect();
    s
}

fn cantilever() -> Model {
    let mut m = bar(3);
    m.add_load("tip", Load::Force { target: Target::Set("end".into()), vector: [0.0, 0.0, -100.0], distribution: Distribution::Total }).unwrap();
    m.add_step(static_step("Bend", &["clamp"], &["tip"]), None).unwrap();
    m
}

/// A 2 × 2 quad4 plate in z = 0, 10 × 10, clamped along x = 0, under pressure and gravity.
fn plate() -> Model {
    let mut m = Model::default();
    for j in 0..3u32 {
        for i in 0..3u32 {
            m.set_node(1 + i + 3 * j, [5.0 * f64::from(i), 5.0 * f64::from(j), 0.0]).unwrap();
        }
    }
    for (e, (i, j)) in [(0, 0), (1, 0), (0, 1), (1, 1)].into_iter().enumerate() {
        let n = 1 + i + 3 * j;
        m.set_element(e as u32 + 1, ElementKind::Quad4, vec![n, n + 1, n + 4, n + 3]).unwrap();
    }
    m.add_set("all", Set { elements: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    m.add_set("edge", Set { nodes: vec![1, 4, 7], ..Set::default() }).unwrap();
    m.add_set("top", Set { faces: (1..=4).map(|element| ElementFace { element, face: 1 }).collect(), ..Set::default() }).unwrap();
    m.add_material(
        "Aluminium",
        Material { elastic: Some(Elastic::Isotropic { young: 70_000.0, poisson: 0.33 }), density: Some(2.7e-9), ..Material::default() },
    )
    .unwrap();
    let shell = SectionKind::Shell { thickness: 2.0, offset: 0.0, points: 5 };
    m.add_section("skin", Section { set: "all".into(), material: "Aluminium".into(), kind: shell, formulation: Formulation::default() }).unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("edge".into()) }).unwrap();
    m.add_load("press", Load::Pressure { set: "top".into(), magnitude: 0.01 }).unwrap();
    m.add_load("weight", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: None }).unwrap();
    m.add_step(static_step("Sag", &["clamp"], &["press", "weight"]), None).unwrap();
    m
}

/// Reference corners of each solid kind, slightly distorted, with straight edges (mid-side nodes
/// at edge midpoints, by the model's own edge table).
fn solid(kind: ElementKind) -> Vec<[f64; 3]> {
    let corners: Vec<[f64; 3]> = match kind {
        ElementKind::Tet4 | ElementKind::Tet10 => vec![[0.0, 0.0, 0.0], [2.0, 0.1, 0.0], [0.2, 1.5, 0.1], [0.3, 0.2, 1.7]],
        ElementKind::Wedge6 => vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.1], [0.1, 1.5, 0.0], [0.1, 0.1, 1.2], [2.1, 0.0, 1.3], [0.0, 1.6, 1.1]],
        _ => vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [2.2, 1.4, 0.1],
            [-0.1, 1.5, 0.0],
            [0.1, 0.0, 1.0],
            [2.0, 0.2, 1.1],
            [2.1, 1.6, 1.2],
            [0.0, 1.4, 0.9],
        ],
    };
    let mut xyz = corners.clone();
    xyz.resize(kind.node_count(), [0.0; 3]);
    for e in kind.edges() {
        if let [a, b, mid] = **e {
            xyz[mid] = [0, 1, 2].map(|i| (corners[a][i] + corners[b][i]) / 2.0);
        }
    }
    xyz
}

/// One element of `kind` with its nodes labelled in reverse (so labels and positions differ), a
/// section, and a static step.
fn one_element(kind: ElementKind) -> Model {
    let mut m = Model::default();
    let xyz = solid(kind);
    let n = xyz.len() as u32;
    for (i, p) in xyz.iter().enumerate() {
        m.set_node(100 - i as u32, *p).unwrap();
    }
    m.set_element(7, kind, (0..n).map(|i| 100 - i).collect()).unwrap();
    m.add_set("all", Set { elements: vec![7], ..Set::default() }).unwrap();
    m.add_set("base", Set { faces: vec![ElementFace { element: 7, face: 1 }], ..Set::default() }).unwrap();
    m.add_material("Steel", steel()).unwrap();
    m.add_section("s", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    m.add_load("push", Load::Pressure { set: "base".into(), magnitude: 1.0 }).unwrap();
    m.add_step(static_step("Load", &["hold"], &["push"]), None).unwrap();
    m
}

/// A hex20 block whose top face is coupled to a reference point, clamped at the base, in a modes
/// step.
fn modes() -> Model {
    let mut m = one_element(ElementKind::Hex20);
    m.add_set("top", Set { faces: vec![ElementFace { element: 7, face: 2 }], ..Set::default() }).unwrap();
    m.add_reference_point("RP", ReferencePoint { xyz: [1.0, 0.75, 2.0], label: 0 }).unwrap();
    m.add_coupling("lid", Coupling { reference_point: "RP".into(), set: "top".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
    m.remove_step("Load").unwrap();
    let mut s = Step::new("Modes", Procedure::Frequency { modes: 6, min_hz: None, max_hz: Some(5000.0) });
    s.boundary_conditions.insert("hold".into());
    m.add_step(s, None).unwrap();
    m
}

fn golden(name: &str, deck: &str) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/nastran").join(format!("{name}.bdf"));
    if std::env::var_os("FIELDCRAFT_BLESS").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, deck).unwrap();
    }
    let want = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e} (bless with FIELDCRAFT_BLESS=1)", path.display()));
    assert_eq!(deck, want, "{name}: the deck changed (bless with FIELDCRAFT_BLESS=1 after checking it)");
}

fn unsupported(m: &Model) -> String {
    match nastran::write(m) {
        Err(Error::Unsupported { format: "Nastran", what }) => what,
        other => panic!("expected unsupported, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------------------------

#[test]
fn golden_decks() {
    for (name, model) in [("cantilever", cantilever()), ("plate", plate()), ("modes", modes())] {
        let deck = nastran::write(&model).unwrap();
        assert_eq!(nastran::write(&model).unwrap(), deck, "deterministic");
        golden(name, &deck);
    }
}

#[test]
fn cantilever_reads_back() {
    let m = cantilever();
    let r = read(&nastran::write(&m).unwrap());
    assert_eq!(r.sol, 101);
    assert_eq!(r.subcases.len(), 1);
    let case = &r.subcases[0];
    assert_eq!(case["SUBTITLE"], "Bend");
    for request in ["DISPLACEMENT", "SPCFORCES", "STRESS", "STRAIN"] {
        assert_eq!(case[request], "ALL");
    }
    // Grids are exact.
    let grids = r.grids();
    assert_eq!(grids.len(), m.nodes().len());
    for (label, node) in m.nodes() {
        assert_eq!(grids[label], node.xyz);
    }
    // Elements in Nastran order map back to the model's.
    let elements = r.elements();
    for (label, el) in m.elements() {
        let (card, nodes) = &elements[label];
        assert_eq!(card, "CHEXA");
        assert_eq!(*nodes, el.nodes);
    }
    // The clamp holds every degree of freedom of the root face's nodes.
    let spc = r.constraints(int(&case["SPC"]));
    assert_eq!(spc.keys().copied().collect::<Vec<_>>(), vec![1, 2, 3, 4]);
    assert!(spc.values().all(|d| d.len() == 6 && d.values().all(|v| *v == 0.0)));
    // The tip force totals -100 in Z, a quarter on each end node.
    let sid = int(&case["LOAD"]);
    let mut total = [0.0; 3];
    for f in r.cards("FORCE").filter(|f| int(&f[0]) == sid) {
        assert!((13..=16).contains(&int(&f[1])));
        for k in 0..3 {
            total[k] += real(&f[3]) * real(&f[4 + k]);
        }
    }
    assert_eq!(total, [0.0, 0.0, -100.0]);
    // Two-point full integration by default.
    let psolid = r.cards("PSOLID").next().unwrap();
    assert_eq!((psolid[3].as_str(), psolid[5].as_str()), ("TWO", "FULL"));
}

#[test]
fn plate_reads_back() {
    let r = read(&nastran::write(&plate()).unwrap());
    let case = &r.subcases[0];
    // Pressure and gravity combine through a LOAD entry.
    let load = r.cards("LOAD").next().unwrap();
    assert_eq!(int(&load[0]), int(&case["LOAD"]));
    let parts: Vec<u32> = load[2..].chunks(2).filter(|c| !c[1].is_empty()).map(|c| int(&c[1])).collect();
    assert_eq!(parts.len(), 2);
    let grav = r.cards("GRAV").next().unwrap();
    assert!(parts.contains(&int(&grav[0])));
    assert_eq!([real(&grav[2]) * real(&grav[3]), real(&grav[4]), real(&grav[5])], [0.0, 0.0, -9810.0]);
    // Pressure on the side the normal points to pushes against the normal.
    let pload: Vec<_> = r.cards("PLOAD4").collect();
    assert_eq!(pload.len(), 4);
    assert!(pload.iter().all(|f| real(&f[2]) == -0.01 && parts.contains(&int(&f[0]))));
    let pshell = r.cards("PSHELL").next().unwrap();
    assert_eq!(real(&pshell[2]), 2.0);
    assert_eq!(r.elements().values().filter(|(c, _)| c == "CQUAD4").count(), 4);
    assert_eq!(r.constraints(int(&case["SPC"])).keys().copied().collect::<Vec<_>>(), vec![1, 4, 7]);
}

#[test]
fn modes_reads_back() {
    let m = modes();
    let r = read(&nastran::write(&m).unwrap());
    assert_eq!(r.sol, 103);
    let case = &r.subcases[0];
    assert!(!case.contains_key("LOAD"));
    assert_eq!(case["DISPLACEMENT"], "ALL");
    let eigrl = r.cards("EIGRL").next().unwrap();
    assert_eq!(int(&eigrl[0]), int(&case["METHOD"]));
    assert_eq!((eigrl[1].as_str(), real(&eigrl[2]), int(&eigrl[3])), ("", 5000.0, 6));
    // The reference point is a grid; the coupling a rigid element on the top face's eight nodes.
    let rp = m.reference_point("RP").unwrap().label;
    assert!(r.grids().contains_key(&rp));
    let rbe2 = r.cards("RBE2").next().unwrap();
    assert_eq!((int(&rbe2[0]), int(&rbe2[1]), int(&rbe2[2])), (8, rp, 123456));
    let mut dependent: Vec<u32> = rbe2[3..].iter().filter(|s| !s.is_empty()).map(|s| int(s)).collect();
    dependent.sort_unstable();
    assert_eq!(dependent, m.resolve_set("top").unwrap().nodes);
    // The CHEXA's mid-side blocks are swapped and map back.
    let (_, nodes) = &r.elements()[&7];
    let el = m.element(7).unwrap();
    let back: Vec<u32> = (0..20).map(|i| nodes[nastran::order(ElementKind::Hex20).iter().position(|&p| p == i).unwrap()]).collect();
    assert_eq!(back, el.nodes);
    assert_ne!(*nodes, el.nodes);
    assert!(r.cards("PARAM").any(|f| f[0] == "COUPMASS"));
}

// ---------------------------------------------------------------------------------------------
// Node order against Nastran's own conventions.

/// Nastran's natural coordinates of each grid, from its corner layout and its edge table (which
/// grids sit between which corners), written independently of the writer's order table.
fn nastran_natural(card: &str, n: usize) -> Vec<[f64; 3]> {
    let (corners, edges): (Vec<[f64; 3]>, Vec<(usize, usize)>) = match card {
        "CHEXA" => (
            vec![
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ],
            // G9–G12 bottom, G13–G16 vertical, G17–G20 top.
            vec![(1, 2), (2, 3), (3, 4), (4, 1), (1, 5), (2, 6), (3, 7), (4, 8), (5, 6), (6, 7), (7, 8), (8, 5)],
        ),
        // (ξ, η) in the triangle, ζ through the thickness; G7–G9 one triangle, G10–G12 vertical.
        "CPENTA" => (
            vec![[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 1.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
            vec![(1, 2), (2, 3), (3, 1), (1, 4), (2, 5), (3, 6), (4, 5), (5, 6), (6, 4)],
        ),
        // G5–G10 on edges 1-2, 2-3, 3-1, 1-4, 2-4, 3-4.
        _ => (vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], vec![(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)]),
    };
    let mut out = corners.clone();
    for (a, b) in edges.into_iter().take(n - corners.len()) {
        out.push([0, 1, 2].map(|i| (corners[a - 1][i] + corners[b - 1][i]) / 2.0));
    }
    out
}

/// Linear (corner) shape functions of each Nastran solid.
fn corner_shape(card: &str, p: [f64; 3]) -> Vec<f64> {
    let [x, y, z] = p;
    match card {
        "CHEXA" => nastran_natural("CHEXA", 8).iter().map(|c| (1.0 + x * c[0]) * (1.0 + y * c[1]) * (1.0 + z * c[2]) / 8.0).collect(),
        "CPENTA" => {
            let l = [1.0 - x - y, x, y];
            [(1.0 - z) / 2.0, (1.0 + z) / 2.0].iter().flat_map(|h| l.map(|li| li * h)).collect()
        }
        _ => vec![1.0 - x - y - z, x, y, z],
    }
}

/// The Jacobian determinant of the corner map at natural point `p`, by central differences.
fn jacobian(card: &str, corners: &[[f64; 3]], p: [f64; 3]) -> f64 {
    let h = 1e-6;
    let at = |q: [f64; 3]| -> [f64; 3] {
        let n = corner_shape(card, q);
        [0, 1, 2].map(|k| n.iter().zip(corners).map(|(w, c)| w * c[k]).sum())
    };
    let d: Vec<[f64; 3]> = (0..3)
        .map(|axis| {
            let (mut a, mut b) = (p, p);
            a[axis] += h;
            b[axis] -= h;
            let (fa, fb) = (at(a), at(b));
            [0, 1, 2].map(|k| (fa[k] - fb[k]) / (2.0 * h))
        })
        .collect();
    d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) - d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0])
        + d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0])
}

/// Check a written solid against Nastran's conventions: positive Jacobian at points near every
/// corner and at the centre, and every mid-side grid at the midpoint of its Nastran edge.
fn check_nastran_solid(card: &str, xyz: &[[f64; 3]]) -> Result<(), String> {
    let natural = nastran_natural(card, xyz.len());
    let ncorner = corner_shape(card, [0.0; 3]).len();
    let corners = &xyz[..ncorner];
    let centre: [f64; 3] = [0, 1, 2].map(|k| natural[..ncorner].iter().map(|c| c[k]).sum::<f64>() / ncorner as f64);
    for c in natural[..ncorner].iter().chain(std::iter::once(&centre)) {
        let p = [0, 1, 2].map(|k| centre[k] + 0.9 * (c[k] - centre[k]));
        let j = jacobian(card, corners, p);
        if !(j > 0.0) {
            return Err(format!("{card}: Jacobian {j} at {p:?}"));
        }
    }
    for (i, nat) in natural.iter().enumerate().skip(ncorner) {
        let n = corner_shape(card, *nat);
        let want: [f64; 3] = [0, 1, 2].map(|k| n.iter().zip(corners).map(|(w, c)| w * c[k]).sum());
        if (0..3).any(|k| (want[k] - xyz[i][k]).abs() > 1e-12) {
            return Err(format!("{card}: grid G{} is not on its edge", i + 1));
        }
    }
    Ok(())
}

#[test]
fn solids_keep_nastran_node_order_and_positive_jacobians() {
    for kind in [ElementKind::Tet4, ElementKind::Tet10, ElementKind::Wedge6, ElementKind::Hex8, ElementKind::Hex20] {
        let m = one_element(kind);
        let r = read(&nastran::write(&m).unwrap());
        let grids = r.grids();
        let (card, nodes) = &r.elements()[&7];
        assert_eq!(card, nastran::element_card(kind));
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|n| grids[n]).collect();
        check_nastran_solid(card, &xyz).unwrap_or_else(|e| panic!("{kind:?}: {e}"));
        // Mapped back, the connectivity is the model's, and so is its positive Jacobian.
        let order = nastran::order(kind);
        let back: Vec<u32> = (0..nodes.len()).map(|i| nodes[order.iter().position(|&p| p == i).unwrap()]).collect();
        assert_eq!(back, m.element(7).unwrap().nodes, "{kind:?}");

        // The face pressure names a corner on the face, and the corner diagonally opposite (hex),
        // the corner off the face (tet), or nothing (a wedge's triangle).
        let p = r.cards("PLOAD4").next().unwrap();
        assert_eq!(real(&p[2]), 1.0);
        let face: Vec<u32> = m.face_nodes(ElementFace { element: 7, face: 1 }).unwrap();
        assert!(face.contains(&int(&p[6])));
        match kind {
            ElementKind::Tet4 | ElementKind::Tet10 => assert!(!face.contains(&int(&p[7])) && nodes[..4].contains(&int(&p[7]))),
            ElementKind::Wedge6 => assert_eq!(p.get(7).map(String::as_str).unwrap_or(""), ""),
            _ => assert_eq!(int(&p[7]), face[2]),
        }
    }
    // Without the reordering, a hex20 breaks Nastran's convention: the check catches it.
    let m = one_element(ElementKind::Hex20);
    let model_order: Vec<[f64; 3]> = m.element(7).unwrap().nodes.iter().map(|n| m.node(*n).unwrap().xyz).collect();
    assert!(check_nastran_solid("CHEXA", &model_order).unwrap_err().contains("G13"));
}

// ---------------------------------------------------------------------------------------------
// The rest of the model.

#[test]
fn beams_rods_orthotropy_couplings_and_prescribed_motion() {
    let mut m = cantilever();
    // An orthotropic solid material.
    let ortho = Material {
        elastic: Some(Elastic::EngineeringConstants { young: [140e3, 10e3, 12e3], poisson: [0.3, 0.28, 0.4], shear: [5e3, 4.5e3, 3.8e3] }),
        density: Some(1.6e-9),
        ..Material::default()
    };
    m.add_material("Carbon", ortho.clone()).unwrap();
    m.replace_section(
        "bar",
        Section { set: "all".into(), material: "Carbon".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() },
    )
    .unwrap();
    // A beam and a rod hanging off the tip, and a shell lying in z = 0 with the same material.
    m.set_node(50, [5.0, 0.0, 0.0]).unwrap();
    m.set_node(51, [4.0, 1.0, 0.0]).unwrap();
    m.set_element(10, ElementKind::Beam2, vec![13, 50]).unwrap();
    m.set_element(11, ElementKind::Rod2, vec![15, 51]).unwrap();
    m.set_element(12, ElementKind::Tri3, vec![13, 50, 51]).unwrap();
    m.add_set("beam", Set { elements: vec![10], ..Set::default() }).unwrap();
    m.add_set("rod", Set { elements: vec![11], ..Set::default() }).unwrap();
    m.add_set("web", Set { elements: vec![12], ..Set::default() }).unwrap();
    let i = BeamProfile::I { height: 10.0, width: 6.0, web: 1.0, flange: 1.5 };
    m.add_section(
        "beam",
        Section {
            set: "beam".into(),
            material: "Steel".into(),
            kind: SectionKind::Beam { profile: i, orientation: [0.0, 0.0, 1.0] },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    m.add_section(
        "rod",
        Section { set: "rod".into(), material: "Steel".into(), kind: SectionKind::Rod { area: 3.5 }, formulation: Formulation::default() },
    )
    .unwrap();
    m.add_section(
        "web",
        Section {
            set: "web".into(),
            material: "Carbon".into(),
            kind: SectionKind::Shell { thickness: 1.0, offset: 0.5, points: 3 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    // A distributing coupling from a reference point to the end face, loaded with a moment.
    m.add_reference_point("Hub", ReferencePoint { xyz: [3.0, 0.5, 0.5], label: 0 }).unwrap();
    m.add_coupling("spider", Coupling { reference_point: "Hub".into(), set: "end".into(), kind: CouplingKind::Distributing, dofs: Dofs::all() })
        .unwrap();
    m.add_load("twist", Load::Moment { target: Target::ReferencePoint("Hub".into()), vector: [1000.0, 0.0, 0.0], distribution: Distribution::Total })
        .unwrap();
    // Prescribed motion of the root in a second step.
    m.add_boundary_condition(
        "lift",
        BoundaryCondition::Displacement { target: Target::Set("root".into()), values: [Some(0.0), Some(0.0), Some(0.25), None, None, None] },
    )
    .unwrap();
    m.add_step(static_step("Lift", &["lift"], &["twist"]), None).unwrap();

    let deck = nastran::write(&m).unwrap();
    let r = read(&deck);
    assert_eq!(r.subcases.len(), 2);

    // MAT9 carries the stiffness of the engineering constants; MAT8 the shell's constants.
    let mat9 = r.cards("MAT9").next().unwrap();
    let c = nastran::orthotropic_stiffness([140e3, 10e3, 12e3], [0.3, 0.28, 0.4], [5e3, 4.5e3, 3.8e3]).unwrap();
    let close = |s: &str, v: f64| (real(s) - v).abs() <= 1e-13 * v.abs();
    assert!(close(&mat9[1], c[0][0]) && close(&mat9[7], c[1][0]) && close(&mat9[12], c[1][2]));
    assert_eq!((real(&mat9[16]), real(&mat9[19]), real(&mat9[21])), (5e3, 3.8e3, 4.5e3));
    assert_eq!(real(&mat9[22]), 1.6e-9);
    let mat8 = r.cards("MAT8").next().unwrap();
    assert_eq!(mat8[1..7].iter().map(|s| real(s)).collect::<Vec<_>>(), vec![140e3, 10e3, 0.3, 5e3, 4.5e3, 3.8e3]);
    // The shell's material axes are the basic ones, and its offset puts the mid-surface below the nodes.
    let tri = r.cards("CTRIA3").next().unwrap();
    assert_eq!((int(&tri[5]), real(&tri[6])), (1, -0.5));
    assert!(r.cards("CORD2R").next().is_some());

    // The beam: its I section's flanges lie along the orientation (Nastran's element z axis).
    let pbeaml = r.cards("PBEAML").next().unwrap();
    assert_eq!(pbeaml[3], "I");
    assert_eq!(pbeaml[8..14].iter().map(|s| real(s)).collect::<Vec<_>>(), vec![10.0, 6.0, 6.0, 1.0, 1.5, 1.5]);
    let cbeam = r.cards("CBEAM").next().unwrap();
    let v = [real(&cbeam[4]), real(&cbeam[5]), real(&cbeam[6])];
    let axis = [2.0, 0.0, 0.0]; // node 13 (3, 0, 0) to node 50 (5, 0, 0)
    let z = [axis[1] * v[2] - axis[2] * v[1], axis[2] * v[0] - axis[0] * v[2], axis[0] * v[1] - axis[1] * v[0]];
    assert!(z[0] == 0.0 && z[1] == 0.0 && z[2] != 0.0, "element z along the orientation: {z:?}");
    let prod = r.cards("PROD").next().unwrap();
    assert_eq!(real(&prod[2]), 3.5);

    // The distributing coupling: the hub depends on the end face's translations, equally weighted.
    let rbe3 = r.cards("RBE3").next().unwrap();
    let hub = m.reference_point("Hub").unwrap().label;
    assert_eq!((rbe3[1].as_str(), int(&rbe3[2]), int(&rbe3[3]), real(&rbe3[4]), int(&rbe3[5])), ("", hub, 123456, 1.0, 123));
    let moment = r.cards("MOMENT").next().unwrap();
    assert_eq!((int(&moment[1]), real(&moment[4])), (hub, 1000.0));

    // Step 2 prescribes Z = 0.25 at the root and holds X and Y.
    let spc = r.constraints(int(&r.subcases[1]["SPC"]));
    assert_eq!(spc.len(), 4);
    for dofs in spc.values() {
        assert_eq!(dofs.iter().map(|(d, v)| (*d, *v)).collect::<Vec<_>>(), vec![(1, 0.0), (2, 0.0), (3, 0.25)]);
    }
    assert_ne!(r.subcases[0]["SPC"], r.subcases[1]["SPC"]);
}

#[test]
fn steps_are_chosen_and_mixed_solutions_are_refused() {
    let mut m = cantilever();
    let mut s = Step::new("Modes", Procedure::Frequency { modes: 4, min_hz: Some(1.0), max_hz: None });
    s.boundary_conditions.insert("clamp".into());
    m.add_step(s, None).unwrap();
    assert!(unsupported(&m).contains("separate decks"));

    let only = |names: &[&str]| names.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    let r = read(&nastran::write(&with_steps(&m, &only(&["Modes"])).unwrap()).unwrap());
    assert_eq!((r.sol, r.subcases.len()), (103, 1));
    let eigrl = r.cards("EIGRL").next().unwrap();
    assert_eq!((real(&eigrl[1]), eigrl[2].as_str(), int(&eigrl[3])), (1.0, "", 4));
    // A modes step holds the clamp but carries no loads.
    assert!(r.cards("FORCE").next().is_none());
    assert_eq!(read(&nastran::write(&with_steps(&m, &only(&["Bend"])).unwrap()).unwrap()).sol, 101);
    let nastran = fieldcraft_decks::format("nastran").unwrap();
    assert_eq!(
        nastran.write_steps(&common::si(&m), Some(&only(&["Bend"])), fieldcraft_model::units::Preset::MmTonneS).unwrap().text,
        nastran::write(&with_steps(&m, &only(&["Bend"])).unwrap()).unwrap()
    );

    for bad in [only(&[]), only(&["Nope"]), only(&[&"x".repeat(100_000)]), only(&["", "\u{0}"])] {
        assert!(matches!(with_steps(&m, &bad), Err(Error::Model(_))), "{bad:?}");
    }
}

#[test]
fn what_nastran_cannot_express_is_named() {
    let solid = |f: Formulation| Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: f };
    let mut m = cantilever();
    m.replace_section("bar", solid(Formulation { hybrid: true, ..Formulation::default() })).unwrap();
    assert!(unsupported(&m).contains("hybrid"));
    m.replace_section("bar", solid(Formulation { integration: Integration::Reduced, ..Formulation::default() })).unwrap();
    assert!(unsupported(&m).contains("reduced integration of hex8"));
    m.replace_section("bar", solid(Formulation { incompatible_modes: true, ..Formulation::default() })).unwrap();
    let r = read(&nastran::write(&m).unwrap());
    assert_eq!(r.cards("PSOLID").next().unwrap()[3], "BUBBLE");

    // Several problems are reported together.
    let mut m = cantilever();
    m.add_load("heavy", Load::BodyForce { set: "all".into(), vector: [0.0, 0.0, -1.0] }).unwrap();
    m.add_set("first", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_load("partial", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: Some("first".into()) }).unwrap();
    m.add_reference_point("RP", ReferencePoint { xyz: [3.0, 0.5, 0.5], label: 0 }).unwrap();
    m.add_coupling("grip", Coupling { reference_point: "RP".into(), set: "root".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
    m.replace_step(static_step("Bend", &["clamp"], &["tip", "heavy", "partial"])).unwrap();
    let what = unsupported(&m);
    for needle in ["body force", "gravity on part", "coupling `grip`"] {
        assert!(what.contains(needle), "{needle}: {what}");
    }

    // Orthotropic beams, tilted orthotropic shells, one-point shells.
    let mut m = plate();
    let ortho = Elastic::EngineeringConstants { young: [140e3, 10e3, 10e3], poisson: [0.3, 0.3, 0.4], shear: [5e3, 5e3, 3.8e3] };
    m.replace_material("Aluminium", Material { elastic: Some(ortho), density: Some(1.6e-9), ..Material::default() }).unwrap();
    m.set_node(9, [10.0, 10.0, 0.5]).unwrap();
    assert!(unsupported(&m).contains("planes of constant Z"));
    m.set_node(9, [10.0, 10.0, 0.0]).unwrap();
    nastran::write(&m).unwrap();
    m.replace_section(
        "skin",
        Section {
            set: "all".into(),
            material: "Aluminium".into(),
            kind: SectionKind::Shell { thickness: 2.0, offset: 0.0, points: 1 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    assert!(unsupported(&m).contains("one integration point"));
    let mut m = plate();
    m.set_element(9, ElementKind::Rod2, vec![1, 9]).unwrap();
    m.add_set("bars", Set { elements: vec![9], ..Set::default() }).unwrap();
    m.add_material("Fibre", Material { elastic: Some(ortho), density: Some(1.6e-9), ..Material::default() }).unwrap();
    m.add_section(
        "bars",
        Section { set: "bars".into(), material: "Fibre".into(), kind: SectionKind::Rod { area: 1.0 }, formulation: Formulation::default() },
    )
    .unwrap();
    assert!(unsupported(&m).contains("beams or rods"));

    // A model that does not resolve is the model's error.
    assert!(matches!(nastran::write(&Model::default()), Err(Error::Model(_))));
}

#[test]
fn conflicting_values_and_ids_at_the_top_of_the_range() {
    let mut m = cantilever();
    m.add_boundary_condition(
        "raise",
        BoundaryCondition::Displacement { target: Target::Set("root".into()), values: [None, None, Some(1.0), None, None, None] },
    )
    .unwrap();
    m.replace_step(static_step("Bend", &["clamp", "raise"], &["tip"])).unwrap();
    assert!(unsupported(&m).contains("different values"));

    // An element at the top of the label range leaves no id for a coupling's rigid element.
    let mut m = cantilever();
    m.set_node(99_999_998, [9.0, 9.0, 9.0]).unwrap();
    m.set_element(99_999_999, ElementKind::Rod2, vec![1, 99_999_998]).unwrap();
    m.replace_set("all", Set { elements: vec![1, 2, 3], ..Set::default() }).unwrap();
    m.add_set("bar2", Set { elements: vec![99_999_999], ..Set::default() }).unwrap();
    m.add_section(
        "bar2",
        Section { set: "bar2".into(), material: "Steel".into(), kind: SectionKind::Rod { area: 1.0 }, formulation: Formulation::default() },
    )
    .unwrap();
    let deck = nastran::write(&m).unwrap();
    assert!(deck.contains("CROD    99999999"));
    m.add_reference_point("RP", ReferencePoint { xyz: [3.0, 0.5, 0.5], label: 1000 }).unwrap();
    m.add_coupling("pad", Coupling { reference_point: "RP".into(), set: "end".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
    assert!(matches!(nastran::write(&m), Err(Error::Limit(_))));
}

/// Plasticity: `MATS1` beside `MAT1` (bilinear by slope, a table through `TABLES1`); Johnson–Cook,
/// damage and erosion are reported.
#[test]
fn plasticity_is_written_as_mats1() {
    use fieldcraft_model::{Erosion, Hardening, HardeningPoint, Plastic};
    let with = |p: Plastic| {
        let mut m = cantilever();
        m.replace_material("Steel", Material { plastic: Some(p), ..steel() }).unwrap();
        m
    };
    let deck = nastran::write(&with(Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 1000.0, hardening: Hardening::Kinematic })).unwrap();
    let r = read(&deck);
    let (_, f) = r.cards.iter().find(|(n, _)| n == "MATS1").unwrap();
    assert_eq!((f[0].as_str(), f[1].as_str(), f[2].as_str(), int(&f[4]), int(&f[5])), ("1", "", "PLASTIC", 1, 2), "{f:?}");
    assert_eq!((real(&f[3]), real(&f[6])), (1000.0, 250.0), "{f:?}");
    let table = Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: 210.0 }, HardeningPoint { strain: 0.1, stress: 315.0 }] };
    let deck = nastran::write(&with(table)).unwrap();
    let r = read(&deck);
    let (_, mats1) = r.cards.iter().find(|(n, _)| n == "MATS1").unwrap();
    let (_, table) = r.cards.iter().find(|(n, _)| n == "TABLES1").unwrap();
    assert_eq!(mats1[1], table[0], "MATS1 names its table");
    // The curve after the header line: the origin, then total strain σ/E + ε_p against σ.
    let data: Vec<&String> = table[8..].iter().filter(|f| !f.is_empty()).collect();
    assert_eq!(data.last().map(|s| s.as_str()), Some("ENDT"));
    let values: Vec<f64> = data[..data.len() - 1].iter().map(|s| real(s)).collect();
    let want = [0.0, 0.0, 0.001, 210.0, 0.1015, 315.0];
    assert!(values.len() == 6 && values.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-12), "{values:?}");
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
    let mut m = with(jc);
    let mut steel = m.material("Steel").unwrap().clone();
    steel.erosion = Some(Erosion { plastic_strain: 1.0 });
    m.replace_material("Steel", steel).unwrap();
    let Err(Error::Unsupported { what, .. }) = nastran::write(&m) else { panic!() };
    assert!(what.contains("Johnson–Cook plasticity") && what.contains("damage and erosion"), "{what}");
}

#[test]
fn nonlinear_static_steps_are_refused_by_name() {
    let mut m = cantilever();
    let name = m.steps()[0].name.clone();
    let mut s = m.steps()[0].clone();
    s.set_procedure(Procedure::NonlinearStatic(fieldcraft_model::NonlinearStatic::default()));
    m.replace_step(s).unwrap();
    assert!(unsupported(&m).contains(&format!("step `{name}`: nonlinear static steps (Nastran's SOL 400")), "{}", unsupported(&m));
}
