//! Abaqus writer tests: golden decks of our own small models, a re-read of the written mesh and
//! sets by a minimal test-only reader, element orientation (Jacobian sign) from the written
//! connectivity, number round trips, line limits and everything reported as unsupported.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::Error;
use fieldcraft_model::{
    BeamProfile, BoundaryCondition, Coupling, CouplingKind, Distribution, Dofs, Elastic, ElementFace, ElementKind, Formulation, Integration, Load,
    Material, Model, OutputVariable, Procedure, ReferencePoint, Section, SectionKind, Set, Step, Symmetry, Target,
};

// ---------------------------------------------------------------------------------------------
// Models

fn steel() -> Material {
    Material { elastic: Some(Elastic::Isotropic { young: 210_000.0, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() }
}

fn faces(list: &[(u32, u8)]) -> Set {
    Set { faces: list.iter().map(|&(element, face)| ElementFace { element, face }).collect(), ..Set::default() }
}

fn elements(list: impl IntoIterator<Item = u32>) -> Set {
    Set { elements: list.into_iter().collect(), ..Set::default() }
}

fn step(name: &str, procedure: Procedure, bcs: &[&str], loads: &[&str]) -> Step {
    let mut s = Step::new(name, procedure);
    s.boundary_conditions = bcs.iter().map(|s| s.to_string()).collect();
    s.loads = loads.iter().map(|s| s.to_string()).collect();
    s
}

/// A 100 × 10 × 10 cantilever of four hex8 elements along X: node (i, j, k) at
/// (25 i, 10 j, 10 k) is labelled 1 + 4 i + 2 j + k.
fn cantilever_mesh() -> Model {
    let mut m = Model::default();
    let n = |i: u32, j: u32, k: u32| 1 + 4 * i + 2 * j + k;
    for i in 0..=4 {
        for j in 0..=1 {
            for k in 0..=1 {
                m.set_node(n(i, j, k), [25.0 * f64::from(i), 10.0 * f64::from(j), 10.0 * f64::from(k)]).unwrap();
            }
        }
    }
    for i in 0..4 {
        let nodes = vec![n(i, 0, 0), n(i + 1, 0, 0), n(i + 1, 1, 0), n(i, 1, 0), n(i, 0, 1), n(i + 1, 0, 1), n(i + 1, 1, 1), n(i, 1, 1)];
        m.set_element(i + 1, ElementKind::Hex8, nodes).unwrap();
    }
    m.add_set("all", elements(1..=4)).unwrap();
    // Face 6 of element 1 lies on x = 0, face 4 of element 4 on x = 100.
    m.add_set("root", faces(&[(1, 6)])).unwrap();
    m.add_set("tip", faces(&[(4, 4)])).unwrap();
    m.add_material("Steel", steel()).unwrap();
    let formulation = Formulation { incompatible_modes: true, ..Formulation::default() };
    m.add_section("beam", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation }).unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    m
}

/// The cantilever with a total end force of 100 down.
fn cantilever() -> Model {
    let mut m = cantilever_mesh();
    m.add_load("end", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, -100.0], distribution: Distribution::Total }).unwrap();
    m.add_step(step("Bend", Procedure::Static {}, &["clamp"], &["end"]), None).unwrap();
    m
}

/// A 200 × 200 plate of 2 × 2 quad4 shells in z = 0, clamped along x = 0, under pressure on its
/// top face and its own weight, with a symmetry plane along x = 100.
fn plate() -> Model {
    let mut m = Model::default();
    let n = |i: u32, j: u32| 1 + i + 3 * j;
    for j in 0..=2 {
        for i in 0..=2 {
            m.set_node(n(i, j), [100.0 * f64::from(i), 100.0 * f64::from(j), 0.0]).unwrap();
        }
    }
    let mut e = 1;
    for j in 0..2 {
        for i in 0..2 {
            m.set_element(e, ElementKind::Quad4, vec![n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1)]).unwrap();
            e += 1;
        }
    }
    m.add_set("skin", elements(1..=4)).unwrap();
    m.add_set("top", faces(&[(1, 1), (2, 1), (3, 1), (4, 1)])).unwrap();
    m.add_set("edge", Set { nodes: vec![1, 4, 7], ..Set::default() }).unwrap();
    m.add_set("mid", Set { nodes: vec![2, 5, 8], ..Set::default() }).unwrap();
    let aluminium = Material { elastic: Some(Elastic::Isotropic { young: 70_000.0, poisson: 0.33 }), density: Some(2.7e-9), ..Material::default() };
    m.add_material("Aluminium", aluminium).unwrap();
    let kind = SectionKind::Shell { thickness: 2.0, offset: 0.5, points: 5 };
    let formulation = Formulation { integration: Integration::Reduced, ..Formulation::default() };
    m.add_section("plate", Section { set: "skin".into(), material: "Aluminium".into(), kind, formulation }).unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("edge".into()) }).unwrap();
    m.add_boundary_condition("sym", BoundaryCondition::Symmetry { target: Target::Set("mid".into()), plane: Symmetry::X }).unwrap();
    m.add_load("wind", Load::Pressure { set: "top".into(), magnitude: 0.01 }).unwrap();
    m.add_load("weight", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: None }).unwrap();
    m.add_step(step("Load", Procedure::Static {}, &["clamp", "sym"], &["wind", "weight"]), None).unwrap();
    m
}

/// The cantilever with a reference point coupled to its tip, three static steps that change the
/// active boundary conditions and loads, and a frequency step in between.
fn frequency() -> Model {
    let mut m = cantilever_mesh();
    m.add_reference_point("Tip", ReferencePoint { xyz: [110.0, 5.0, 5.0], label: 0 }).unwrap();
    m.add_coupling("tie", Coupling { reference_point: "Tip".into(), set: "tip".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
    let rp = || Target::ReferencePoint("Tip".into());
    m.add_load("twist", Load::Moment { target: rp(), vector: [5000.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    m.add_load("swell", Load::BodyForce { set: "all".into(), vector: [0.0, 1e-3, 0.0] }).unwrap();
    m.add_boundary_condition("lift", BoundaryCondition::Displacement { target: rp(), values: [None, None, Some(0.5), None, None, None] }).unwrap();
    m.add_step(step("Twist", Procedure::Static {}, &["clamp"], &["twist"]), None).unwrap();
    m.add_step(step("Modes", Procedure::Frequency { modes: 6, min_hz: Some(1.0), max_hz: None }, &["clamp", "lift"], &[]), None).unwrap();
    m.add_step(step("Lift", Procedure::Static {}, &["clamp", "lift"], &["twist", "swell"]), None).unwrap();
    let mut release = step("Release", Procedure::Static {}, &["clamp"], &[]);
    release.output = [OutputVariable::Displacement].into_iter().collect();
    m.add_step(release, None).unwrap();
    m
}

// ---------------------------------------------------------------------------------------------
// A minimal reader of what the writer emits (full deck readers are a later milestone)

#[derive(Default, Debug)]
struct Parsed {
    nodes: BTreeMap<u32, [f64; 3]>,
    elements: BTreeMap<u32, (String, Vec<u32>)>,
    nsets: BTreeMap<String, Vec<u32>>,
    elsets: BTreeMap<String, Vec<u32>>,
    surfaces: BTreeMap<String, Vec<(u32, String)>>,
}

fn param(keyword: &str, name: &str) -> Option<String> {
    keyword.split(',').skip(1).find_map(|p| {
        let (k, v) = p.split_once('=')?;
        (k.trim().eq_ignore_ascii_case(name)).then(|| v.trim().to_string())
    })
}

fn read(deck: &str) -> Parsed {
    let mut p = Parsed::default();
    let mut keyword = String::new();
    let mut pending: Vec<String> = Vec::new();
    for line in deck.lines() {
        if line.starts_with("**") {
            continue;
        }
        if line.starts_with('*') {
            keyword = line.to_string();
            continue;
        }
        let head = keyword.split(',').next().unwrap_or_default().trim().to_ascii_uppercase();
        let fields: Vec<String> = line.split(',').map(|f| f.trim().to_string()).collect();
        let label = |s: &str| s.parse::<u32>().unwrap();
        match head.as_str() {
            "*NODE" => {
                let xyz = [1, 2, 3].map(|i| fields[i].parse::<f64>().unwrap());
                p.nodes.insert(label(&fields[0]), xyz);
            }
            "*ELEMENT" => {
                pending.extend(fields.iter().filter(|f| !f.is_empty()).cloned());
                if !line.trim_end().ends_with(',') {
                    let ids: Vec<u32> = pending.drain(..).map(|f| label(&f)).collect();
                    p.elements.insert(ids[0], (param(&keyword, "TYPE").unwrap(), ids[1..].to_vec()));
                }
            }
            "*NSET" => p.nsets.entry(param(&keyword, "NSET").unwrap()).or_default().extend(fields.iter().map(|f| label(f))),
            "*ELSET" => p.elsets.entry(param(&keyword, "ELSET").unwrap()).or_default().extend(fields.iter().map(|f| label(f))),
            "*SURFACE" if param(&keyword, "TYPE").as_deref() == Some("ELEMENT") => {
                p.surfaces.entry(param(&keyword, "NAME").unwrap()).or_default().push((label(&fields[0]), fields[1].clone()));
            }
            _ => {}
        }
    }
    p
}

// ---------------------------------------------------------------------------------------------
// Golden decks

fn golden(name: &str, deck: &str) {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/abaqus").join(format!("{name}.inp"));
    if std::env::var_os("FIELDCRAFT_BLESS").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, deck).unwrap();
    }
    let expected = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e} (run with FIELDCRAFT_BLESS=1 to write it)", path.display()));
    assert!(expected == deck, "{name}.inp differs from the written deck (FIELDCRAFT_BLESS=1 rewrites it):\n{deck}");
}

fn assert_lines_fit(deck: &str) {
    for (i, line) in deck.lines().enumerate() {
        assert!(line.len() <= abaqus::MAX_LINE, "line {} has {} characters", i + 1, line.len());
        assert!(line.is_ascii(), "line {} is not ASCII", i + 1);
    }
}

#[test]
fn golden_decks() {
    for (name, model) in [("cantilever", cantilever()), ("plate", plate()), ("frequency", frequency())] {
        let deck = abaqus::write(&model).unwrap();
        assert_lines_fit(&deck);
        golden(name, &deck);
        // Deterministic: the same model, and the model after a file round trip, give the same bytes.
        assert_eq!(abaqus::write(&model).unwrap(), deck);
        let reread: Model = serde_json::from_str(&serde_json::to_string(&model).unwrap()).unwrap();
        assert_eq!(abaqus::write(&reread).unwrap(), deck, "{name}");
    }
}

#[test]
fn the_end_force_is_shared_by_area_and_keeps_its_total() {
    let deck = abaqus::write(&cantilever()).unwrap();
    let cload: Vec<&str> = deck.lines().skip_while(|l| !l.starts_with("*CLOAD")).skip(1).take_while(|l| !l.starts_with('*')).collect();
    // Four equal corners of a flat square face.
    assert_eq!(cload, ["17, 3, -25.0", "18, 3, -25.0", "19, 3, -25.0", "20, 3, -25.0"]);
    let total: f64 = cload.iter().map(|l| l.rsplit(", ").next().unwrap().parse::<f64>().unwrap()).sum();
    assert_eq!(total, -100.0);
}

#[test]
fn the_written_mesh_and_sets_read_back() {
    for model in [cantilever(), plate(), frequency()] {
        let deck = abaqus::write(&model).unwrap();
        let p = read(&deck);
        let rps: BTreeMap<u32, [f64; 3]> = model.reference_points().values().map(|r| (r.label, r.xyz)).collect();
        let mut nodes: BTreeMap<u32, [f64; 3]> = model.nodes().iter().map(|(l, n)| (*l, n.xyz)).collect();
        nodes.extend(rps);
        assert_eq!(p.nodes, nodes);
        let connectivity: BTreeMap<u32, Vec<u32>> = p.elements.iter().map(|(l, (_, n))| (*l, n.clone())).collect();
        assert_eq!(connectivity, model.elements().iter().map(|(l, e)| (*l, e.nodes.clone())).collect());
        let r = model.resolve().unwrap();
        for (name, s) in &r.sets {
            assert_eq!(p.nsets.get(name).cloned().unwrap_or_default(), s.nodes, "nset {name}");
            assert_eq!(p.elsets.get(name).cloned().unwrap_or_default(), s.elements, "elset {name}");
            let faces: Vec<(u32, String)> =
                s.faces.iter().map(|f| (f.element, abaqus::face_id(model.element(f.element).unwrap().kind, f.face).unwrap().to_string())).collect();
            assert_eq!(p.surfaces.get(name).cloned().unwrap_or_default(), faces, "surface {name}");
        }
    }
}

#[test]
fn production_reader_preserves_supported_written_model_data() {
    use fieldcraft_model::units::Preset;
    let mut models = vec![cantilever(), plate(), frequency()];
    models.extend(
        [
            ElementKind::Tet4,
            ElementKind::Tet10,
            ElementKind::Wedge6,
            ElementKind::Hex8,
            ElementKind::Hex20,
            ElementKind::Tri3,
            ElementKind::Quad4,
            ElementKind::Quad8,
        ]
        .map(one_element),
    );
    let model_data = |text: &str| {
        let mut out = Vec::new();
        let mut included = false;
        for line in text.lines().filter(|line| !line.starts_with("**")) {
            if line.starts_with("*STEP") {
                break;
            }
            if line.starts_with('*') {
                included = [
                    "*NODE",
                    "*ELEMENT",
                    "*NSET",
                    "*ELSET",
                    "*SURFACE, TYPE=ELEMENT",
                    "*MATERIAL",
                    "*ELASTIC",
                    "*DENSITY",
                    "*SOLID SECTION",
                    "*SHELL SECTION",
                    "*BEAM SECTION",
                ]
                .iter()
                .any(|prefix| line == *prefix || line.starts_with(&format!("{prefix},")));
                // An explicit reference node is an ordinary input NODE; its separate block can merge.
                if line == "*NODE" {
                    continue;
                }
            }
            if included {
                out.push(line.to_string());
            }
        }
        out
    };
    for original in models {
        let written = abaqus::write(&original).unwrap();
        let mut imported = fieldcraft_decks::abaqus::read_in(&written, Preset::MmTonneS).unwrap();
        let si = common::si(&original);
        let mut expected_nodes = si.nodes().clone();
        expected_nodes.extend(si.reference_points().values().map(|p| (p.label, fieldcraft_model::Node { xyz: p.xyz })));
        assert_eq!(&expected_nodes, imported.model.nodes());
        assert_eq!(si.elements(), imported.model.elements());
        assert_eq!(si.materials(), imported.model.materials());
        for section in si.sections().values() {
            assert!(imported.model.sections().values().any(|other| other == section));
        }
        // History is deliberately unread; add a new ordinary step for the writer's model check.
        imported.model.add_step(Step::new("Imported", Procedure::Static {}), None).unwrap();
        let rewritten = fieldcraft_decks::abaqus::write(&imported.model, Preset::MmTonneS).unwrap();
        assert_eq!(model_data(&written), model_data(&rewritten));
    }
}

// ---------------------------------------------------------------------------------------------
// Node order: the element Jacobian of the written connectivity, in Abaqus's own conventions

/// Natural coordinates of each node, by Abaqus element family.
fn natural(kind: ElementKind) -> Vec<[f64; 3]> {
    let hex = [[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.], [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]];
    match kind {
        ElementKind::Hex8 => hex.to_vec(),
        ElementKind::Hex20 => {
            let mut v = hex.to_vec();
            // Mid-sides: the bottom face's edges, the top face's, then the vertical ones.
            v.extend([[0., -1., -1.], [1., 0., -1.], [0., 1., -1.], [-1., 0., -1.]]);
            v.extend([[0., -1., 1.], [1., 0., 1.], [0., 1., 1.], [-1., 0., 1.]]);
            v.extend([[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.]]);
            v
        }
        ElementKind::Tet4 | ElementKind::Tet10 => {
            let mut v = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
            if kind == ElementKind::Tet10 {
                v.extend([[0.5, 0., 0.], [0.5, 0.5, 0.], [0., 0.5, 0.], [0., 0., 0.5], [0.5, 0., 0.5], [0., 0.5, 0.5]]);
            }
            v
        }
        ElementKind::Wedge6 => vec![[0., 0., -1.], [1., 0., -1.], [0., 1., -1.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]],
        ElementKind::Quad4 => vec![[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.]],
        ElementKind::Quad8 => {
            vec![[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.], [0., -1., 0.], [1., 0., 0.], [0., 1., 0.], [-1., 0., 0.]]
        }
        ElementKind::Tri3 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
        _ => Vec::new(),
    }
}

/// Abaqus's shape functions of an element family at natural point `q`.
fn shape(kind: ElementKind, q: [f64; 3]) -> Vec<f64> {
    let nodes = natural(kind);
    let [x, y, z] = q;
    match kind {
        ElementKind::Hex8 => nodes.iter().map(|p| (1. + x * p[0]) * (1. + y * p[1]) * (1. + z * p[2]) / 8.).collect(),
        ElementKind::Hex20 => nodes
            .iter()
            .map(|p| {
                let (a, b, c) = (1. + x * p[0], 1. + y * p[1], 1. + z * p[2]);
                if p[0] == 0. {
                    (1. - x * x) * b * c / 4.
                } else if p[1] == 0. {
                    a * (1. - y * y) * c / 4.
                } else if p[2] == 0. {
                    a * b * (1. - z * z) / 4.
                } else {
                    a * b * c * (x * p[0] + y * p[1] + z * p[2] - 2.) / 8.
                }
            })
            .collect(),
        ElementKind::Tet4 => vec![1. - x - y - z, x, y, z],
        ElementKind::Tet10 => {
            let l = [1. - x - y - z, x, y, z];
            let mut n: Vec<f64> = l.iter().map(|li| li * (2. * li - 1.)).collect();
            for (a, b) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
                n.push(4. * l[a] * l[b]);
            }
            n
        }
        ElementKind::Wedge6 => {
            let t = [1. - x - y, x, y];
            t.iter().map(|ti| ti * (1. - z) / 2.).chain(t.iter().map(|ti| ti * (1. + z) / 2.)).collect()
        }
        ElementKind::Quad4 => nodes.iter().map(|p| (1. + x * p[0]) * (1. + y * p[1]) / 4.).collect(),
        ElementKind::Quad8 => nodes
            .iter()
            .map(|p| {
                if p[0] == 0. {
                    (1. - x * x) * (1. + y * p[1]) / 2.
                } else if p[1] == 0. {
                    (1. + x * p[0]) * (1. - y * y) / 2.
                } else {
                    (1. + x * p[0]) * (1. + y * p[1]) * (x * p[0] + y * p[1] - 1.) / 4.
                }
            })
            .collect(),
        ElementKind::Tri3 => vec![1. - x - y, x, y],
        _ => Vec::new(),
    }
}

/// Columns ∂X/∂ξ, ∂X/∂η, ∂X/∂ζ of the mapping at `q` (central differences).
fn jacobian(kind: ElementKind, xyz: &[[f64; 3]], q: [f64; 3]) -> [[f64; 3]; 3] {
    let map = |q: [f64; 3]| {
        let n = shape(kind, q);
        let mut p = [0.0; 3];
        for (ni, x) in n.iter().zip(xyz) {
            for k in 0..3 {
                p[k] += ni * x[k];
            }
        }
        p
    };
    let h = 1e-6;
    [0, 1, 2].map(|d| {
        let (mut a, mut b) = (q, q);
        a[d] += h;
        b[d] -= h;
        let (pa, pb) = (map(a), map(b));
        [0, 1, 2].map(|k| (pa[k] - pb[k]) / (2. * h))
    })
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The smallest Jacobian determinant over the nodes and the centroid (solids), or the smallest
/// normal component along `up` (shells).
fn min_jacobian(kind: ElementKind, xyz: &[[f64; 3]], up: [f64; 3]) -> f64 {
    let nodes = natural(kind);
    let centroid = [0, 1, 2].map(|k| nodes.iter().map(|p| p[k]).sum::<f64>() / nodes.len() as f64);
    nodes
        .iter()
        .chain([&centroid])
        .map(|q| {
            let [a, b, c] = jacobian(kind, xyz, *q);
            if kind.dimension() == 3 { dot(cross(a, b), c) } else { dot(cross(a, b), up) }
        })
        .fold(f64::INFINITY, f64::min)
}

/// A one-element model of `kind` with a slanted, stretched, positive-volume placement of the
/// model's reference element (mid-side nodes at edge midpoints).
fn one_element(kind: ElementKind) -> Model {
    let corners: Vec<[f64; 3]> = match kind {
        ElementKind::Tri3 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
        ElementKind::Quad4 | ElementKind::Quad8 => vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
        ElementKind::Tet4 | ElementKind::Tet10 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        ElementKind::Wedge6 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]],
        _ => vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]],
    };
    let mut xyz = corners.clone();
    xyz.resize(kind.node_count(), [0.0; 3]);
    for e in kind.edges() {
        if let [a, b, m] = **e {
            xyz[m] = [0, 1, 2].map(|i| (corners[a][i] + corners[b][i]) / 2.0);
        }
    }
    // A shear and stretch with a positive determinant, then a shift.
    let place = |p: [f64; 3]| [3. * p[0] + 0.5 * p[1] + 10., 2. * p[1] + 0.3 * p[2] - 4., 1.5 * p[2] + 0.2 * p[0] + 7.];
    let mut m = Model::default();
    for (i, p) in xyz.iter().enumerate() {
        m.set_node(i as u32 + 1, place(*p)).unwrap();
    }
    m.set_element(1, kind, (1..=kind.node_count() as u32).collect()).unwrap();
    m.add_set("all", elements([1])).unwrap();
    m.add_set("base", Set { nodes: vec![1], ..Set::default() }).unwrap();
    m.add_material("Steel", steel()).unwrap();
    let (section, formulation) = if kind.dimension() == 3 {
        (SectionKind::Solid {}, Formulation::default())
    } else {
        (SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 }, Formulation { integration: Integration::Reduced, ..Formulation::default() })
    };
    m.add_section("s", Section { set: "all".into(), material: "Steel".into(), kind: section, formulation }).unwrap();
    m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    m.add_step(step("Modes", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None }, &["fix"], &[]), None).unwrap();
    m
}

#[test]
fn written_connectivity_has_a_positive_jacobian_in_abaqus_conventions() {
    let kinds = [ElementKind::Tet4, ElementKind::Tet10, ElementKind::Wedge6, ElementKind::Hex8, ElementKind::Hex20];
    let shells = [ElementKind::Tri3, ElementKind::Quad4, ElementKind::Quad8];
    for kind in kinds.into_iter().chain(shells) {
        let model = one_element(kind);
        let p = read(&abaqus::write(&model).unwrap());
        let (_, nodes) = &p.elements[&1];
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|n| p.nodes[n]).collect();
        // A shell's positive normal is the side its face 1 (SPOS) looks to: the placed +Z here.
        let up = cross([3., 0., 0.2], [0.5, 2., 0.]);
        let j = min_jacobian(kind, &xyz, up);
        assert!(j > 0.0, "{kind:?}: smallest Jacobian {j}");

        // The test has teeth: a mirrored linear element shows up (and a misplaced quadratic
        // mid-side block, below).
        if kind.node_count() == natural(kind).len()
            && matches!(kind, ElementKind::Tet4 | ElementKind::Wedge6 | ElementKind::Hex8 | ElementKind::Tri3 | ElementKind::Quad4)
        {
            let mut mirrored = xyz.clone();
            match kind {
                ElementKind::Tet4 | ElementKind::Tri3 => mirrored.swap(1, 2),
                ElementKind::Wedge6 => {
                    mirrored.swap(1, 2);
                    mirrored.swap(4, 5);
                }
                ElementKind::Quad4 => mirrored.swap(1, 3),
                _ => {
                    mirrored.swap(1, 3);
                    mirrored.swap(5, 7);
                }
            }
            assert!(min_jacobian(kind, &mirrored, up) < 0.0, "{kind:?}: mirrored");
        }
    }
    // Hex20 in the order with the vertical edges' mid-sides before the top face's (Nastran's).
    let model = one_element(ElementKind::Hex20);
    let p = read(&abaqus::write(&model).unwrap());
    let mut xyz: Vec<[f64; 3]> = p.elements[&1].1.iter().map(|n| p.nodes[n]).collect();
    let (top, vertical) = (xyz[12..16].to_vec(), xyz[16..20].to_vec());
    xyz[12..16].copy_from_slice(&vertical);
    xyz[16..20].copy_from_slice(&top);
    assert!(min_jacobian(ElementKind::Hex20, &xyz, [0.; 3]) < 0.0, "a misplaced mid-side block is caught");
}

#[test]
fn model_faces_are_abaqus_faces() {
    // Each face's node positions (1-based, Abaqus element library, our transcription).
    let table: [(ElementKind, &[&[usize]]); 4] = [
        (ElementKind::Hex8, &[&[1, 2, 3, 4], &[5, 8, 7, 6], &[1, 5, 6, 2], &[2, 6, 7, 3], &[3, 7, 8, 4], &[4, 8, 5, 1]]),
        (ElementKind::Tet4, &[&[1, 2, 3], &[1, 4, 2], &[2, 4, 3], &[3, 4, 1]]),
        (ElementKind::Wedge6, &[&[1, 2, 3], &[4, 6, 5], &[1, 4, 5, 2], &[2, 5, 6, 3], &[3, 6, 4, 1]]),
        (ElementKind::Pyramid5, &[&[1, 2, 3, 4], &[1, 5, 2], &[2, 5, 3], &[3, 5, 4], &[4, 5, 1]]),
    ];
    for (kind, faces) in table {
        assert_eq!(kind.faces().len(), faces.len());
        for (i, (ours, theirs)) in kind.faces().iter().zip(faces).enumerate() {
            let ours: BTreeSet<usize> = ours.iter().map(|p| p + 1).collect();
            let theirs: BTreeSet<usize> = theirs.iter().copied().collect();
            assert_eq!(ours, theirs, "{kind:?} face {}", i + 1);
            assert_eq!(abaqus::face_id(kind, i as u8 + 1), Some(format!("S{}", i + 1).as_str()));
        }
        assert_eq!(abaqus::face_id(kind, 0), None);
        assert_eq!(abaqus::face_id(kind, faces.len() as u8 + 1), None);
    }
    assert_eq!(abaqus::face_id(ElementKind::Quad4, 1), Some("SPOS"));
    assert_eq!(abaqus::face_id(ElementKind::Tri3, 2), Some("SNEG"));
    assert_eq!(abaqus::face_id(ElementKind::Beam2, 1), None);
}

// ---------------------------------------------------------------------------------------------
// Numbers and limits

#[test]
fn numbers_read_back_exactly() {
    for (v, s) in [(0.1, "0.1"), (210_000.0, "210000.0"), (7.85e-9, "7.85E-9"), (1e300, "1.0E300"), (-0.5, "-0.5"), (1e16, "1.0E16")] {
        assert_eq!(abaqus::number(v).unwrap(), s);
    }
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..100_000 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        let v = f64::from_bits(x);
        if !v.is_finite() {
            assert!(abaqus::number(v).is_err());
            continue;
        }
        let s = abaqus::number(v).unwrap();
        assert_eq!(s.parse::<f64>().unwrap().to_bits(), v.to_bits(), "{v:e} → {s}");
        assert!(s.len() <= 24, "{s}");
        assert!(!s.contains('e'), "{s}");
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(abaqus::number(bad), Err(Error::Unsupported { .. })));
    }
}

#[test]
fn largest_labels_longest_numbers_and_names_fit_a_line() {
    let mut m = Model::default();
    let base = fieldcraft_model::MAX_LABEL - 30;
    let ugly = -1.234_567_890_123_456_7e-300;
    let mut xyz = one_element(ElementKind::Hex20).nodes().values().map(|n| n.xyz).collect::<Vec<_>>();
    xyz[0] = [ugly, ugly, ugly];
    for (i, p) in xyz.iter().enumerate() {
        m.set_node(base + i as u32, *p).unwrap();
    }
    m.set_element(fieldcraft_model::MAX_LABEL, ElementKind::Hex20, (0..20).map(|i| base + i).collect()).unwrap();
    let long = format!("A{}", "x".repeat(fieldcraft_model::MAX_NAME_LEN - 1));
    m.add_set(&long, Set { elements: vec![fieldcraft_model::MAX_LABEL], nodes: (0..20).map(|i| base + i).collect(), ..Set::default() }).unwrap();
    m.add_material(
        &long,
        Material { elastic: Some(Elastic::Isotropic { young: ugly.abs(), poisson: 0.3 }), density: Some(1.0), ..Material::default() },
    )
    .unwrap();
    m.add_section(&long, Section { set: long.clone(), material: long.clone(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.add_reference_point(&long, ReferencePoint { xyz: [ugly; 3], label: 0 }).unwrap();
    m.add_coupling(&long, Coupling { reference_point: long.clone(), set: long.clone(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
    m.add_boundary_condition(&long, BoundaryCondition::Fixed { target: Target::Set(long.clone()) }).unwrap();
    m.add_step(step(&long, Procedure::Frequency { modes: 100_000, min_hz: Some(ugly.abs()), max_hz: Some(1e300) }, &[&long], &[]), None).unwrap();
    let deck = abaqus::write(&m).unwrap();
    assert_lines_fit(&deck);
    let p = read(&deck);
    assert_eq!(p.elements[&fieldcraft_model::MAX_LABEL].1.len(), 20, "a hex20 continues onto a second line");
    assert_eq!(p.nodes[&base], [ugly; 3]);
}

// ---------------------------------------------------------------------------------------------
// Step activity, load summing

fn block<'a>(deck: &'a str, step: &str) -> Vec<&'a str> {
    let start = format!("*STEP, NAME={step}");
    deck.lines().skip_while(|l| !l.starts_with(&start)).take_while(|l| *l != "*END STEP").filter(|l| !l.starts_with("**")).collect()
}

#[test]
fn steps_write_the_activity_table() {
    let deck = abaqus::write(&frequency()).unwrap();
    // `clamp` holds at zero in every step: the initial state.
    let initial: Vec<&str> = deck.lines().skip_while(|l| !l.starts_with("** Initial state")).skip(1).take(3).collect();
    assert_eq!(initial, ["*BOUNDARY", "** boundary condition clamp", "root, ENCASTRE"]);
    let twist = block(&deck, "Twist");
    assert!(twist.contains(&"*CLOAD, OP=NEW") && twist.contains(&"99999999, 4, 5000.0"), "{twist:?}");
    assert!(!twist.iter().any(|l| l.starts_with("*BOUNDARY")), "{twist:?}");
    // The frequency step adds `lift` to its base state and writes no loads.
    let modes = block(&deck, "Modes");
    assert!(modes.contains(&"*BOUNDARY") && modes.contains(&"99999999, 3, 3, 0.5"), "{modes:?}");
    assert!(!modes.iter().any(|l| l.contains("LOAD")), "{modes:?}");
    assert!(modes.contains(&"6, 1.0"), "{modes:?}");
    // `Lift` adds `lift` again: the frequency step's change did not carry on.
    let lift = block(&deck, "Lift");
    assert!(lift.contains(&"*BOUNDARY") && lift.contains(&"99999999, 3, 3, 0.5"), "{lift:?}");
    assert!(lift.contains(&"*DLOAD, OP=NEW") && lift.contains(&"all, BY, 0.001"), "{lift:?}");
    // `Release` drops `lift` and every load.
    let release = block(&deck, "Release");
    assert!(release.contains(&"*BOUNDARY, OP=NEW") && release.contains(&"root, ENCASTRE"), "{release:?}");
    for k in ["*CLOAD, OP=NEW", "*DLOAD, OP=NEW"] {
        let at = release.iter().position(|l| *l == k).unwrap_or_else(|| panic!("{k} in {release:?}"));
        assert!(release.get(at + 1).is_none_or(|l| l.starts_with('*')), "{k} restates nothing: {release:?}");
    }
    assert!(!release.contains(&"*DSLOAD, OP=NEW"), "no pressure before, none now");
    assert!(release.contains(&"U") && !release.contains(&"*ELEMENT OUTPUT"), "{release:?}");
}

#[test]
fn loads_sharing_nodes_elements_or_surfaces_add_up() {
    let mut m = plate();
    m.add_load("gust", Load::Pressure { set: "top".into(), magnitude: 0.02 }).unwrap();
    m.add_set("half", elements([1, 2])).unwrap();
    m.add_load("extra", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: Some("half".into()) }).unwrap();
    m.add_load("p1", Load::Force { target: Target::Set("mid".into()), vector: [3.0, 0.0, 0.0], distribution: Distribution::PerNode }).unwrap();
    m.add_load("p2", Load::Force { target: Target::Set("mid".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::PerNode }).unwrap();
    for l in ["gust", "extra", "p1", "p2"] {
        m.set_load_active("Load", l, true).unwrap();
    }
    let deck = abaqus::write(&m).unwrap();
    let lines = block(&deck, "Load");
    assert!(lines.contains(&"top, P, 0.03"), "{lines:?}");
    assert!(lines.contains(&"5, 1, 4.0"), "{lines:?}");
    // Two gravity loads on elements 1 and 2: summed per element.
    assert!(lines.contains(&"1, GRAV, 19620.0, 0.0, 0.0, -1.0"), "{lines:?}");
    assert!(lines.contains(&"3, GRAV, 9810.0, 0.0, 0.0, -1.0"), "{lines:?}");
    assert!(!lines.iter().any(|l| l.starts_with("skin, GRAV")), "{lines:?}");
}

// ---------------------------------------------------------------------------------------------
// What the format cannot express, and models that are not ready

fn unsupported(m: &Model) -> String {
    match abaqus::write(m) {
        Err(Error::Unsupported { format: "Abaqus", what }) => what,
        other => panic!("expected unsupported, got {other:?}"),
    }
}

fn with_section(kind: ElementKind, section: SectionKind, integration: Integration) -> Model {
    let mut m = one_element(kind);
    let formulation = Formulation { integration, ..Formulation::default() };
    m.replace_section("s", Section { set: "all".into(), material: "Steel".into(), kind: section, formulation }).unwrap();
    m
}

#[test]
fn what_abaqus_cannot_express_is_named() {
    let shell = SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 };
    let mut tri6 = Model::default();
    for (i, p) in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0.5, 0., 0.], [0.5, 0.5, 0.], [0., 0.5, 0.]].iter().enumerate() {
        tri6.set_node(i as u32 + 1, *p).unwrap();
    }
    tri6.set_element(1, ElementKind::Tri6, (1..=6).collect()).unwrap();
    tri6.add_set("all", elements([1])).unwrap();
    tri6.add_material("Steel", steel()).unwrap();
    tri6.add_section("s", Section { set: "all".into(), material: "Steel".into(), kind: shell, formulation: Formulation::default() }).unwrap();
    tri6.add_step(step("Modes", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None }, &[], &[]), None).unwrap();
    assert!(unsupported(&tri6).contains("tri6 shells"));

    assert!(unsupported(&with_section(ElementKind::Quad8, shell, Integration::Full)).contains("full integration of quad8"));
    assert!(unsupported(&with_section(ElementKind::Tet10, SectionKind::Solid {}, Integration::Reduced)).contains("reduced integration of tet10"));
    assert!(unsupported(&with_section(ElementKind::Wedge6, SectionKind::Solid {}, Integration::Reduced)).contains("wedge6"));
    // A linear tetrahedron integrates at one point either way.
    assert!(abaqus::write(&with_section(ElementKind::Tet4, SectionKind::Solid {}, Integration::Reduced)).unwrap().contains("TYPE=C3D4\n"));

    let mut case = cantilever();
    case.add_set("ROOT", Set { nodes: vec![1], ..Set::default() }).unwrap();
    let what = unsupported(&case);
    assert!(what.contains("sets `ROOT` and `root`") || what.contains("sets `root` and `ROOT`"), "{what}");

    let mut ortho = plate();
    let carbon = Elastic::EngineeringConstants { young: [140e3, 10e3, 10e3], poisson: [0.3, 0.3, 0.4], shear: [5e3, 5e3, 3.8e3] };
    ortho.replace_material("Aluminium", Material { elastic: Some(carbon), density: Some(1.6e-9), ..Material::default() }).unwrap();
    assert!(unsupported(&ortho).contains("orthotropic material `Aluminium` in a shell section"));

    let mut loose = frequency();
    let dofs = Dofs::try_from(vec![1, 2, 4, 5, 6]).unwrap();
    loose.replace_coupling("tie", Coupling { reference_point: "Tip".into(), set: "tip".into(), kind: CouplingKind::Distributing, dofs }).unwrap();
    assert!(unsupported(&loose).contains("coupling `tie`"));

    // Several problems are named together.
    let mut both = with_section(ElementKind::Quad8, shell, Integration::Full);
    both.add_set("ALL", Set { nodes: vec![1], ..Set::default() }).unwrap();
    assert!(abaqus::write(&both).unwrap_err().to_string().contains("full integration"));
}

#[test]
fn orthotropic_solids_beams_and_rods_are_written() {
    let mut m = cantilever();
    let carbon = Elastic::EngineeringConstants { young: [140e3, 10e3, 10e3], poisson: [0.3, 0.3, 0.4], shear: [5e3, 5e3, 3.8e3] };
    m.replace_material("Steel", Material { elastic: Some(carbon), density: None, ..Material::default() }).unwrap();
    let deck = abaqus::write(&m).unwrap();
    assert!(deck.contains("*ELASTIC, TYPE=ENGINEERING CONSTANTS\n140000.0, 10000.0, 10000.0, 0.3, 0.3, 0.4, 5000.0, 5000.0\n3800.0\n"), "{deck}");

    // A frame of one beam and one rod.
    let mut f = Model::default();
    for (i, p) in [[0., 0., 0.], [1000., 0., 0.], [1000., 0., 500.]].iter().enumerate() {
        f.set_node(i as u32 + 1, *p).unwrap();
    }
    f.set_element(1, ElementKind::Beam2, vec![1, 2]).unwrap();
    f.set_element(2, ElementKind::Rod2, vec![2, 3]).unwrap();
    f.add_set("beam", elements([1])).unwrap();
    f.add_set("rod", elements([2])).unwrap();
    f.add_set("base", Set { nodes: vec![1, 3], ..Set::default() }).unwrap();
    f.add_material("Steel", steel()).unwrap();
    let i = BeamProfile::I { height: 100.0, width: 50.0, web: 5.0, flange: 8.0 };
    let beam = SectionKind::Beam { profile: i, orientation: [0.0, 1.0, 0.0] };
    f.add_section("b", Section { set: "beam".into(), material: "Steel".into(), kind: beam, formulation: Formulation::default() }).unwrap();
    let rod = SectionKind::Rod { area: 50.0 };
    f.add_section("r", Section { set: "rod".into(), material: "Steel".into(), kind: rod, formulation: Formulation::default() }).unwrap();
    f.add_boundary_condition("pin", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    f.add_step(step("Modes", Procedure::Frequency { modes: 3, min_hz: None, max_hz: Some(500.0) }, &["pin"], &[]), None).unwrap();
    let deck = abaqus::write(&f).unwrap();
    for expected in [
        "*ELEMENT, TYPE=B31\n1, 1, 2\n",
        "*ELEMENT, TYPE=T3D2\n2, 2, 3\n",
        "*BEAM SECTION, SECTION=I, ELSET=beam, MATERIAL=Steel\n50.0, 100.0, 50.0, 50.0, 8.0, 8.0, 5.0\n0.0, 1.0, 0.0\n",
        "*SOLID SECTION, ELSET=rod, MATERIAL=Steel\n50.0\n",
        "3, , 500.0\n",
    ] {
        assert!(deck.contains(expected), "{expected}\n{deck}");
    }
}

#[test]
fn a_model_that_is_not_ready_is_refused() {
    let mut m = cantilever();
    m.remove_section("beam").unwrap();
    assert!(matches!(abaqus::write(&m), Err(Error::Model(_))));
    assert!(matches!(abaqus::write(&Model::default()), Err(Error::Model(_))));
}

/// The cantilever with its steel given each plasticity in turn.
fn plastic_cantilever(p: fieldcraft_model::Plastic) -> Model {
    let mut m = cantilever();
    m.replace_material("Steel", Material { plastic: Some(p), ..steel() }).unwrap();
    m
}

/// Plasticity is written as `*PLASTIC` (and Johnson–Cook's rate and heat keywords); damage and
/// erosion, failure models of Abaqus/Explicit, are reported.
#[test]
fn plasticity_is_written_and_failure_reported() {
    use fieldcraft_model::{Damage, Erosion, Hardening, HardeningPoint, Plastic};
    let section = |deck: &str| deck.split("*MATERIAL, NAME=Steel").nth(1).and_then(|s| s.split("*SOLID SECTION").next()).unwrap().to_string();
    let deck =
        abaqus::write(&plastic_cantilever(Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 1000.0, hardening: Hardening::Isotropic }))
            .unwrap();
    assert!(section(&deck).contains("*PLASTIC\n250.0, 0.0\n10250.0, 10.0\n"), "{}", section(&deck));
    let deck =
        abaqus::write(&plastic_cantilever(Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 1000.0, hardening: Hardening::Kinematic }))
            .unwrap();
    assert!(section(&deck).contains("*PLASTIC, HARDENING=KINEMATIC\n250.0, 0.0\n1250.0, 1.0\n"), "{}", section(&deck));
    let table = Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: 200.0 }, HardeningPoint { strain: 0.1, stress: 300.0 }] };
    let deck = abaqus::write(&plastic_cantilever(table)).unwrap();
    assert!(section(&deck).contains("*PLASTIC\n200.0, 0.0\n300.0, 0.1\n"), "{}", section(&deck));
    let jc = Plastic::JohnsonCook {
        a: 90.0,
        b: 292.0,
        n: 0.31,
        c: 0.025,
        m: 1.09,
        reference_strain_rate: 1.0,
        room_temperature: 294.0,
        melt_temperature: 1356.0,
        specific_heat: Some(3.83e8),
        taylor_quinney: 0.9,
    };
    let deck = abaqus::write(&plastic_cantilever(jc.clone())).unwrap();
    let s = section(&deck);
    assert!(
        s.contains("*PLASTIC, HARDENING=JOHNSON COOK\n90.0, 292.0, 0.31, 1.09, 1356.0, 294.0\n*RATE DEPENDENT, TYPE=JOHNSON COOK\n0.025, 1.0\n*SPECIFIC HEAT\n383000000.0\n*INELASTIC HEAT FRACTION\n0.9\n"),
        "{s}"
    );
    let mut failing = plastic_cantilever(jc);
    let mut steel = failing.material("Steel").unwrap().clone();
    steel.damage = Some(Damage::JohnsonCook { d1: 0.54, d2: 4.89, d3: -3.03, d4: 0.014, d5: 1.12 });
    steel.erosion = Some(Erosion { plastic_strain: 1.0 });
    failing.replace_material("Steel", steel).unwrap();
    let Err(Error::Unsupported { what, .. }) = abaqus::write(&failing) else { panic!("damage is not written") };
    assert!(what.contains("material `Steel`: damage and erosion"), "{what}");
}

/// A nonlinear static step is a general step with NLGEOM, its increments and its line search; a
/// neo-Hookean material is then `*HYPERELASTIC` and a linear static step beside it is refused.
#[test]
fn nonlinear_steps_and_hyperelastic_materials_are_written() {
    use fieldcraft_model::{Hyperelastic, NonlinearStatic};
    let mut m = cantilever();
    let controls = NonlinearStatic {
        initial_increment: 0.05,
        min_increment: 1e-4,
        max_increment: 0.5,
        max_increments: 40,
        line_search: true,
        ..NonlinearStatic::default()
    };
    m.replace_step(step("Bend", Procedure::NonlinearStatic(controls.clone()), &["clamp"], &["end"])).unwrap();
    let deck = abaqus::write(&m).unwrap();
    let bend = block(&deck, "Bend");
    assert_eq!(bend.first(), Some(&"*STEP, NAME=Bend, NLGEOM=YES, INC=40"), "{bend:?}");
    let at = bend.iter().position(|l| *l == "*STATIC").unwrap();
    assert_eq!(bend.get(at + 1..at + 4), Some(&["0.05, 1.0, 0.0001, 0.5", "*CONTROLS, PARAMETERS=LINE SEARCH", "5,"][..]), "{bend:?}");
    assert!(deck.contains("*ELASTIC, TYPE=ISOTROPIC\n210000.0, 0.3\n"), "St Venant–Kirchhoff is the elastic material");

    let small = NonlinearStatic { nonlinear_geometry: false, line_search: false, ..controls };
    m.replace_step(step("Bend", Procedure::NonlinearStatic(small), &["clamp"], &["end"])).unwrap();
    m.replace_material("Steel", Material { hyperelastic: Some(Hyperelastic::NeoHookean {}), ..steel() }).unwrap();
    let deck = abaqus::write(&m).unwrap();
    let bend = block(&deck, "Bend");
    assert_eq!(bend.first(), Some(&"*STEP, NAME=Bend, NLGEOM=NO, INC=40"), "{bend:?}");
    assert!(!bend.iter().any(|l| l.starts_with("*CONTROLS")), "{bend:?}");
    // μ = E/(2(1 + ν)) = 80769.23…, κ = E/(3(1 − 2ν)) = 175000: C10 = μ/2, D1 = 2/κ.
    let (c10, d1) = (210_000.0 / 2.6 / 2.0, 2.0 / 175_000.0);
    let line = deck.split("*HYPERELASTIC, NEO HOOKE\n").nth(1).and_then(|s| s.lines().next()).unwrap();
    let values: Vec<f64> = line.split(", ").map(|v| v.parse().unwrap()).collect();
    assert!((values[0] - c10).abs() < 1e-9 * c10 && (values[1] - d1).abs() < 1e-9 * d1, "{line}");
    assert!(!deck.contains("*ELASTIC"), "{deck}");

    // Mooney–Rivlin and Yeoh: C10 = μ/2 − C01, C01, D1 and C10, C20, C30, D1.
    m.replace_material("Steel", Material { hyperelastic: Some(Hyperelastic::MooneyRivlin { c01: 10000.0 }), ..steel() }).unwrap();
    let deck = abaqus::write(&m).unwrap();
    let line = deck.split("*HYPERELASTIC, MOONEY-RIVLIN\n").nth(1).and_then(|s| s.lines().next()).unwrap();
    let values: Vec<f64> = line.split(", ").map(|v| v.parse().unwrap()).collect();
    assert!((values[0] - (c10 - 10000.0)).abs() < 1e-9 * c10 && values[1] == 10000.0 && (values[2] - d1).abs() < 1e-9 * d1, "{line}");
    m.replace_material("Steel", Material { hyperelastic: Some(Hyperelastic::Yeoh { c20: -50.0, c30: 2.0 }), ..steel() }).unwrap();
    let deck = abaqus::write(&m).unwrap();
    let line = deck.split("*HYPERELASTIC, YEOH\n").nth(1).and_then(|s| s.lines().next()).unwrap();
    let values: Vec<f64> = line.split(", ").map(|v| v.parse().unwrap()).collect();
    assert!((values[0] - c10).abs() < 1e-9 * c10 && values[1..3] == [-50.0, 2.0] && (values[3] - d1).abs() < 1e-9 * d1, "{line}");
    m.replace_material("Steel", Material { hyperelastic: Some(Hyperelastic::NeoHookean {}), ..steel() }).unwrap();

    m.add_step(step("Linear", Procedure::Static {}, &["clamp"], &["end"]), None).unwrap();
    assert!(unsupported(&m).contains("step `Linear`: a linear static step beside nonlinear steps of hyperelastic material `Steel`"));
    m.remove_step("Bend").unwrap();
    assert!(abaqus::write(&m).unwrap().contains("*ELASTIC, TYPE=ISOTROPIC\n"), "with linear steps only, the elastic constants are the whole of it");
}

/// A plastic material in nonlinear steps is `*PLASTIC` in a general step; what Abaqus would
/// read differently is refused: a load path, and several general steps (Abaqus carries the
/// plastic state between them).
#[test]
fn nonlinear_plastic_steps_and_what_is_refused() {
    use fieldcraft_model::{Hardening, NonlinearStatic, Plastic};
    let mut m = cantilever();
    let plastic = Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 2000.0, hardening: Hardening::Isotropic };
    m.replace_material("Steel", Material { plastic: Some(plastic), ..steel() }).unwrap();
    m.replace_step(step("Bend", Procedure::NonlinearStatic(NonlinearStatic::default()), &["clamp"], &["end"])).unwrap();
    let deck = abaqus::write(&m).unwrap();
    assert!(deck.contains("*PLASTIC\n250.0, 0.0\n") && deck.contains("NLGEOM=YES"), "{deck}");
    let path = NonlinearStatic { load_path: Some(vec![[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]), ..NonlinearStatic::default() };
    m.replace_step(step("Bend", Procedure::NonlinearStatic(path), &["clamp"], &["end"])).unwrap();
    assert!(unsupported(&m).contains("step `Bend`: a load path"));
    let own = NonlinearStatic { load_paths: [("end".to_string(), vec![[0.0, 0.0], [1.0, 1.0]])].into(), ..NonlinearStatic::default() };
    m.replace_step(step("Bend", Procedure::NonlinearStatic(own), &["clamp"], &["end"])).unwrap();
    assert!(unsupported(&m).contains("step `Bend`: load `end`'s own load path"));
    m.replace_step(step("Bend", Procedure::NonlinearStatic(NonlinearStatic::default()), &["clamp"], &["end"])).unwrap();
    m.add_step(step("Again", Procedure::NonlinearStatic(NonlinearStatic::default()), &["clamp"], &["end"]), None).unwrap();
    assert!(unsupported(&m).contains("plastic material `Steel` in a deck of 2 general steps"));
}

/// A buckling step is `*BUCKLE` with its loads as the perturbation loads; after a general step
/// (Abaqus's base state would be that step's) and with a shift it is refused by name.
#[test]
fn buckling_steps() {
    let mut m = cantilever_mesh();
    m.add_load("end", Load::Force { target: Target::Set("tip".into()), vector: [-100.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    m.add_step(step("Buckle", Procedure::Buckle { modes: 4, shift: None }, &["clamp"], &["end"]), None).unwrap();
    let deck = abaqus::write(&m).unwrap();
    assert!(deck.contains("*BUCKLE, EIGENSOLVER=LANCZOS\n4\n") && deck.contains("*CLOAD"), "{deck}");
    let mut shifted = m.clone();
    shifted.replace_step(step("Buckle", Procedure::Buckle { modes: 4, shift: Some(2.0) }, &["clamp"], &["end"])).unwrap();
    assert!(abaqus::write(&shifted).unwrap_err().to_string().contains("shift"));
    m.add_step(step("Bend", Procedure::Static {}, &["clamp"], &["end"]), Some("Buckle")).unwrap();
    assert!(abaqus::write(&m).unwrap_err().to_string().contains("after a general step"));
    assert!(lsdyna::write(&m).is_err() && nastran::write(&m).is_err());
}

/// An arc-length step is `*STATIC, RIKS` with the arcs, the maximum load factor and the
/// displacement that end it; a spherical constraint is refused, and LS-DYNA refuses the step.
#[test]
fn arc_length_steps() {
    use fieldcraft_model::{ArcConstraint, ArcLength, DisplacementLimit, NonlinearStatic};
    let mut m = cantilever_mesh();
    m.add_load("end", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, -100.0], distribution: Distribution::Total }).unwrap();
    let arc = ArcLength {
        max_load_factor: 3.0,
        displacement_limit: Some(DisplacementLimit { node: 2, direction: 3, displacement: -5.0 }),
        ..ArcLength::default()
    };
    let settings = NonlinearStatic { arc_length: Some(arc.clone()), ..NonlinearStatic::default() };
    m.add_step(step("Snap", Procedure::NonlinearStatic(settings), &["clamp"], &["end"]), None).unwrap();
    let deck = abaqus::write(&m).unwrap();
    assert!(deck.contains("*STATIC, RIKS\n0.1, 1.0, 1.0E-5, 1.0, 3.0, 2, 3, -5.0\n"), "{deck}");
    assert!(lsdyna::write(&m).is_err());
    let spherical = NonlinearStatic { arc_length: Some(ArcLength { constraint: ArcConstraint::Spherical, ..arc }), ..NonlinearStatic::default() };
    m.replace_step(step("Snap", Procedure::NonlinearStatic(spherical), &["clamp"], &["end"])).unwrap();
    assert!(abaqus::write(&m).is_err());
}
