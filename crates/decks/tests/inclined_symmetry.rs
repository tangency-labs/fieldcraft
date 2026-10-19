//! Symmetry about a plane not normal to X, Y or Z in the three decks (#242): Abaqus `*EQUATION`s,
//! Nastran `MPC` sets and an LS-DYNA `*BOUNDARY_SPC` in a local coordinate system, on a small
//! shell plate of our own whose diagonal x = y is the plane.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::Error;
use fieldcraft_model::{
    BoundaryCondition, Elastic, ElementFace, ElementKind, Formulation, Load, Material, Model, Procedure, Section, SectionKind, Set, Step, Symmetry,
    Target,
};

/// A 2 × 2 quad4 shell plate in the XY plane (nodes 1–9, x and y in 0, 1, 2), clamped on x = 0,
/// the inclined plane (normal (1, −1, 0)) on the diagonal nodes 1, 5 and 9, under pressure;
/// steps `Load` (static) and `Modes`, `in_modes` saying whether the plane is active in both.
fn plate(in_modes: bool) -> Model {
    let mut m = Model::default();
    for j in 0..3u32 {
        for i in 0..3u32 {
            m.set_node(1 + i + 3 * j, [f64::from(i), f64::from(j), 0.0]).unwrap();
        }
    }
    for (e, n) in [(1, [1, 2, 5, 4]), (2, [2, 3, 6, 5]), (3, [4, 5, 8, 7]), (4, [5, 6, 9, 8])] {
        m.set_element(e, ElementKind::Quad4, n.to_vec()).unwrap();
    }
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 210e3, poisson: 0.3 }), density: Some(7.85e-9), ..Material::default() },
    )
    .unwrap();
    m.add_set("all", Set { elements: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    let kind = SectionKind::Shell { thickness: 0.1, offset: 0.0, points: 5 };
    m.add_section("plate", Section { set: "all".into(), material: "Steel".into(), kind, formulation: Formulation::default() }).unwrap();
    m.add_set("edge", Set { nodes: vec![1, 4, 7], ..Set::default() }).unwrap();
    m.add_set("diagonal", Set { nodes: vec![1, 5, 9], ..Set::default() }).unwrap();
    let faces = (1..=4).map(|element| ElementFace { element, face: 1 }).collect();
    m.add_set("top", Set { faces, ..Set::default() }).unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("edge".into()) }).unwrap();
    let plane = Symmetry::from_normal([1.0, -1.0, 0.0]).unwrap();
    m.add_boundary_condition("sym", BoundaryCondition::Symmetry { target: Target::Set("diagonal".into()), plane }).unwrap();
    m.add_load("q", Load::Pressure { set: "top".into(), magnitude: 0.01 }).unwrap();
    let mut load = Step::new("Load", Procedure::Static {});
    load.boundary_conditions = ["clamp", "sym"].map(String::from).into();
    load.loads = ["q".to_string()].into();
    m.add_step(load, None).unwrap();
    let mut modes = Step::new("Modes", Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
    modes.boundary_conditions = if in_modes { ["clamp", "sym"].map(String::from).into() } else { ["clamp".to_string()].into() };
    m.add_step(modes, None).unwrap();
    m
}

/// A Nastran real: an exponent may come without `E`.
fn real(s: &str) -> f64 {
    let split = s.char_indices().skip(1).find(|(i, c)| (*c == '+' || *c == '-') && !s[..*i].ends_with(['E', 'e']));
    match split {
        Some((i, _)) => format!("{}E{}", &s[..i], &s[i..]).parse().unwrap(),
        None => s.parse().unwrap(),
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| a[k] * b[k]).sum()
}

/// Every equation is one of the plane's at a free node: its terms on one node, its coefficients
/// along the normal (translations) or along an in-plane axis (rotations), unit length; nodes 5
/// and 9 get three each (n · u, t₁ · θ, t₂ · θ); node 1 is clamped, so its equations are left out;
/// no dependent degree of freedom repeats.
fn check(equations: &[(u32, Vec<(u8, f64)>)]) {
    let n = [1.0 / 2f64.sqrt(), -1.0 / 2f64.sqrt(), 0.0];
    let mut dependents = std::collections::BTreeSet::new();
    let mut per_node = std::collections::BTreeMap::<u32, usize>::new();
    for (node, terms) in equations {
        assert!(dependents.insert((*node, terms[0].0)), "{node} {terms:?}");
        *per_node.entry(*node).or_default() += 1;
        let mut v = [0.0; 3];
        let rotation = terms[0].0 > 3;
        for (dof, c) in terms {
            assert_eq!(*dof > 3, rotation, "{terms:?}");
            v[usize::from((dof - 1) % 3)] = *c;
        }
        assert!((dot(v, v) - 1.0).abs() < 1e-12, "{terms:?}");
        if rotation {
            assert!(dot(v, n).abs() < 1e-12, "{terms:?}");
        } else {
            assert!((dot(v, n).abs() - 1.0).abs() < 1e-12, "{terms:?}");
        }
    }
    assert_eq!(per_node, [(5, 3), (9, 3)].into());
}

#[test]
fn abaqus_writes_equations_in_the_model_data() {
    let deck = abaqus::write(&plate(true)).unwrap();
    let lines: Vec<&str> = deck.lines().collect();
    let mut equations = Vec::new();
    for (i, l) in lines.iter().enumerate() {
        if *l == "*EQUATION" {
            let count: usize = lines[i + 1].trim().parse().unwrap();
            let f: Vec<&str> = lines[i + 2].split(',').map(str::trim).collect();
            assert_eq!(f.len(), 3 * count, "{}", lines[i + 2]);
            let node: u32 = f[0].parse().unwrap();
            let terms = f.chunks(3).map(|t| (t[1].parse().unwrap(), t[2].parse().unwrap())).collect();
            assert!(f.chunks(3).all(|t| t[0] == f[0]));
            equations.push((node, terms));
        }
    }
    check(&equations);
    // Model data, before the first step; the plane is no `*BOUNDARY` line.
    let first_step = deck.find("*STEP").unwrap();
    assert!(deck.rfind("*EQUATION").unwrap() < first_step);
    assert!(!deck.contains("SYMM"), "{deck}");

    let e = abaqus::write(&plate(false)).unwrap_err();
    assert!(matches!(&e, Error::Unsupported { what, .. } if what.contains("boundary condition `sym`") && what.contains("some steps only")), "{e}");
}

#[test]
fn nastran_writes_one_mpc_set_per_step() {
    // A static deck and a modes deck (one solution per run).
    for (steps, sol) in [("Load", 101), ("Modes", 103)] {
        let m = fieldcraft_decks::with_steps(&plate(true), &[steps.to_string()]).unwrap();
        let deck = nastran::write(&m).unwrap();
        assert!(deck.contains(&format!("SOL {sol}")));
        let sid = deck.lines().find_map(|l| l.trim().strip_prefix("MPC = ")).unwrap().to_string();
        // MPC entries: SID, then (grid, component, coefficient) triples; the continuation after
        // two triples starts with a blank field. Large-field lines are 16 columns per field.
        let mut equations = Vec::new();
        let mut fields: Vec<String> = Vec::new();
        let flush = |fields: &mut Vec<String>, equations: &mut Vec<(u32, Vec<(u8, f64)>)>| {
            if fields.is_empty() {
                return;
            }
            assert_eq!(fields[0], sid);
            let data: Vec<&String> = fields[1..].iter().enumerate().filter(|(i, _)| !matches!(i % 8, 6 | 7)).map(|(_, f)| f).collect();
            let triples: Vec<(u32, u8, f64)> = data
                .chunks(3)
                .filter(|t| t.len() == 3 && !t[0].is_empty())
                .map(|t| (t[0].parse().unwrap(), t[1].parse().unwrap(), real(t[2])))
                .collect();
            assert!(triples.iter().all(|t| t.0 == triples[0].0));
            equations.push((triples[0].0, triples.iter().map(|t| (t.1, t.2)).collect()));
            fields.clear();
        };
        for line in deck.lines() {
            let (head, rest) = line.split_at(8.min(line.len()));
            let row: Vec<String> = rest.as_bytes().chunks(16).map(|c| String::from_utf8_lossy(c).trim().to_string()).collect();
            let mut row = row;
            row.resize(4, String::new());
            match head.trim() {
                "MPC*" => {
                    flush(&mut fields, &mut equations);
                    fields.extend(row);
                }
                "*" if !fields.is_empty() => fields.extend(row),
                _ => flush(&mut fields, &mut equations),
            }
        }
        flush(&mut fields, &mut equations);
        check(&equations);
    }
}

#[test]
fn lsdyna_holds_the_plane_in_a_local_system() {
    let deck = lsdyna::write(&plate(true)).unwrap().text;
    let block = |keyword: &str| -> Vec<Vec<f64>> {
        let at = deck.find(keyword).unwrap();
        deck[at..]
            .lines()
            .skip(1)
            .filter(|l| !l.starts_with('$'))
            .take_while(|l| !l.starts_with('*'))
            .map(|l| l.as_bytes().chunks(20).filter_map(|c| String::from_utf8_lossy(c).trim().parse().ok()).collect())
            .collect()
    };
    let system = block("*DEFINE_COORDINATE_SYSTEM");
    let (cid, x, p) = (system[0][0], [system[0][4], system[0][5], system[0][6]], [system[1][0], system[1][1], system[1][2]]);
    assert!((x[0] + x[1]).abs() < 1e-12 && x[2] == 0.0 && (dot(x, x) - 1.0).abs() < 1e-12, "{x:?}");
    assert!(dot(x, p).abs() < 1e-12, "{p:?}");
    // The SPC after the system: its id card, then nsid, cid and the held local directions.
    let at = deck.find("*DEFINE_COORDINATE_SYSTEM").unwrap();
    let spc = deck[at..].find("*BOUNDARY_SPC_SET_ID").unwrap() + at;
    let card: Vec<f64> = deck[spc..]
        .lines()
        .skip(2)
        .find(|l| !l.starts_with('$'))
        .unwrap()
        .as_bytes()
        .chunks(20)
        .map(|c| String::from_utf8_lossy(c).trim().parse().unwrap())
        .collect();
    assert_eq!(card[1], cid);
    assert_eq!(&card[2..], [1.0, 0.0, 0.0, 0.0, 1.0, 1.0]);
}
