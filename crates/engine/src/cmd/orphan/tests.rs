use super::*;
use fieldcraft_model::units::Preset;
use fieldcraft_model::{DisplaySystem, ElementKind, Model, ReferencePoint, Set};

struct Input {
    dir: std::path::PathBuf,
    path: String,
}
impl Input {
    fn new(name: &str, bytes: &[u8]) -> Self {
        let dir = std::env::temp_dir().join(format!("fieldcraft-msh-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("input.msh").to_string_lossy().into_owned();
        std::fs::write(&path, bytes).unwrap();
        Self { dir, path }
    }
    fn write(&self, text: impl AsRef<[u8]>) {
        std::fs::write(&self.path, text).unwrap();
    }
}
impl Drop for Input {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
        let _ = std::fs::remove_dir(&self.dir);
    }
}

const TET: &str = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n$Elements\n1 1 50 50\n3 1 4 1\n50 1 2 3 4\n$EndElements\n";

/// Original small test writer, not a product exporter. Reorder the model's mids back onto
/// the edges in the manual's diagrams; this independently constructs Gmsh-order records.
fn writer(m: &Model, units: &DisplaySystem, boundaries: bool) -> String {
    let kinds: BTreeMap<_, _> = m.elements().iter().fold(BTreeMap::<ElementKind, Vec<u32>>::new(), |mut map, (&l, e)| {
        map.entry(e.kind).or_default().push(l);
        map
    });
    let mut blocks = Vec::<(u8, u32, u32, Vec<(u32, Vec<u32>)>)>::new();
    for (kind, labels) in kinds {
        let code = match kind {
            ElementKind::Tet4 => 4,
            ElementKind::Tet10 => 11,
            ElementKind::Wedge6 => 6,
            ElementKind::Hex8 => 5,
            ElementKind::Hex20 => 17,
            _ => panic!("test writer kind"),
        };
        let gmsh_edges: &[[usize; 2]] = match kind {
            ElementKind::Tet10 => &[[0, 1], [1, 2], [2, 0], [0, 3], [2, 3], [1, 3]],
            ElementKind::Hex20 => &[[0, 1], [0, 3], [0, 4], [1, 2], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 7], [5, 6], [6, 7]],
            _ => &[],
        };
        let records = labels
            .into_iter()
            .map(|label| {
                let e = m.element(label).unwrap();
                let mut nodes = e.nodes[..kind.corner_count()].to_vec();
                for &[a, b] in gmsh_edges {
                    let edge = kind.edges().iter().find(|edge| (edge[0] == a && edge[1] == b) || (edge[0] == b && edge[1] == a)).unwrap();
                    nodes.push(e.nodes[edge[2]]);
                }
                (label, nodes)
            })
            .collect();
        blocks.push((3, 1, code, records));
    }
    let mut next = m.elements().keys().next_back().copied().unwrap() + 1;
    if boundaries {
        for (i, name) in ["root", "tip"].iter().enumerate() {
            let resolved = m.resolve_set(name).unwrap();
            let mut records = Vec::new();
            let mut code = 0;
            for side in resolved.faces {
                let e = m.element(side.element).unwrap();
                let positions = e.kind.faces()[usize::from(side.face) - 1];
                code = match positions.len() {
                    3 => 2,
                    4 => 3,
                    6 => 9,
                    8 => 16,
                    _ => panic!("test face"),
                };
                records.push((next, positions.iter().map(|&p| e.nodes[p]).collect()));
                next += 1;
            }
            blocks.push((2, i as u32 + 1, code, records));
        }
    }
    let mut text = String::from("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");
    if boundaries {
        text.push_str("$PhysicalNames\n3\n3 1 \"all\"\n2 1 \"root\"\n2 2 \"tip\"\n$EndPhysicalNames\n$Entities\n0 0 2 1\n1 0 0 0 10 1 1 1 1 0\n2 0 0 0 10 1 1 1 2 0\n1 0 0 0 10 1 1 1 1 2 1 2\n$EndEntities\n");
    }
    text.push_str(&format!(
        "$Nodes\n1 {} {} {}\n3 1 0 {}\n",
        m.nodes().len(),
        m.nodes().keys().next().unwrap(),
        m.nodes().keys().next_back().unwrap(),
        m.nodes().len()
    ));
    for label in m.nodes().keys() {
        text.push_str(&format!("{label}\n"));
    }
    for node in m.nodes().values() {
        let p = node.xyz.map(|x| units.express(Quantity::Length, x).unwrap());
        text.push_str(&format!("{} {} {}\n", p[0], p[1], p[2]));
    }
    let labels: Vec<_> = blocks.iter().flat_map(|b| b.3.iter().map(|r| r.0)).collect();
    text.push_str(&format!(
        "$EndNodes\n$Elements\n{} {} {} {}\n",
        blocks.len(),
        labels.len(),
        labels.iter().min().unwrap(),
        labels.iter().max().unwrap()
    ));
    for (dim, entity, code, records) in blocks {
        text.push_str(&format!("{dim} {entity} {code} {}\n", records.len()));
        for (label, nodes) in records {
            text.push_str(&label.to_string());
            for n in nodes {
                text.push_str(&format!(" {n}"));
            }
            text.push('\n');
        }
    }
    text.push_str("$EndElements\n");
    text
}

#[test]
fn orphan_import_replaces_atomically_preserves_definitions_and_undoes() {
    let input = Input::new("replace", TET.as_bytes());
    let mut s = Session::default();
    s.execute("mesh.node.add", &json!({"label":99,"xyz":[2,3,4]})).unwrap();
    s.execute("model.material.create", &json!({"name":"Steel","elastic":{"type":"isotropic","young":"200 GPa","poisson":0.3}})).unwrap();
    let before = s.model.clone();
    let journal = s.journal.len();
    let out = s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap();
    assert_eq!(out["replaced"], json!({"nodes":1,"elements":0}));
    assert_eq!(out["coordinate_unit"], "mm");
    assert_eq!(out["elements"], json!({"tet4":1}));
    assert_eq!(s.model.without_mesh(), before.without_mesh());
    assert!(s.model.node(99).is_none());
    assert_eq!(s.journal.len(), journal + 1);
    assert_eq!(s.execute("mesh.quality", &json!({})).unwrap()["valid"], 1);
    s.execute("edit.undo", &json!({})).unwrap();
    assert_eq!(s.model, before);
    s.execute("edit.redo", &json!({})).unwrap();
    assert_eq!(s.model.elements().len(), 1);
}

#[test]
fn orphan_import_respects_si_and_imperial_coordinate_units() {
    let input = Input::new("units", TET.as_bytes());
    for preset in [Preset::MKgS, Preset::InLbfS] {
        let mut s = Session { default_units: DisplaySystem::new(preset), ..Session::default() };
        let unit = s.units().clone();
        let out = s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap();
        assert_eq!(out["coordinate_unit"], unit.symbol(Quantity::Length));
        assert_eq!(s.model.node(2).unwrap().xyz, [unit.to_si(Quantity::Length, 1.).unwrap(), 0., 0.]);
    }
}

#[test]
fn orphan_refusals_preserve_model_and_journal_even_when_labels_match() {
    let input = Input::new("refuse", TET.as_bytes());
    let mut s = Session::default();
    s.model.set_node(1, [2., 2., 2.]).unwrap();
    for text in [
        TET.replace("4.1 0 8", "2.2 0 8"),
        TET.replace("50 1 2 3 4", "50 1 2 3 99"),
        TET.replace("$EndElements", ""),
        TET.replace("3 1 4 1", "3 1 7 1"),
        TET.replace("0 0 0\n1 0 0", "NaN 0 0\n1 0 0"),
    ] {
        input.write(text);
        let before = (s.model.clone(), s.journal.clone());
        assert!(s.execute("mesh.import.orphan", &json!({"path":input.path})).is_err());
        assert_eq!((s.model.clone(), s.journal.clone()), before);
    }
    input.write(TET);
    s.model.add_set("old", Set { nodes: vec![1], ..Set::default() }).unwrap();
    let before = (s.model.clone(), s.journal.clone());
    let e = s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap_err().to_string();
    assert!(e.contains("label-based set `old`"), "{e}");
    assert_eq!((s.model.clone(), s.journal.clone()), before);
    s.model.remove_set("old").unwrap();
    s.model.add_set("old", Set::of_geometry(vec![fieldcraft_model::GeometryRef::body("Body").unwrap()])).unwrap();
    assert!(s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap_err().to_string().contains("geometry-based set"));
    s.model.remove_set("old").unwrap();
    s.model.add_reference_point("RP", ReferencePoint { xyz: [0.; 3], label: 2 }).unwrap();
    assert!(s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap_err().to_string().contains("reference point `RP`"));
    input.write(b"$MeshFormat\n4.1 1 8\n\xff");
    assert!(s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap_err().to_string().contains("binary MSH not supported yet"));
    input.write(TET);
    let mut capped = Session { limits: crate::Limits { max_import_nodes: 3, ..crate::Limits::default() }, ..Session::default() };
    assert!(capped.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap_err().to_string().contains("node count"));
}

fn beam(kind: &str, divisions: [usize; 3]) -> Session {
    let mut s = Session::default();
    s.execute("mesh.block",&json!({"body":"Beam","kind":kind,"corners":[[0,0,0],[10,0,0],[10,1,0],[0,1,0],[0,0,1],[10,0,1],[10,1,1],[0,1,1]],"divisions":divisions,"faces":{"6":"root","4":"tip"}})).unwrap();
    for (name, face) in [("all", None), ("root", Some("root")), ("tip", Some("tip"))] {
        let geometry = match face {
            Some(f) => json!([{"body":"Beam","face":f}]),
            None => json!([{"body":"Beam"}]),
        };
        s.execute("model.set.create", &json!({"name":name,"geometry":geometry})).unwrap();
    }
    s
}

#[test]
fn orphan_roundtrips_all_volume_kinds_and_physical_end_faces() {
    let input = Input::new("roundtrip", b"");
    for kind in ["tet4", "tet10", "wedge6", "hex8", "hex20"] {
        let original = beam(kind, [1, 1, 1]);
        input.write(writer(&original.model, original.units(), true));
        let mut s = Session::default();
        let out = s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap();
        assert_eq!(s.model.nodes(), original.model.nodes(), "{kind}");
        assert_eq!(s.model.elements(), original.model.elements(), "{kind}");
        assert!(s.model.classification().is_empty());
        for name in ["root", "tip"] {
            assert_eq!(s.model.resolve_set(name).unwrap(), original.model.resolve_set(name).unwrap(), "{kind} {name}");
        }
        assert_eq!(s.model.resolve_set("all").unwrap().elements, original.model.resolve_set("all").unwrap().elements);
        assert_eq!(out["sets"].as_array().unwrap().len(), 3);
        let quality = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(quality["valid"].as_u64().unwrap() as usize, s.model.elements().len(), "{kind}: {quality}");
        assert_eq!(quality["invalid"]["count"], 0);
    }
}

#[test]
fn orphan_hex20_physical_faces_solve_the_analytical_cantilever() {
    let original = beam("hex20", [20, 2, 2]);
    let input = Input::new("cantilever", writer(&original.model, original.units(), true).as_bytes());
    let mut s = Session::default();
    s.execute("mesh.import.orphan", &json!({"path":input.path})).unwrap();
    let recipe: Value = serde_json::from_str(include_str!("../../../../../tests/oracle/cantilever-hex20/recipe.json")).unwrap();
    for entry in recipe.as_array().unwrap().iter().skip(4) {
        s.execute(entry["command"].as_str().unwrap(), &entry["params"]).unwrap();
    }
    let out = s.execute("solve.probe", &json!({"step":"Bend","at":[[10.,0.5,0.5]]})).unwrap();
    let tip = out["nodes"][0]["displacement"][2].as_f64().unwrap();
    let expected: Value = serde_json::from_str(include_str!("../../../../../tests/oracle/cantilever-hex20/expected.json")).unwrap();
    let theory = expected["quantities"][0]["value"].as_f64().unwrap();
    assert!((tip - theory).abs() <= theory.abs() * 0.005, "{tip} vs {theory}");
    let accepted: Value = serde_json::from_str(include_str!("../../../../../tests/oracle/cantilever-hex20/fixture.json")).unwrap();
    let fixture = accepted["values"]["tip_uz"]["value"].as_f64().unwrap();
    assert!((tip - fixture).abs() < 1e-8, "{tip} vs accepted {fixture}");
}
