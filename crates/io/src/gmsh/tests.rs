//! Original handwritten fixtures, independently derived from the manual's reference elements.
use super::*;

type Block = (u8, u32, u32, Vec<(u32, Vec<u32>)>);

fn file(nodes: &[(u32, [f64; 3])], blocks: &[Block], prefix: &str) -> String {
    let tags: Vec<_> = nodes.iter().map(|n| n.0).collect();
    let mut text = format!(
        "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n{prefix}$Nodes\n1 {} {} {}\n3 1 0 {}\n",
        nodes.len(),
        tags.iter().min().unwrap(),
        tags.iter().max().unwrap(),
        nodes.len()
    );
    for (tag, _) in nodes {
        text.push_str(&format!("{tag}\n"));
    }
    for (_, [x, y, z]) in nodes {
        text.push_str(&format!("{x} {y} {z}\n"));
    }
    let labels: Vec<_> = blocks.iter().flat_map(|b| b.3.iter().map(|e| e.0)).collect();
    text.push_str(&format!(
        "$EndNodes\n$Elements\n{} {} {} {}\n",
        blocks.len(),
        labels.len(),
        labels.iter().min().unwrap(),
        labels.iter().max().unwrap()
    ));
    for (dim, entity, code, records) in blocks {
        text.push_str(&format!("{dim} {entity} {code} {}\n", records.len()));
        for (tag, nodes) in records {
            text.push_str(&tag.to_string());
            for node in nodes {
                text.push_str(&format!(" {node}"));
            }
            text.push('\n');
        }
    }
    text.push_str("$EndElements\n");
    text
}

fn read(text: &str) -> Imported {
    read_msh(text, Limits { nodes: 1000, elements: 1000, members: 10000 }, Unmatched::Refuse).unwrap()
}
fn tet() -> String {
    file(&[(10, [0., 0., 0.]), (20, [1., 0., 0.]), (30, [0., 1., 0.]), (40, [0., 0., 1.])], &[(3, 1, 4, vec![(77, vec![10, 20, 30, 40])])], "")
}

#[test]
fn sparse_labels_and_ignored_sections_survive() {
    let text = tet().replace("$Nodes", "$Comments\noriginal handwritten fixture\n$EndComments\n$Nodes");
    let out = read(&text);
    assert_eq!(out.model.nodes().keys().copied().collect::<Vec<_>>(), [10, 20, 30, 40]);
    assert_eq!(out.model.element(77).unwrap().nodes, [10, 20, 30, 40]);
    assert_eq!(out.skipped_sections, ["Comments"]);
}

#[test]
fn all_supported_kinds_preserve_the_reference_midpoint_edges() {
    let line = vec![[0., 0., 0.], [1., 0., 0.]];
    let tri = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    let quad = vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]];
    let tet = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let wedge = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]];
    let hex = vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    // Gmsh diagram edges, independent of the production permutation tables.
    let cases = [
        (1, ElementKind::Rod2, line.clone(), vec![]),
        (8, ElementKind::Rod3, line, vec![[0, 1]]),
        (2, ElementKind::Tri3, tri.clone(), vec![]),
        (9, ElementKind::Tri6, tri, vec![[0, 1], [1, 2], [2, 0]]),
        (3, ElementKind::Quad4, quad.clone(), vec![]),
        (16, ElementKind::Quad8, quad, vec![[0, 1], [1, 2], [2, 3], [3, 0]]),
        (4, ElementKind::Tet4, tet.clone(), vec![]),
        (11, ElementKind::Tet10, tet, vec![[0, 1], [1, 2], [2, 0], [0, 3], [2, 3], [1, 3]]),
        (6, ElementKind::Wedge6, wedge, vec![]),
        (5, ElementKind::Hex8, hex.clone(), vec![]),
        (17, ElementKind::Hex20, hex, vec![[0, 1], [0, 3], [0, 4], [1, 2], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 7], [5, 6], [6, 7]]),
    ];
    for (code, kind, mut xyz, mids) in cases {
        let corners = xyz.len();
        for [a, b] in mids {
            xyz.push(std::array::from_fn(|i| (xyz[a][i] + xyz[b][i]) / 2.));
        }
        let nodes: Vec<_> = xyz.iter().enumerate().map(|(i, &p)| (100 + 3 * i as u32, p)).collect();
        let out = read(&file(&nodes, &[(kind.dimension(), 1, code, vec![(900, nodes.iter().map(|n| n.0).collect())])], ""));
        let e = out.model.element(900).unwrap();
        assert_eq!(e.kind, kind);
        assert_eq!(&e.nodes[..corners], &nodes[..corners].iter().map(|n| n.0).collect::<Vec<_>>());
        for edge in kind.edges().iter().filter(|e| e.len() == 3) {
            let at = |i: usize| out.model.node(e.nodes[edge[i]]).unwrap().xyz;
            assert_eq!(at(2), std::array::from_fn(|i| (at(0)[i] + at(1)[i]) / 2.), "{kind:?} {edge:?}");
        }
    }
}

fn physical_tet() -> String {
    let prefix = "$PhysicalNames\n4\n3 1 \"solid\"\n2 2 \"wall\"\n1 3 \"rim\"\n0 4 \"vertex\"\n$EndPhysicalNames\n$Entities\n1 1 1 1\n1 0 0 0 1 4\n1 0 0 0 1 0 0 1 3 0\n1 0 0 0 1 1 0 1 2 0\n1 0 0 0 1 1 1 1 1 0\n$EndEntities\n";
    file(
        &[(1, [0., 0., 0.]), (2, [1., 0., 0.]), (3, [0., 1., 0.]), (4, [0., 0., 1.])],
        &[
            (3, 1, 4, vec![(10, vec![1, 2, 3, 4])]),
            (2, 1, 2, vec![(20, vec![1, 3, 2])]),
            (1, 1, 1, vec![(30, vec![1, 2])]),
            (0, 1, 15, vec![(40, vec![1])]),
        ],
        prefix,
    )
}

#[test]
fn physical_boundaries_select_sides_without_adding_stiffness() {
    let m = read(&physical_tet()).model;
    assert_eq!(m.elements().len(), 1);
    assert_eq!(m.set("solid").unwrap().elements, [10]);
    assert_eq!(m.set("wall").unwrap().faces, [ElementFace { element: 10, face: 1 }]);
    assert_eq!(m.set("rim").unwrap().edges, [ElementEdge { element: 10, edge: 1 }]);
    assert_eq!(m.set("vertex").unwrap().nodes, [1]);
    for name in ["solid", "wall", "rim", "vertex"] {
        assert!(m.resolve_set(name).is_ok());
    }
    let unnamed = physical_tet().replace("$PhysicalNames\n4\n3 1 \"solid\"\n2 2 \"wall\"\n1 3 \"rim\"\n0 4 \"vertex\"\n$EndPhysicalNames\n", "");
    assert!(read(&unnamed).model.set("physical_2_2").is_some());
}

#[test]
fn an_internal_physical_surface_selects_both_adjacent_faces() {
    let prefix =
        "$PhysicalNames\n1\n2 2 \"interface\"\n$EndPhysicalNames\n$Entities\n0 0 1 1\n1 0 0 0 1 1 0 1 2 0\n1 0 0 -1 1 1 1 0 0\n$EndEntities\n";
    let text = file(
        &[(1, [0., 0., 0.]), (2, [1., 0., 0.]), (3, [0., 1., 0.]), (4, [0., 0., 1.]), (5, [0., 0., -1.])],
        &[(3, 1, 4, vec![(10, vec![1, 2, 3, 4]), (11, vec![1, 3, 2, 5])]), (2, 1, 2, vec![(20, vec![1, 2, 3])])],
        prefix,
    );
    let m = read(&text).model;
    assert_eq!(m.elements().len(), 2);
    assert_eq!(m.set("interface").unwrap().faces, [ElementFace { element: 10, face: 1 }, ElementFace { element: 11, face: 1 }]);
    assert_eq!(m.resolve_set("interface").unwrap().nodes, [1, 2, 3]);
}

#[test]
fn late_entities_are_resolved_and_missing_entity_references_refuse() {
    let text = physical_tet();
    let start = text.find("$Entities\n").unwrap();
    let end = text.find("$EndEntities\n").unwrap() + "$EndEntities\n".len();
    let entities = &text[start..end];
    let late = text.replace(entities, "") + entities;
    assert_eq!(read(&late).model, read(&text).model);
    let missing = late.replace("3 1 0 4\n", "3 9 0 4\n");
    assert!(read_msh(&missing, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("missing entity (3, 9)"));
    let empty = tet().replace("$Nodes", "$Entities\n0 0 0 0\n$EndEntities\n$Nodes");
    assert!(read_msh(&empty, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("missing entity"));
}

#[test]
fn a_surface_mesh_keeps_its_shell_and_maps_physical_curves_to_edges() {
    let prefix = "$PhysicalNames\n3\n2 1 \"plate\"\n1 1 \"edge\"\n0 1 \"corner\"\n$EndPhysicalNames\n$Entities\n1 1 1 0\n1 0 0 0 1 1\n1 0 0 0 1 0 0 1 1 0\n1 0 0 0 1 1 0 1 1 0\n$EndEntities\n";
    let text = file(
        &[(1, [0., 0., 0.]), (2, [1., 0., 0.]), (3, [1., 1., 0.]), (4, [0., 1., 0.])],
        &[(2, 1, 3, vec![(10, vec![1, 2, 3, 4])]), (1, 1, 1, vec![(20, vec![1, 2])]), (0, 1, 15, vec![(30, vec![1])])],
        prefix,
    )
    .replace("3 1 0 4\n", "2 1 0 4\n");
    let m = read(&text).model;
    assert_eq!(m.elements().len(), 1);
    assert_eq!(m.element(10).unwrap().kind, ElementKind::Quad4);
    assert_eq!(m.set("plate").unwrap().elements, [10]);
    assert_eq!(m.set("edge").unwrap().edges, [ElementEdge { element: 10, edge: 1 }]);
    assert_eq!(m.set("corner").unwrap().nodes, [1]);
    let missing = text.replace("20 1 2\n", "20 1 3\n");
    assert!(read_msh(&missing, Limits::default(), Unmatched::Shell).unwrap_err().to_string().contains("matches no element edge"));
    let duplicate = tet().replace("1 1 77 77\n3 1 4 1\n77 10 20 30 40", "1 2 77 77\n3 1 4 2\n77 10 20 30 40\n77 10 20 30 40");
    assert!(read_msh(&duplicate, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("duplicate element tag 77"));
}

#[test]
fn parametric_coordinates_are_consumed() {
    let text = tet()
        .replace("3 1 0 4\n", "3 1 1 4\n")
        .replace("0 0 0\n1 0 0\n0 1 0\n0 0 1\n", "0 0 0 .1 .2 .3\n1 0 0 .4 .5 .6\n0 1 0 .7 .8 .9\n0 0 1 .1 .2 .3\n");
    assert_eq!(read(&text).model, read(&tet()).model);
}

#[test]
fn hostile_and_truncated_files_never_panic() {
    let good = physical_tet();
    let complete = good.trim_end();
    for end in 0..complete.len() {
        let result = std::panic::catch_unwind(|| read_msh(&complete[..end], Limits::default(), Unmatched::Refuse));
        assert!(result.is_ok_and(|r| r.is_err()), "byte {end}");
    }
    for (from, to, message) in [
        ("4.1 0 8", "2.2 0 8", "version 2.2"),
        ("4.1 0 8", "4.1 1 8", "binary MSH not supported yet"),
        ("1 4 10 40", "1 1000000000000000000 10 40", "node count"),
        ("1 4 10 40", "1 -4 10 40", "node count"),
        ("10\n20\n", "0\n20\n", "node tag"),
        ("10\n20\n", "4294967296\n20\n", "node tag"),
        ("0 0 0\n1 0 0", "NaN 0 0\n1 0 0", "coordinate"),
        ("0 0 0\n1 0 0", "1e400 0 0\n1 0 0", "coordinate"),
        ("10\n20\n", "20\n20\n", "duplicate node"),
        ("77 10 20 30 40", "77 10 20 30 99", "missing node 99"),
    ] {
        let e = read_msh(&tet().replace(from, to), Limits::default(), Unmatched::Refuse).unwrap_err().to_string();
        assert!(e.contains(message) && e.contains("line"), "{message}: {e}");
    }
    let text = tet().replace("$Nodes", &format!("$Comments\n{}\n$EndComments\n$Nodes", "x".repeat(10 << 20)));
    assert!(read_msh(&text, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("line exceeds"));
    assert!(read_msh(&tet(), Limits { nodes: 3, ..Limits::default() }, Unmatched::Refuse).is_err());
    assert!(read_msh(&good, Limits { members: 1, ..Limits::default() }, Unmatched::Refuse).is_err());
}

#[test]
fn unsupported_types_and_unsafe_boundaries_are_named() {
    for (code, name) in [(7, "pyramid"), (10, "quadrangle"), (12, "hexahedron"), (13, "prism"), (18, "prism"), (19, "pyramid"), (999, "unknown")] {
        let text = tet().replace("3 1 4 1\n", &format!("3 1 {code} 1\n"));
        let e = read_msh(&text, Limits::default(), Unmatched::Refuse).unwrap_err().to_string();
        assert!(e.contains(&format!("type {code}")) && e.contains(name), "{e}");
    }
    let text = tet().replace("$Nodes", "$PartitionedEntities\n0\n$EndPartitionedEntities\n$Nodes");
    assert!(read_msh(&text, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("partitioned"));
    for name in ["has space", "2bad", "bad/name"] {
        let text = physical_tet().replace("\"wall\"", &format!("\"{name}\""));
        assert!(read_msh(&text, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains(name));
    }
    let text = physical_tet().replace("2 1 2 1\n20 1 3 2", "2 1 3 1\n20 1 2 3 4");
    assert!(read_msh(&text, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("matches no volume face"));
    let out = read_msh(&text, Limits::default(), Unmatched::Shell).unwrap();
    assert_eq!(out.model.element(20).unwrap().kind, ElementKind::Quad4);
    assert_eq!(out.model.set("wall").unwrap().elements, [20]);
    let text = file(&[(1, [0., 0., 0.]), (2, [1., 0., 0.])], &[(1, 1, 1, vec![(1, vec![1, 2])]), (0, 1, 15, vec![(2, vec![1])])], "");
    assert!(read_msh(&text, Limits::default(), Unmatched::Refuse).unwrap_err().to_string().contains("0-dimensional physical group"));
}
