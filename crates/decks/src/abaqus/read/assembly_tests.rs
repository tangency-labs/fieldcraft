use super::*;

const PART: &str = "*PART, NAME=P\n*NODE, NSET=All\n1,1,0,0\n2,2,0,0\n3,1,1,0\n4,1,0,1\n*ELEMENT, TYPE=C3D4, ELSET=Solid\n7,1,2,3,4\n*NSET,NSET=Tip,GENERATE\n2,4,1\n*SOLID SECTION, ELSET=Solid, MATERIAL=M\n*END PART\n";
const MAT: &str = "*MATERIAL,NAME=M\n*ELASTIC\n200000,0.3\n*DENSITY\n7.8e-9\n";

fn assembled(extra: &str) -> String {
    format!(
        "{PART}*ASSEMBLY,NAME=A\n*INSTANCE,NAME=First,PART=P\n1,0,0\n0,0,0,0,0,1,90\n*END INSTANCE\n*INSTANCE,NAME=Second,PART=P\n*END INSTANCE\n{extra}*END ASSEMBLY\n{MAT}"
    )
}

#[test]
fn instance_labels_sets_transform_and_qualified_references_are_resolved() {
    let input = assembled(
        "*NSET,NSET=Selected,INSTANCE=First,GENERATE\n1,4,1\n*NSET,NSET=FromElements,ELSET=Solid,INSTANCE=First\n*NSET,NSET=Across\nA.First.Tip, Second.1\n*ELSET,ELSET=Both\nFirst.7, A.Second.Solid\n*SURFACE,NAME=Face,TYPE=ELEMENT\nFirst.7,S1\n",
    );
    let got = read_in(&input, Preset::MmTonneS).unwrap();
    assert_eq!(got.model.nodes().len(), 8);
    assert_eq!(got.model.elements().len(), 2);
    let first = &got.instances[0];
    let second = &got.instances[1];
    let node = got.model.node(first.nodes[&1]).unwrap();
    // Translate (1,0,0) by (1,0,0), then rotate 90 degrees around global z.
    assert!(node.xyz[0].abs() < 1e-15);
    assert!((node.xyz[1] - 0.002).abs() < 1e-15);
    assert_eq!(got.model.node(second.nodes[&1]).unwrap().xyz, [0.001, 0.0, 0.0]);
    assert_eq!(got.model.set("SELECTED").unwrap().nodes, first.nodes.values().copied().collect::<Vec<_>>());
    assert_eq!(got.model.set("FROMELEMENTS").unwrap().nodes, first.nodes.values().copied().collect::<Vec<_>>());
    assert_eq!(got.model.set("ACROSS").unwrap().nodes, vec![first.nodes[&2], first.nodes[&3], first.nodes[&4], second.nodes[&1]]);
    assert_eq!(got.model.set("BOTH").unwrap().elements, vec![first.elements[&7], second.elements[&7]]);
    assert_eq!(got.model.set("FACE").unwrap().faces[0].element, first.elements[&7]);
    assert_ne!(first.names["SOLID"], second.names["SOLID"]);
    assert_eq!(read_in(&input, Preset::MmTonneS).unwrap().model, got.model);
}

#[test]
fn assembly_does_not_consume_analysis_history() {
    let history = "*STEP,NAME=Load\n*STATIC\n0.1,1\n*BOUNDARY\nFirst.All,ENCASTRE\n*CLOAD\nSecond.2,2,-10\n*END STEP\n";
    let input = format!("{}{history}", assembled(""));
    let got = read(&input).unwrap();
    assert!(got.model.steps().is_empty());
    let line = input.lines().position(|s| s.starts_with("*STEP")).unwrap() + 1;
    assert!(got.unknown.iter().any(|u| u.keyword == "STEP" && u.line == line));
    let doc = source::SourceDocument::parse(&input).unwrap();
    assert!(doc.original().ends_with(history));
}

#[test]
fn malformed_scopes_references_and_limits_are_named_without_panics() {
    for input in [
        "*END PART\n".to_string(),
        "*PART,NAME=P\n".to_string(),
        "*ASSEMBLY,NAME=A\n*INSTANCE,NAME=X,PART=Missing\n*END INSTANCE\n*END ASSEMBLY\n".to_string(),
        assembled("*NSET,NSET=Bad\nFirst.999\n"),
        assembled("*NSET,NSET=Bad,INSTANCE=Missing\n1\n"),
        assembled("*NSET,NSET=Bad,GENERATE\n1,9,0\n"),
        assembled("*NSET,NSET=Bad\nFirst.Missing\n"),
        assembled("").replace("NAME=Second", "NAME=First"),
        assembled("").replace("0,0,0,0,0,1,90", "0,0,0,0,0,0,90"),
        assembled("").replace("*END INSTANCE\n", ""),
    ] {
        let error = read(&input).unwrap_err().to_string();
        assert!(error.contains("line"), "{error}");
    }
    assert!(read_limited(&assembled(""), Preset::MKgS, Limits { nodes: 7, elements: 10, members: 100 }).is_err());
    assert!(read_limited(&assembled(""), Preset::MKgS, Limits { nodes: 100, elements: 1, members: 100 }).is_err());
    assert!(read_limited(&assembled(""), Preset::MKgS, Limits { nodes: 100, elements: 100, members: 1 }).is_err());
}

#[test]
fn mesh_in_instance_and_continued_connectivity_use_local_namespace() {
    let input = "*PART,NAME=P\n*END PART\n*ASSEMBLY,NAME=A\n*INSTANCE,NAME=I,PART=P\n*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n*ELEMENT,TYPE=C3D4\n9,1,2,\n3,4\n*END INSTANCE\n*END ASSEMBLY\n";
    let got = read(input).unwrap();
    assert_eq!(got.model.elements().len(), 1);
    assert_eq!(got.instances[0].elements[&9], 1);
}

#[test]
fn parse_errors_keep_original_physical_lines_after_expansion() {
    let input = assembled("").replace("200000,0.3", "BAD,0.3");
    let line = input.lines().position(|s| s == "BAD,0.3").unwrap() + 1;
    let error = read(&input).unwrap_err().to_string();
    assert!(error.contains(&format!("line {line}:")), "{error}");
}

#[test]
fn local_orientation_is_transformed_with_the_instance() {
    let input = assembled("").replace(
        "*SOLID SECTION, ELSET=Solid, MATERIAL=M",
        "*ORIENTATION,NAME=Axes\n1,0,0,0,1,0\n*SOLID SECTION, ELSET=Solid, MATERIAL=M, ORIENTATION=Axes",
    );
    let got = read(&input).unwrap();
    assert_eq!(got.model.sections().len(), 2);
    assert!(got.model.sections().values().all(|s| s.orientation.is_some()));
    assert_ne!(got.model.sections().values().next().unwrap().orientation, got.model.sections().values().last().unwrap().orientation);
}
