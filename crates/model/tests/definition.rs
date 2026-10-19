//! A whole model definition on a tagged mesh: sets by geometry name survive a remesh, the check
//! names what is missing, and the model round-trips through JSON.

use fieldcraft_model::*;

/// A bar of `n` unit hexahedra along X, labels starting after `offset`, tagged as body `Bar` with
/// faces `left` (x = 0) and `right` (x = n), and the vertex `tip` at (n, 0, 0).
fn bar(m: &mut Model, n: u32, offset: u32) {
    m.clear_geometry_tags();
    for e in m.elements().keys().copied().collect::<Vec<_>>() {
        m.remove_element(e).unwrap();
    }
    for v in m.nodes().keys().copied().collect::<Vec<_>>() {
        m.remove_node(v).unwrap();
    }
    let node = |i: u32, j: u32, k: u32| offset + 1 + i * 4 + j * 2 + k;
    for i in 0..=n {
        for j in 0..2 {
            for k in 0..2 {
                m.set_node(node(i, j, k), [f64::from(i), f64::from(j), f64::from(k)]).unwrap();
            }
        }
    }
    let body = m.geometry_tag(&GeometryRef::body("Bar").unwrap()).unwrap();
    let left = m.geometry_tag(&GeometryRef::face("Bar", "left").unwrap()).unwrap();
    let right = m.geometry_tag(&GeometryRef::face("Bar", "right").unwrap()).unwrap();
    let tip = m.geometry_tag(&GeometryRef::new("Bar", Entity::Vertex("tip".into())).unwrap()).unwrap();
    for i in 0..n {
        let label = offset + 1 + i;
        let nodes = vec![
            node(i, 0, 0),
            node(i + 1, 0, 0),
            node(i + 1, 1, 0),
            node(i, 1, 0),
            node(i, 0, 1),
            node(i + 1, 0, 1),
            node(i + 1, 1, 1),
            node(i, 1, 1),
        ];
        m.set_element(label, ElementKind::Hex8, nodes).unwrap();
        m.tag_element(label, body).unwrap();
        if i == 0 {
            m.tag_element_face(label, 6, left).unwrap();
        }
        if i + 1 == n {
            m.tag_element_face(label, 4, right).unwrap();
        }
    }
    m.tag_node(node(n, 0, 0), tip).unwrap();
}

fn defined() -> Model {
    let mut m = Model::default();
    bar(&mut m, 2, 0);
    m.add_set("all", Set::of_geometry(vec![GeometryRef::body("Bar").unwrap()])).unwrap();
    m.add_set("root", Set::of_geometry(vec![GeometryRef::face("Bar", "left").unwrap()])).unwrap();
    m.add_set("end", Set::of_geometry(vec![GeometryRef::face("Bar", "right").unwrap()])).unwrap();
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 2.1e11, poisson: 0.3 }), density: Some(7850.0), ..Material::default() },
    )
    .unwrap();
    m.add_section("bar", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.add_reference_point("RP", ReferencePoint { xyz: [3.0, 0.5, 0.5], label: 0 }).unwrap();
    m.add_coupling("load-pad", Coupling { reference_point: "RP".into(), set: "end".into(), kind: CouplingKind::Distributing, dofs: Dofs::all() })
        .unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    m.add_load("push", Load::Force { target: Target::ReferencePoint("RP".into()), vector: [0.0, 0.0, -100.0], distribution: Distribution::Total })
        .unwrap();
    m.add_load("press", Load::Pressure { set: "end".into(), magnitude: 1.0 }).unwrap();
    let mut s = Step::new("Bend", Procedure::Static {});
    s.boundary_conditions.insert("clamp".into());
    s.loads.extend(["push".to_string(), "press".to_string()]);
    m.add_step(s, None).unwrap();
    let mut f = Step::new("Modes", Procedure::Frequency { modes: 6, min_hz: None, max_hz: None });
    f.boundary_conditions.insert("clamp".into());
    m.add_step(f, None).unwrap();
    m
}

#[test]
fn a_complete_model_resolves_and_survives_a_remesh() {
    let mut m = defined();
    assert_eq!(m.check(), vec![]);
    let r = m.resolve().unwrap();
    assert_eq!(r.sets["all"].elements, vec![1, 2]);
    assert_eq!(r.sets["all"].nodes.len(), 12);
    assert_eq!(r.sets["all"].faces, vec![ElementFace { element: 1, face: 6 }, ElementFace { element: 2, face: 4 }], "a body's skin");
    assert_eq!(r.sets["root"].faces, vec![ElementFace { element: 1, face: 6 }]);
    assert_eq!(r.sets["root"].nodes, vec![1, 2, 3, 4]);
    assert!(r.sets["root"].elements.is_empty(), "a face holds no solid elements");
    assert_eq!(r.sets["end"].nodes, vec![9, 10, 11, 12]);
    assert_eq!(r.reference_point_nodes["RP"], MAX_LABEL, "from the top of the label range");
    assert_eq!(r.nodal_loads["push"], vec![(MAX_LABEL, [0.0, 0.0, -100.0])]);
    assert_eq!(r.target_nodes(&Target::ReferencePoint("RP".into())), [MAX_LABEL]);

    // The pressure acts on the end face's outward-oriented nodes (x = 2).
    let face = m.face_nodes(ElementFace { element: 2, face: 4 }).unwrap();
    assert!(face.iter().all(|n| m.node(*n).unwrap().xyz[0] == 2.0));

    // Remesh with four elements and other labels: the same definitions resolve on the new mesh.
    bar(&mut m, 4, 100);
    assert_eq!(m.check(), vec![]);
    let r = m.resolve().unwrap();
    assert_eq!(r.sets["all"].elements, vec![101, 102, 103, 104]);
    assert_eq!(r.sets["end"].faces, vec![ElementFace { element: 104, face: 4 }]);
    assert_eq!(r.sets["end"].nodes, vec![117, 118, 119, 120]);
    assert_eq!(r.reference_point_nodes["RP"], MAX_LABEL, "a remesh does not renumber reference points");

    // Vertex members resolve to their node.
    m.add_set("tip", Set::of_geometry(vec![GeometryRef::new("Bar", Entity::Vertex("tip".into())).unwrap()])).unwrap();
    assert_eq!(m.resolve_set("tip").unwrap().nodes, vec![117]);
}

#[test]
fn the_check_names_what_is_missing() {
    let mut m = defined();
    // A set naming a face the mesh does not carry, and explicit labels not in the mesh.
    m.add_set("ghost", Set { geometry: vec![GeometryRef::face("Bar", "top").unwrap()], nodes: vec![999], ..Set::default() }).unwrap();
    let e = m.resolve_set("ghost").unwrap_err().to_string();
    assert!(e.contains("face top of body Bar is not on the mesh") && e.contains("node 999"), "{e}");
    m.add_set("side", Set { faces: vec![ElementFace { element: 1, face: 3 }], ..Set::default() }).unwrap();
    m.add_load("ghost-load", Load::Pressure { set: "ghost".into(), magnitude: 1.0 }).unwrap();
    m.add_load("twist", Load::Moment { target: Target::Set("end".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    m.add_material("Soft", Material { elastic: Some(Elastic::Isotropic { young: 2.1e17, poisson: 0.3 }), density: None, ..Material::default() })
        .unwrap();
    m.add_section("again", Section { set: "side".into(), material: "Soft".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.add_section(
        "shell",
        Section {
            set: "root".into(),
            material: "Steel".into(),
            kind: SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();

    let diagnostics = m.check();
    let text: Vec<String> = diagnostics.iter().map(ToString::to_string).collect();
    let has = |needle: &str| text.iter().any(|t| t.contains(needle));
    assert!(has("error: set `ghost`: face top of body Bar is not on the mesh"), "{text:#?}");
    assert!(has("error: set `ghost`: node 999 is not in the mesh"), "{text:#?}");
    assert!(has("error: load `ghost-load`: set `ghost` has no element faces"), "{text:#?}");
    assert!(has("error: load `twist`: set `end` has no nodes with rotations"), "{text:#?}");
    assert!(has("error: section `shell`: set `root` has no elements"), "{text:#?}");
    assert!(has("error: section `again`: set `side` has no elements"), "{text:#?}");
    assert!(has("error: material `Soft`: has no density, and step `Modes` extracts frequencies"), "{text:#?}");
    assert!(has("warning: material `Soft`: Young's modulus 2.1e11 MPa"), "{text:#?}");
    assert!(has("warning: load `twist`: is not active in any step"), "{text:#?}");
    assert_eq!(diagnostics.first().map(|d| d.severity), Some(Severity::Error), "errors first");
    assert!(m.resolve().is_err());

    // Two sections over the same elements, and an element with none.
    let mut m = defined();
    m.add_section("twice", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    m.set_node(50, [9.0, 9.0, 9.0]).unwrap();
    m.set_node(51, [10.0, 9.0, 9.0]).unwrap();
    m.set_element(9, ElementKind::Rod2, vec![50, 51]).unwrap();
    let text: Vec<String> = m.check().iter().map(ToString::to_string).collect();
    assert!(text.iter().any(|t| t.contains("section `twice`: elements 1, 2 are also in section `bar`")), "{text:#?}");
    assert!(text.iter().any(|t| t.contains("model: elements 9 have no section")), "{text:#?}");

    // Incompatible modes on tetrahedra, a beam oriented along its axis.
    let mut m = defined();
    m.replace_section(
        "bar",
        Section {
            set: "all".into(),
            material: "Steel".into(),
            kind: SectionKind::Solid {},
            formulation: Formulation { incompatible_modes: true, ..Formulation::default() },
        },
    )
    .unwrap();
    assert_eq!(m.check().iter().filter(|d| d.severity == Severity::Error).count(), 0, "hex8 takes incompatible modes");
    m.set_element(2, ElementKind::Tet4, vec![5, 6, 7, 9]).unwrap();
    let body = m.geometry_tag(&GeometryRef::body("Bar").unwrap()).unwrap();
    m.tag_element(2, body).unwrap();
    assert!(m.check().iter().any(|d| d.message.contains("incompatible modes apply to hex8 elements only, not elements 2")));
    m.set_element(30, ElementKind::Beam2, vec![1, 5]).unwrap();
    m.add_set("beam", Set { elements: vec![30], ..Set::default() }).unwrap();
    let profile = BeamProfile::Rectangle { width: 1.0, height: 2.0 };
    m.add_section(
        "beam",
        Section {
            set: "beam".into(),
            material: "Steel".into(),
            kind: SectionKind::Beam { profile, orientation: [-3.0, 0.0, 0.0] },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    assert!(m.check().iter().any(|d| d.subject == "section `beam`" && d.message.contains("parallel to elements 30")));
}

#[test]
fn definitions_refer_only_to_existing_items_and_cannot_be_deleted_while_used() {
    let mut m = defined();
    let section = |set: &str, material: &str| Section {
        set: set.into(),
        material: material.into(),
        kind: SectionKind::Solid {},
        formulation: Formulation::default(),
    };
    assert!(matches!(m.add_section("x", section("nope", "Steel")), Err(Error::NotFound(_))));
    assert!(matches!(m.add_section("x", section("all", "Nope")), Err(Error::NotFound(_))));
    assert!(
        m.add_load("x", Load::Force { target: Target::ReferencePoint("nope".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total })
            .is_err()
    );
    assert!(
        m.add_coupling("x", Coupling { reference_point: "RP".into(), set: "end".into(), kind: CouplingKind::Kinematic, dofs: Dofs::none() }).is_err()
    );

    let e = m.remove_set("end").unwrap_err().to_string();
    assert!(e.contains("coupling `load-pad`") && e.contains("load `press`"), "{e}");
    assert!(m.remove_material("Steel").unwrap_err().to_string().contains("section `bar`"));
    assert!(m.remove_reference_point("RP").unwrap_err().to_string().contains("load `push`"));
    assert!(m.remove_load("push").unwrap_err().to_string().contains("step `Bend`"));
    m.set_load_active("Bend", "push", false).unwrap();
    m.remove_load("push").unwrap();
    m.remove_coupling("load-pad").unwrap();
    m.remove_reference_point("RP").unwrap();

    // Removing a tagged element drops its tags; a new element under the label starts untagged.
    m.remove_element(1).unwrap();
    assert!(!m.classification().elements().contains_key(&1));
}

#[test]
fn the_model_round_trips_through_json() {
    let m = defined();
    let json = serde_json::to_string(&m).unwrap();
    let back: Model = serde_json::from_str(&json).unwrap();
    assert_eq!(back, m);
    assert_eq!(serde_json::to_string(&back).unwrap(), json, "deterministic");
    // Reading checks what commands check: a dangling reference, a bad tag, a typo.
    let mut v = serde_json::to_value(&m).unwrap();
    v["sections"]["bar"]["material"] = serde_json::json!("Nope");
    assert!(serde_json::from_value::<Model>(v).unwrap_err().to_string().contains("material `Nope`"));
    let mut v = serde_json::to_value(&m).unwrap();
    v["classification"]["elements"]["1"]["faces"]["6"] = serde_json::json!(0);
    assert!(serde_json::from_value::<Model>(v).is_err(), "tag 0 is the body, not a face");
    let mut v = serde_json::to_value(&m).unwrap();
    v["materials"]["Steel"]["elastic"]["young"] = serde_json::json!(-1.0);
    assert!(serde_json::from_value::<Model>(v).is_err());
    let mut v = serde_json::to_value(&m).unwrap();
    v["stepz"] = serde_json::json!([]);
    assert!(serde_json::from_value::<Model>(v).is_err());

    // A model file from before the definitions existed still reads.
    let old: Model = serde_json::from_value(serde_json::json!({"nodes": {}, "elements": {}})).unwrap();
    assert_eq!(old, Model::default());
}

#[test]
fn tags_are_checked() {
    let mut m = defined();
    let face = m.geometry_tag(&GeometryRef::face("Bar", "left").unwrap()).unwrap();
    let body = m.geometry_tag(&GeometryRef::body("Bar").unwrap()).unwrap();
    assert!(m.tag_element(1, face).is_err(), "a solid fills a body, not a face");
    assert!(m.tag_element_face(1, 7, face).is_err());
    assert!(m.tag_element_face(1, 0, face).is_err());
    assert!(m.tag_element_face(1, 1, body).is_err());
    assert!(m.tag_element_edge(1, 13, 0).is_err());
    assert!(m.tag_node(999, 0).is_err());
    assert!(m.tag_element(1, 9999).is_err());
    assert!(m.tag_element(77, body).is_err());
    assert_eq!(m.geometry_tag(&GeometryRef::face("Bar", "left").unwrap()).unwrap(), face, "interned once");
}

/// A block `[0, 1] × [0, 1] × [0, 1]` of `n × n × n` hex8 with graded spacing (so area weighting
/// matters), tagged as body `Block` with its face `end` at x = 1.
fn block(m: &mut Model, n: u32) {
    m.clear_geometry_tags();
    for e in m.elements().keys().copied().collect::<Vec<_>>() {
        m.remove_element(e).unwrap();
    }
    for v in m.nodes().keys().copied().collect::<Vec<_>>() {
        m.remove_node(v).unwrap();
    }
    let at = |i: u32| (f64::from(i) / f64::from(n)).powi(2);
    let node = |i: u32, j: u32, k: u32| 1 + i + j * (n + 1) + k * (n + 1) * (n + 1);
    for k in 0..=n {
        for j in 0..=n {
            for i in 0..=n {
                m.set_node(node(i, j, k), [at(i), at(j), at(k)]).unwrap();
            }
        }
    }
    let body = m.geometry_tag(&GeometryRef::body("Block").unwrap()).unwrap();
    let end = m.geometry_tag(&GeometryRef::face("Block", "end").unwrap()).unwrap();
    let mut label = 0;
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                label += 1;
                let nodes = vec![
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j + 1, k),
                    node(i, j + 1, k),
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i + 1, j + 1, k + 1),
                    node(i, j + 1, k + 1),
                ];
                m.set_element(label, ElementKind::Hex8, nodes).unwrap();
                m.tag_element(label, body).unwrap();
                if i + 1 == n {
                    m.tag_element_face(label, 4, end).unwrap();
                }
            }
        }
    }
}

fn sum(values: &[(u32, [f64; 3])]) -> [f64; 3] {
    values.iter().fold([0.0; 3], |s, (_, v)| [s[0] + v[0], s[1] + v[1], s[2] + v[2]])
}

#[test]
fn a_total_force_on_a_face_is_the_same_on_every_mesh() {
    let mut m = Model::default();
    block(&mut m, 2);
    m.add_set("all", Set::of_geometry(vec![GeometryRef::body("Block").unwrap()])).unwrap();
    m.add_set("end", Set::of_geometry(vec![GeometryRef::face("Block", "end").unwrap()])).unwrap();
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 2.1e11, poisson: 0.3 }), density: Some(7850.0), ..Material::default() },
    )
    .unwrap();
    m.add_section("block", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    let force = |distribution| Load::Force { target: Target::Set("end".into()), vector: [100.0, 0.0, -50.0], distribution };
    m.add_load("total", force(Distribution::Total)).unwrap();
    m.add_load("each", force(Distribution::PerNode)).unwrap();
    m.add_load("drag", Load::Traction { set: "end".into(), vector: [0.0, 2.0, 0.0] }).unwrap();
    let mut s = Step::new("Pull", Procedure::Static {});
    s.loads.extend(["total".to_string(), "each".to_string(), "drag".to_string()]);
    m.add_step(s, None).unwrap();

    let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9);
    let coarse = m.resolve().unwrap();
    block(&mut m, 5);
    let fine = m.resolve().unwrap();
    for r in [&coarse, &fine] {
        assert!(close(sum(&r.nodal_loads["total"]), [100.0, 0.0, -50.0]), "{:?}", sum(&r.nodal_loads["total"]));
        // Traction × area (the end face is 1 × 1).
        assert!(close(sum(&r.nodal_loads["drag"]), [0.0, 2.0, 0.0]));
        assert_eq!(r.nodal_loads["total"].len(), r.sets["end"].nodes.len());
    }
    // A per-node force multiplies with the node count, by request.
    assert!(close(sum(&coarse.nodal_loads["each"]), [900.0, 0.0, -450.0]));
    assert!(close(sum(&fine.nodal_loads["each"]), [3600.0, 0.0, -1800.0]));
    let weights: Vec<f64> = fine.nodal_loads["total"].iter().map(|(_, v)| v[0]).collect();
    assert!(weights.iter().any(|w| (w - weights[0]).abs() > 1e-6), "graded faces give unequal shares");

    // A node set without faces splits a total force equally.
    m.add_set("pins", Set { nodes: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    m.add_load("pinned", Load::Force { target: Target::Set("pins".into()), vector: [4.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    assert_eq!(
        m.resolve().unwrap().nodal_loads["pinned"],
        vec![(1, [1.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]), (3, [1.0, 0.0, 0.0]), (4, [1.0, 0.0, 0.0])]
    );

    // A traction needs faces.
    m.add_load("nowhere", Load::Traction { set: "pins".into(), vector: [1.0, 0.0, 0.0] }).unwrap();
    assert!(m.check().iter().any(|d| d.subject == "load `nowhere`" && d.message.contains("no element faces")));
}

#[test]
fn reference_point_labels_are_stable() {
    let mut m = defined();
    m.add_reference_point("A", ReferencePoint { xyz: [0.0; 3], label: 0 }).unwrap();
    m.add_reference_point("Z", ReferencePoint { xyz: [0.0; 3], label: 0 }).unwrap();
    let before = m.reference_point_nodes();
    assert_eq!((before["RP"], before["A"], before["Z"]), (MAX_LABEL, MAX_LABEL - 1, MAX_LABEL - 2), "in creation order");
    // Another point, more nodes, a remesh and an edit renumber nothing.
    m.add_reference_point("B", ReferencePoint { xyz: [1.0; 3], label: 0 }).unwrap();
    m.set_node(5000, [9.0; 3]).unwrap();
    bar(&mut m, 3, 40);
    m.replace_reference_point("A", ReferencePoint { xyz: [2.0; 3], label: 0 }).unwrap();
    let after = m.reference_point_nodes();
    for name in ["RP", "A", "Z"] {
        assert_eq!(after[name], before[name], "{name}");
    }
    // Removing one frees its label for the next; explicit labels are checked.
    m.remove_reference_point("Z").unwrap();
    m.add_reference_point("C", ReferencePoint { xyz: [0.0; 3], label: 0 }).unwrap();
    assert_eq!(m.reference_point_nodes()["C"], MAX_LABEL - 2);
    assert!(m.add_reference_point("D", ReferencePoint { xyz: [0.0; 3], label: MAX_LABEL }).is_err(), "taken");
    assert!(m.add_reference_point("D", ReferencePoint { xyz: [0.0; 3], label: MAX_LABEL + 1 }).is_err());
    m.add_reference_point("D", ReferencePoint { xyz: [0.0; 3], label: 41 }).unwrap();
    assert!(
        m.check().iter().any(|d| d.subject == "reference point `D`" && d.message.contains("also a mesh node")),
        "41 is a node of the remeshed bar"
    );
}

/// Two quad4 in the XY plane (x from `x0`, 2 × 1), tagged as face `area` of body `Plate` with the
/// edges `left` (x = x0) and `right`, and a static step.
fn plate(kind: SectionKind, x0: f64) -> Model {
    let mut m = Model::default();
    for (i, x) in [x0, x0 + 1.0, x0 + 2.0].iter().enumerate() {
        m.set_node(i as u32 + 1, [*x, 0.0, 0.0]).unwrap();
        m.set_node(i as u32 + 4, [*x, 1.0, 0.0]).unwrap();
    }
    m.set_element(1, ElementKind::Quad4, vec![1, 2, 5, 4]).unwrap();
    m.set_element(2, ElementKind::Quad4, vec![2, 3, 6, 5]).unwrap();
    let area = m.geometry_tag(&GeometryRef::face("Plate", "area").unwrap()).unwrap();
    let left = m.geometry_tag(&GeometryRef::new("Plate", Entity::Edge("left".into())).unwrap()).unwrap();
    let right = m.geometry_tag(&GeometryRef::new("Plate", Entity::Edge("right".into())).unwrap()).unwrap();
    m.tag_element(1, area).unwrap();
    m.tag_element(2, area).unwrap();
    m.tag_element_edge(1, 4, left).unwrap();
    m.tag_element_edge(2, 2, right).unwrap();
    m.add_set("all", Set::of_geometry(vec![GeometryRef::face("Plate", "area").unwrap()])).unwrap();
    m.add_set("left", Set::of_geometry(vec![GeometryRef::new("Plate", Entity::Edge("left".into())).unwrap()])).unwrap();
    m.add_set("right", Set::of_geometry(vec![GeometryRef::new("Plate", Entity::Edge("right".into())).unwrap()])).unwrap();
    m.add_material(
        "Steel",
        Material { elastic: Some(Elastic::Isotropic { young: 2.1e11, poisson: 0.3 }), density: Some(7850.0), ..Material::default() },
    )
    .unwrap();
    m.add_section("plate", Section { set: "all".into(), material: "Steel".into(), kind, formulation: Formulation::default() }).unwrap();
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("left".into()) }).unwrap();
    m.add_load("pull", Load::Traction { set: "right".into(), vector: [10.0, 0.0, 0.0] }).unwrap();
    m.add_load("total", Load::Force { target: Target::Set("right".into()), vector: [0.0, 6.0, 0.0], distribution: Distribution::Total }).unwrap();
    let mut s = Step::new("Pull", Procedure::Static {});
    s.boundary_conditions.insert("clamp".into());
    s.loads.extend(["pull".to_string(), "total".to_string()]);
    m.add_step(s, None).unwrap();
    m
}

#[test]
fn two_dimensional_models_resolve_edges_and_check_their_rules() {
    // Edge sets resolve to element edges; tractions spread over length × thickness.
    let m = plate(SectionKind::PlaneStress { thickness: 2.0 }, 0.0);
    let r = m.resolve().unwrap();
    assert!(r.is_planar());
    assert_eq!(r.planar_edges("right"), vec![ElementEdge { element: 2, edge: 2 }]);
    assert_eq!(r.sets["left"].nodes, vec![1, 4]);
    let pull: f64 = r.nodal_loads["pull"].iter().map(|(_, v)| v[0]).sum();
    assert!((pull - 10.0 * 1.0 * 2.0).abs() < 1e-12, "{pull}");
    assert_eq!(r.nodal_loads["total"], vec![(3, [0.0, 3.0, 0.0]), (6, [0.0, 3.0, 0.0])]);
    // Axisymmetric: the right edge at r = 2 is a cylinder of area 2π·2·1.
    let r = plate(SectionKind::Axisymmetric {}, 0.0).resolve().unwrap();
    let pull: f64 = r.nodal_loads["pull"].iter().map(|(_, v)| v[0]).sum();
    assert!((pull - 10.0 * 4.0 * std::f64::consts::PI).abs() < 1e-9, "{pull}");

    let errors = |m: &Model| -> Vec<String> { m.check().into_iter().filter(|d| d.severity == Severity::Error).map(|d| d.to_string()).collect() };
    assert!(errors(&plate(SectionKind::Axisymmetric {}, -1.0)).iter().any(|e| e.contains("negative x")));
    let mut m = plate(SectionKind::PlaneStrain { thickness: 1.0 }, 0.0);
    m.set_node(6, [2.0, 1.0, 0.5]).unwrap();
    m.set_element(1, ElementKind::Quad4, vec![1, 4, 5, 2]).unwrap();
    let area = m.geometry_tag(&GeometryRef::face("Plate", "area").unwrap()).unwrap();
    m.tag_element(1, area).unwrap();
    m.add_load("lift", Load::Force { target: Target::Set("right".into()), vector: [0.0, 0.0, 1.0], distribution: Distribution::Total }).unwrap();
    m.add_load("spin", Load::Centrifugal { angular_velocity: 10.0, point: [0.0; 3], axis: [0.0, 1.0, 0.0], set: None }).unwrap();
    let e = errors(&m);
    assert!(e.iter().any(|e| e.contains("off the XY plane")), "{e:#?}");
    assert!(e.iter().any(|e| e.contains("clockwise")), "{e:#?}");
    assert!(e.iter().any(|e| e.contains("load `lift`") && e.contains("Z component")), "{e:#?}");
    assert!(e.iter().any(|e| e.contains("load `spin`") && e.contains("about an axis along Z")), "{e:#?}");
    // A 2D section with a solid one; quad4 takes incompatible modes in a 2D section.
    let mut m = plate(SectionKind::PlaneStress { thickness: 1.0 }, 0.0);
    m.replace_section(
        "plate",
        Section {
            set: "all".into(),
            material: "Steel".into(),
            kind: SectionKind::PlaneStress { thickness: 1.0 },
            formulation: Formulation { incompatible_modes: true, ..Formulation::default() },
        },
    )
    .unwrap();
    assert!(errors(&m).is_empty(), "{:#?}", errors(&m));
    for v in 1..=4 {
        m.set_node(10 + v, [f64::from(v & 1), f64::from(v >> 1 & 1), 5.0]).unwrap();
    }
    m.set_element(9, ElementKind::Tet4, vec![11, 12, 13, 14]).unwrap();
    m.add_set("tet", Set { elements: vec![9], ..Set::default() }).unwrap();
    m.add_section("tet", Section { set: "tet".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    assert!(errors(&m).iter().any(|e| e.contains("cannot be mixed")), "{:#?}", errors(&m));
    // Pressure on an edge set of a 2D model is fine; on nothing it is an error.
    let mut m = plate(SectionKind::PlaneStress { thickness: 1.0 }, 0.0);
    m.add_load("press", Load::Pressure { set: "right".into(), magnitude: 1.0 }).unwrap();
    m.add_load("bad", Load::Pressure { set: "all".into(), magnitude: 1.0 }).unwrap();
    let e = errors(&m);
    assert!(!e.iter().any(|e| e.contains("load `press`")) && e.iter().any(|e| e.contains("load `bad`")), "{e:#?}");
    // Explicit element edges round-trip through JSON; edge 0 is an error.
    let set = Set { edges: vec![ElementEdge { element: 2, edge: 2 }], ..Set::default() };
    let json = serde_json::to_value(&set).unwrap();
    assert_eq!(json, serde_json::json!({"edges": [{"element": 2, "edge": 2}]}));
    assert!(m.add_set("zero", Set { edges: vec![ElementEdge { element: 2, edge: 0 }], ..Set::default() }).is_err());
    assert!(m.add_set("five", Set { edges: vec![ElementEdge { element: 2, edge: 5 }], ..Set::default() }).is_ok());
    assert!(m.resolve_set("five").is_err());
    let json = serde_json::to_string(&m).unwrap();
    assert_eq!(serde_json::from_str::<Model>(&json).unwrap(), m);
    // Hostile centrifugal loads.
    for bad in [
        Load::Centrifugal { angular_velocity: 0.0, point: [0.0; 3], axis: [0.0, 1.0, 0.0], set: None },
        Load::Centrifugal { angular_velocity: 1e200, point: [0.0; 3], axis: [0.0, 1.0, 0.0], set: None },
        Load::Centrifugal { angular_velocity: 1.0, point: [f64::NAN, 0.0, 0.0], axis: [0.0, 1.0, 0.0], set: None },
        Load::Centrifugal { angular_velocity: 1.0, point: [0.0; 3], axis: [0.0; 3], set: None },
    ] {
        assert!(m.add_load("hostile", bad.clone()).is_err(), "{bad:?}");
    }
}

#[test]
fn display_resolves_what_is_on_the_mesh() {
    let mut m = defined();
    // A set with a node the mesh does not have: the check refuses it, the display shows the rest.
    m.add_set("half", Set { nodes: vec![1, 999], ..Set::default() }).unwrap();
    assert!(m.resolve_set("half").is_err());
    assert_eq!(m.resolve_set_partial("half").map(|r| r.nodes), Some(vec![1]));
    assert_eq!(m.resolve_set_partial("nope"), None);
}
