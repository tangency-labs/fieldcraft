//! Independent continuum patch: two unit cubes, E=1e6 Pa, nu=0.
//! q = F (L1+L2)/(EA); U = F q/2. No external solver fixture.
use fieldcraft_model::*;
use fieldcraft_solver::{Limits, StaticResult, solve_static};
mod common;

#[test]
fn reactions_on_a_private_face_are_returned_on_the_source_nodes() {
    let mut source = model(Some(1000.));
    source.add_boundary_condition("CutHeld", BoundaryCondition::Fixed { target: Target::Set("CutFace".into()) }).unwrap();
    let mut step = source.step("Preload").unwrap().clone();
    step.boundary_conditions.insert("CutHeld".into());
    source.replace_step(step).unwrap();
    let r = solve_static(&source, None, &Limits::default()).unwrap();
    let cut: f64 = r.nodes.iter().zip(&r.reaction).filter(|(n, _)| (5..=8).contains(*n)).map(|(_, f)| f[2]).sum();
    assert!((cut + 1000.).abs() < 1e-9, "cut reaction {cut}");
    let sum: f64 = r.reaction.iter().map(|f| f[2]).sum();
    assert!((sum - r.reaction_total[2]).abs() < 1e-9);
    assert!((r.pretension[0].adjustment / 0.001 - 1.).abs() < 1e-12);
}

fn model(force: Option<f64>) -> Model {
    let mut m = Model::default();
    for z in 0..3u32 {
        for (i, [x, y]) in [[0., 0.], [1., 0.], [1., 1.], [0., 1.]].into_iter().enumerate() {
            m.set_node(4 * z + i as u32 + 1, [x, y, f64::from(z)]).unwrap();
        }
    }
    m.set_element(1, ElementKind::Hex8, (1..=8).collect()).unwrap();
    m.set_element(2, ElementKind::Hex8, (5..=12).collect()).unwrap();
    for (name, nodes) in [("Bottom", vec![1, 2, 3, 4]), ("Top", vec![9, 10, 11, 12]), ("Ends", vec![1, 2, 3, 4, 9, 10, 11, 12])] {
        m.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }
    m.add_set("CutFace", Set { faces: vec![ElementFace { element: 1, face: 2 }], ..Set::default() }).unwrap();
    m.add_set("All", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    m.add_material("Elastic", Material { elastic: Some(Elastic::Isotropic { young: 1e6, poisson: 0. }), density: Some(1.), ..Material::default() })
        .unwrap();
    m.add_section(
        "Solid",
        Section {
            set: "All".into(),
            material: "Elastic".into(),
            kind: SectionKind::Solid {},
            orientation: None,
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    m.add_boundary_condition("Ends", BoundaryCondition::Fixed { target: Target::Set("Ends".into()) }).unwrap();
    let mut step = Step::new("Preload", Procedure::Static {});
    step.boundary_conditions.insert("Ends".into());
    if let Some(f) = force {
        m.add_reference_point("Control", ReferencePoint { xyz: [0.5, 0.5, 1.], label: 100 }).unwrap();
        m.add_pretension(
            "Cut",
            PretensionSection {
                first: "Bottom".into(),
                second: "Top".into(),
                face: Some("CutFace".into()),
                normal: [0., 0., 1.],
                control: "Control".into(),
            },
        )
        .unwrap();
        step.pretension.insert("Cut".into(), PretensionControl::Force { value: f });
    }
    m.add_step(step, None).unwrap();
    m
}

fn bits(r: &StaticResult) -> Vec<u64> {
    r.displacement.iter().flatten().chain(r.reaction.iter().flatten()).chain(r.stress.iter().flatten()).map(|x| x.to_bits()).collect()
}

#[test]
fn continuous_cut_zero_force_matches_continuum_bit_for_bit() {
    let plain = solve_static(&model(None), None, &Limits::default()).unwrap();
    let source = model(Some(0.));
    let before = source.clone();
    let cut = solve_static(&source, None, &Limits::default()).unwrap();
    assert_eq!(source, before);
    assert_eq!(cut.nodes, plain.nodes);
    assert_eq!(bits(&cut), bits(&plain));
    assert_eq!(cut.pretension[0].adjustment, 0.);
    assert_eq!(cut.pretension[0].section_force, 0.);
}

#[test]
fn continuous_cut_force_reactions_and_element_energy_balance() {
    let source = model(Some(1000.));
    let before = source.clone();
    let r = solve_static(&source, None, &Limits::default()).unwrap();
    assert_eq!(source, before);
    assert_eq!(r.nodes, (1..=12).collect::<Vec<_>>());
    let q = r.pretension[0].adjustment;
    let force = r.pretension[0].section_force;
    assert!((q / 0.002 - 1.).abs() < 1e-12, "q={q}");
    assert!((force / 1000. - 1.).abs() < 1e-12, "F={force}");
    let bottom: f64 = r.nodes.iter().zip(&r.reaction).filter(|(n, _)| **n <= 4).map(|(_, v)| v[2]).sum();
    let top: f64 = r.nodes.iter().zip(&r.reaction).filter(|(n, _)| **n >= 9).map(|(_, v)| v[2]).sum();
    assert!((bottom + 1000.).abs() < 1e-9, "bottom={bottom}");
    assert!((top - 1000.).abs() < 1e-9, "top={top}");
    assert!(r.applied_total.iter().chain(&r.reaction_total).all(|v| v.abs() < 1e-9));
    // Unit cube, 2x2x2 Gauss weights and Jacobian: dV=1/8 at each point.
    let mut energy = 0.;
    for e in &r.elements {
        assert_eq!(e.points.len(), 8);
        for p in &e.points {
            assert!((p.stress[2] - 1000.).abs() < 1e-9);
            energy += 0.5 * p.stress.iter().zip(p.strain).map(|(s, e)| s * e).sum::<f64>() / 8.;
        }
    }
    assert!((energy - 0.5 * 1000. * q).abs() < 1e-12, "U={energy}");
    // Both outer ends remain fixed. Their relative displacement is NOT q.
    for (n, u) in r.nodes.iter().zip(&r.displacement) {
        if *n <= 4 || *n >= 9 {
            assert_eq!(*u, [0.; 3]);
        }
    }
    println!("q={q:.17e} section_force={force:.17e} bottom={bottom:.17e} top={top:.17e} energy={energy:.17e}");
}

#[test]
fn continuous_cut_rejects_boundary_face_and_nonplanar_axis() {
    let mut boundary = model(Some(1000.));
    boundary.replace_set("CutFace", Set { faces: vec![ElementFace { element: 2, face: 2 }], ..Set::default() }).unwrap();
    let before = boundary.clone();
    let error = solve_static(&boundary, None, &Limits::default()).unwrap_err().to_string();
    assert!(error.contains("exactly two connected solid neighbours"), "{error}");
    assert_eq!(before, boundary);
    let mut tilted = model(Some(1000.));
    let mut section = tilted.pretensions()["Cut"].clone();
    section.normal = [1., 0., 0.];
    tilted.replace_pretension("Cut", section).unwrap();
    let error = solve_static(&tilted, None, &Limits::default()).unwrap_err().to_string();
    assert!(error.contains("not planar"), "{error}");
}

#[test]
fn continuous_cut_small_strain_nonlinear_uses_the_same_equilibrium_and_projects_frames() {
    let mut source = model(Some(1000.));
    let linear = solve_static(&source, None, &Limits::default()).unwrap();
    let mut step = source.step("Preload").unwrap().clone();
    step.procedure = Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: false, ..NonlinearStatic::default() });
    source.replace_step(step).unwrap();
    let r = fieldcraft_solver::solve_nonlinear(&source, None, &Limits::default()).unwrap();
    assert_eq!(r.nodes, linear.nodes);
    assert!((r.pretension[0].adjustment - linear.pretension[0].adjustment).abs() < 1e-12);
    assert!((r.pretension[0].section_force - 1000.).abs() < 1e-8);
    for frame in &r.nonlinear.as_ref().unwrap().frames {
        assert_eq!(frame.displacement.len(), 12);
        assert!(frame.stress.is_empty() || frame.stress.len() == 12);
        assert!(frame.reaction.iter().all(|(i, _)| *i < 12));
    }
}

#[test]
fn continuous_multiface_cuts_include_quadratic_nodes_and_reject_partial_sections() {
    for kind in [ElementKind::Hex8, ElementKind::Tet10] {
        let mut g = common::grid([2, 1, 2], kind, |x| [2. * x[0], x[1], 2. * x[2]]);
        g.define(common::iso(1e6, 0.), Some(1.), Formulation::default(), Procedure::Static {});
        g.node_set("Bottom", |x| common::near(x[2], 0.));
        g.node_set("Top", |x| common::near(x[2], 1.));
        g.fixed("Ends", |x| common::near(x[2], 0.) || common::near(x[2], 1.));
        let faces: Vec<_> = g
            .faces_where(|x| common::near(x[2], 0.5))
            .into_iter()
            .filter(|f| g.model.element(f.element).unwrap().nodes.iter().any(|n| g.model.node(*n).unwrap().xyz[2] < 1. - 1e-9))
            .collect();
        g.model.add_set("CutFace", Set { faces: faces.clone(), ..Set::default() }).unwrap();
        g.model.add_reference_point("Control", ReferencePoint { xyz: [1., 0.5, 1.], label: 1000 }).unwrap();
        g.model
            .add_pretension(
                "Cut",
                PretensionSection {
                    first: "Bottom".into(),
                    second: "Top".into(),
                    face: Some("CutFace".into()),
                    normal: [0., 0., 1.],
                    control: "Control".into(),
                },
            )
            .unwrap();
        let mut step = g.model.step("Step").unwrap().clone();
        step.pretension.insert("Cut".into(), PretensionControl::Force { value: 1000. });
        g.model.replace_step(step).unwrap();
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        assert!((r.pretension[0].adjustment / 0.001 - 1.).abs() < 1e-11, "{kind:?}: {:?}", r.pretension);
        assert!((r.pretension[0].section_force / 1000. - 1.).abs() < 1e-11);
        g.model.replace_set("CutFace", Set { faces: vec![faces[0]], ..Set::default() }).unwrap();
        let e = solve_static(&g.model, None, &Limits::default()).unwrap_err().to_string();
        assert!(e.contains("incomplete cut"), "{kind:?}: {e}");
    }
}

#[test]
fn continuous_cut_uses_geometric_axis_and_rejects_an_inward_normal() {
    let mut rotated = model(Some(1000.));
    let nodes = rotated.nodes().clone();
    for (n, p) in nodes {
        rotated.set_node(n, [p.xyz[2], p.xyz[0], p.xyz[1]]).unwrap();
    }
    let mut cut = rotated.pretensions()["Cut"].clone();
    cut.normal = [1., 0., 0.];
    rotated.replace_pretension("Cut", cut.clone()).unwrap();
    let r = solve_static(&rotated, None, &Limits::default()).unwrap();
    assert!((r.pretension[0].adjustment / 0.002 - 1.).abs() < 1e-12);
    cut.normal = [-1., 0., 0.];
    rotated.replace_pretension("Cut", cut).unwrap();
    let error = solve_static(&rotated, None, &Limits::default()).unwrap_err().to_string();
    assert!(error.contains("normal must point away"), "{error}");
}
