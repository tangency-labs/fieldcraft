use super::*;
use fieldcraft_model::{Elastic, Formulation, Material, NonlinearStatic, Procedure, Section, SectionKind, Set, Step};

fn cubes(finite: bool, friction: Friction) -> (Mesh, Interface) {
    let mut model = Model::default();
    let cube = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    for block in 0..2 {
        for (i, x) in cube.iter().enumerate() {
            model.set_node(8 * block + i as u32 + 1, [x[0], x[1], x[2] + f64::from(block)]).unwrap();
        }
        model.set_element(block + 1, ElementKind::Hex8, (1..=8).map(|i| i + 8 * block).collect()).unwrap();
    }
    model.add_set("all", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    model.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1000., poisson: 0.0 }), ..Material::default() }).unwrap();
    model
        .add_section("solid", Section { set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    let top = ElementKind::Hex8.faces().iter().position(|face| face.iter().all(|i| cube[*i][2] == 1.0)).unwrap();
    model.add_set("main", Set { faces: vec![ElementFace { element: 1, face: top as u8 + 1 }], ..Set::default() }).unwrap();
    model.add_set("secondary", Set { faces: (1..=6).map(|face| ElementFace { element: 2, face }).collect(), ..Set::default() }).unwrap();
    let mut contact = Contact::new(ContactSurfaces::Pair { main: "main".into(), secondary: "secondary".into() });
    contact.friction = friction;
    model.add_contact("Interface", contact).unwrap();
    model.add_step(Step::new("Step", Procedure::NonlinearStatic(NonlinearStatic::default())), None).unwrap();
    let resolved = model.resolve().unwrap();
    let mesh = Mesh::new(&model, &resolved, crate::Limits::default().max_dofs).unwrap();
    let mut warnings = Vec::new();
    let interface = Interface::all(&model, &resolved, &mesh, finite, &mut warnings).unwrap().remove(0);
    assert!(warnings.is_empty(), "{warnings:?}");
    (mesh, interface)
}

fn tilted(mesh: &Mesh) -> Vec<f64> {
    let mut u = vec![0.0; mesh.dofs()];
    for (i, x) in mesh.xyz.iter().enumerate().filter(|(i, _)| mesh.labels[*i] >= 9) {
        let dofs = translations(mesh, i).unwrap();
        u[dofs[0]] = 0.001;
        u[dofs[2]] = -0.02 + 0.005 * x[0];
    }
    u
}

#[test]
fn small_sliding_reports_the_reference_pair_normal() {
    let (mesh, interface) = cubes(false, Friction::NONE);
    let u = tilted(&mesh);
    let reference = vec![0.0; mesh.dofs()];
    let mut internal = reference.clone();
    let mut blocks = Vec::new();
    let state = interface.evaluate(&mesh, &u, &reference, &interface.committed(), 0.0, true, &mut internal, &mut blocks).unwrap();
    let result = interface.result(&mesh, &state, 0);
    assert!(result.normal_force > 0.0 && !blocks.is_empty());
    assert_eq!(result.normal_force, result.force[2]);
    for (i, node) in interface.nodes.iter().enumerate() {
        if mesh.labels[*node] <= 12 {
            assert_eq!(state.report_normal[i], [0.0, 0.0, -1.0]);
            assert_ne!(state.normal[i], state.report_normal[i], "uncontacted sides still belong to the unchanged friction normal");
        } else {
            assert_eq!(state.report_normal[i], [0.0; 3]);
            assert_eq!(result.pressure[i], 0.0);
        }
    }
}

#[test]
fn finite_sliding_reports_the_current_overlap_normal() {
    let (mesh, interface) = cubes(true, Friction::NONE);
    let u = tilted(&mesh);
    let reference = vec![0.0; mesh.dofs()];
    let mut internal = reference.clone();
    let state = interface.evaluate(&mesh, &u, &reference, &interface.committed(), 0.0, false, &mut internal, &mut Vec::new()).unwrap();
    let direction = unit([0.005, 0.0, -1.0]).unwrap();
    for (i, support) in state.support.iter().enumerate().filter(|(_, a)| **a > 0.0) {
        assert!(support.is_finite());
        assert!((dot(state.report_normal[i], direction) - 1.0).abs() < 1e-12);
        assert_ne!(state.report_normal[i], [0.0, 0.0, -1.0]);
    }
}

#[test]
fn reporting_supports_are_weighted_and_missing_or_cancelling_supports_have_no_direction() {
    for scale in [1e-200, 1.0, 1e200] {
        let found = |normal, area| Found { main: 0, gap: [0.0; 4], support: [area, 0.0, 0.0, 0.0], normal };
        let pairs = vec![vec![found([1., 0., 0.], scale), found([0., 1., 0.], 3.0 * scale)]];
        let n = report_normals(&pairs, &[vec![0]], &[4.0 * scale, 0.0]);
        let direction = unit([1., 3., 0.]).unwrap();
        assert!((dot(n[0], direction) - 1.0).abs() < 1e-12);
        assert_eq!(n[1], [0.0; 3]);
        let opposing = vec![vec![found([1., 0., 0.], scale), found([-1., 0., 0.], scale)]];
        assert_eq!(report_normals(&opposing, &[vec![0]], &[2.0 * scale]), vec![[0.0; 3]]);
        assert_eq!(report_normals(&pairs, &[vec![0]], &[0.0]), vec![[0.0; 3]]);
    }
    let (mesh, interface) = cubes(false, Friction::coulomb(0.3));
    let u = tilted(&mesh);
    let reference = vec![0.0; mesh.dofs()];
    let mut internal = reference.clone();
    let mut state = interface.evaluate(&mesh, &u, &reference, &interface.committed(), 0.0, true, &mut internal, &mut Vec::new()).unwrap();
    // A cancelling support sum has no reporting direction. Preserve the full physical force
    // as finite shear rather than manufacturing a normal or dividing by zero.
    state.report_normal.fill([0.0; 3]);
    let result = interface.result(&mesh, &state, 0);
    assert_eq!(result.normal_force, 0.0);
    assert!(result.pressure.iter().all(|p| *p == 0.0));
    assert!(result.shear.iter().flatten().all(|v| v.is_finite()));
    assert!(result.force.iter().any(|f| *f != 0.0));
}
