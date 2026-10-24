#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use super::*;
use fieldcraft_model::{Contact, ContactSurfaces, Elastic, ElementFace, ElementKind, Formulation, Material, Section, Set, Target};

/// Two original unit hexes, the lower one fixed, the upper one's lateral motion fixed and its
/// top displaced. Four free vertical degrees of freedom make a small real contact problem.
fn blocks() -> Model {
    let mut model = Model::default();
    let corners = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    for b in 0..2u32 {
        for (i, p) in corners.iter().enumerate() {
            model.set_node(8 * b + i as u32 + 1, [p[0], p[1], p[2] + b as f64]).unwrap();
        }
        model.set_element(b + 1, ElementKind::Hex8, (8 * b + 1..=8 * b + 8).collect()).unwrap();
    }
    model.add_set("all", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    model.add_set("base", Set { nodes: (1..=8).collect(), ..Set::default() }).unwrap();
    model.add_set("upper", Set { nodes: (9..=16).collect(), ..Set::default() }).unwrap();
    model.add_set("top", Set { nodes: (13..=16).collect(), ..Set::default() }).unwrap();
    model.add_set("main", Set { faces: vec![ElementFace { element: 1, face: 2 }], ..Set::default() }).unwrap();
    model.add_set("secondary", Set { faces: vec![ElementFace { element: 2, face: 1 }], ..Set::default() }).unwrap();
    model.add_material("Elastic", Material { elastic: Some(Elastic::Isotropic { young: 1000., poisson: 0. }), ..Material::default() }).unwrap();
    model
        .add_section(
            "Solid",
            Section { set: "all".into(), material: "Elastic".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() },
        )
        .unwrap();
    model
        .add_step(Step::new("Step", Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: false, ..NonlinearStatic::default() })), None)
        .unwrap();
    for (name, values) in [
        ("base", [Some(0.), Some(0.), Some(0.), None, None, None]),
        ("upper", [Some(0.), Some(0.), None, None, None, None]),
        ("top", [None, None, Some(-0.01), None, None, None]),
    ] {
        model.add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values }).unwrap();
        model.set_boundary_condition_active("Step", name, true).unwrap();
    }
    model.add_contact("Interface", Contact::new(ContactSurfaces::Pair { main: "main".into(), secondary: "secondary".into() })).unwrap();
    model
}

/// Test-only assembly of the normal Run, with no production hook or policy override. This
/// small linear-elastic case needs no thermal, follower, plastic or rotational setup.
fn with_run<T>(model: &Model, f: impl FnOnce(&mut Run<'_>) -> T) -> T {
    let limits = Limits::default();
    let resolved = model.resolve().unwrap();
    let step = model.step("Step").unwrap();
    let Procedure::NonlinearStatic(settings) = &step.procedure else { panic!("test step") };
    let mesh = Mesh::new(model, &resolved, limits.max_dofs).unwrap();
    let (graph, order) = pattern(&mesh);
    let mut warnings = Vec::new();
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let red = reduction(&setup, true, &order, budget(&limits, &mesh, &graph, 1), &mut warnings).unwrap();
    assert_eq!(red.neq(), 4);
    let kinds: Vec<_> = mesh.elements.iter().map(|e| e.kernel().unwrap().dofs().node).collect();
    let problem = Problem {
        mesh: &mesh,
        red: &red,
        dead: loads(model, &resolved, &mesh, step, &red.active).unwrap(),
        followers: Vec::new(),
        springs: None,
        finite: false,
        element_dofs: mesh.elements.iter().zip(&kinds).map(|(e, k)| mesh.layout.element_dofs(&e.nodes, *k).unwrap()).collect(),
        element_kinds: kinds,
        materials: mesh.elements.iter().map(|e| MaterialModel::from_elasticity(&e.d, None)).collect(),
        field: None,
        pattern: graph.zero_matrix(&mesh.layout),
        graph,
        contacts: Interface::all(model, &resolved, &mesh, false, &mut warnings).unwrap(),
    };
    let mut run = Run {
        problem: &problem,
        settings,
        limits: &limits,
        force_scale: 0.,
        force_floor: None,
        moment_scale: 0.,
        length: 6f64.sqrt(),
        prescribed: (red.largest_prescribed(), 0.),
        rotation: vec![false; mesh.dofs()],
        rotation_eq: vec![false; red.neq()],
        iterations: 0,
        factor_entries: 0,
        committed: vec![Vec::new(); mesh.elements.len()],
        path_dependent: false,
        undeformed: None,
        contact: problem.contacts.iter().map(Interface::committed).collect(),
        reference: vec![0.; mesh.dofs()],
        band: 0.,
        augmentations: 0,
        augmentation: 0,
        converged_scale: 0.,
        converged_moment: 0.,
    };
    f(&mut run)
}

fn solve_increment(
    run: &mut Run<'_>,
    q: &[f64],
    from: f64,
    lambda: f64,
    predictor: Option<&[f64]>,
) -> std::result::Result<(Converged, usize), Failure> {
    let mut iterations = Vec::new();
    run.augmentation = 0;
    let size = (lambda - from).abs();
    let c = run.increment(q, predictor, None, lambda, size, &mut iterations)?;
    let mut augmentations = 0;
    let c = run.augment(c, q, lambda, size, &mut iterations, &mut augmentations, &mut Vec::new())?;
    Ok((c, augmentations))
}

fn succeeded<T>(result: std::result::Result<T, Failure>) -> T {
    match result {
        Ok(x) => x,
        Err(Failure::Cutback(s)) => panic!("cutback: {s}"),
        Err(Failure::Fatal(e)) => panic!("fatal: {e}"),
    }
}

fn tangent_error(run: &Run<'_>, q: &[f64], lambda: f64) -> f64 {
    let e = succeeded(run.state(q, lambda, true));
    assert!(e.contact.iter().all(|c| c.pressure.iter().all(|p| *p == 0.)), "the control must be open");
    let [k] = run.problem.red.reduce([e.tangent.as_ref().unwrap()]).unwrap();
    let mut exact = vec![0.; q.len()];
    k.mul(&vec![1.; q.len()], &mut exact);
    let h = 1e-7;
    let a = succeeded(run.state(&q.iter().map(|x| x - h).collect::<Vec<_>>(), lambda, false));
    let b = succeeded(run.state(&q.iter().map(|x| x + h).collect::<Vec<_>>(), lambda, false));
    let difference: Vec<_> = b.internal.iter().zip(&a.internal).map(|(b, a)| (b - a) / (2. * h)).collect();
    let numeric = run.problem.red.restrict(&difference);
    exact.iter().zip(numeric).map(|(a, b)| (a - b).abs()).fold(0., f64::max)
}

#[test]
fn a_committed_contact_increment_does_not_carry_its_chatter_band_into_the_next() {
    for seeded in [false, true] {
        with_run(&blocks(), |run| {
            // Controlled entry to exactly the retry state left by active-set chatter; this avoids
            // relying on platform rounding to create a small naturally chattering contact patch.
            if seeded {
                run.band = contact::ENGAGED_BAND;
                let failed = run.increment(&[0.; 4], None, None, -0.05, 0.05, &mut Vec::new());
                assert!(matches!(failed, Err(Failure::Cutback(_))), "exercise an actual failed retry");
                assert_eq!(run.band, contact::ENGAGED_BAND, "a failed retry must not clear the band");
            }
            let first = succeeded(solve_increment(run, &[0.; 4], 0., 0.5, None));
            assert!(first.1 > 0, "exercise the actual Uzawa loop");
            assert_eq!(run.band, if seeded { contact::ENGAGED_BAND } else { 0. }, "retries and augmentations retain the band");
            let open = vec![0.0005; 4];
            let stale_error = tangent_error(run, &open, -0.05);
            if seeded {
                assert!(stale_error > 100., "prove the open state's stale-band tangent mismatch: {stale_error:e}");
            } else {
                assert!(stale_error < 1e-5);
            }
            let q = first.0.q;
            run.commit_increment(&q, 0.5, &first.0.evaluation);
            let error = tangent_error(run, &open, -0.05);
            eprintln!("after accepted increment: band {} near-open tangent error {error:e}", run.band);
            let (predictor, _) = succeeded(run.tangent_predictor(&q, 0.5, -0.05)).unwrap();
            let second = succeeded(solve_increment(run, &q, 0.5, -0.05, Some(&predictor)));
            assert!(error < 1e-5, "a later open contact has a force-inconsistent tangent: {error:e}");
            assert!(second.0.evaluation.contact.iter().all(|c| c.pressure.iter().all(|p| p.abs() < 1e-8)));
            assert!(second.0.q.iter().all(|x| (*x - 0.0005).abs() < 1e-8), "{:?}", second.0.q);
        });
    }
}
