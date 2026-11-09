mod common;
use common::*;
use fieldcraft_model::design::{FoundationComponent, Parameter};
use fieldcraft_model::harmonic::{Damping, Frequencies, Harmonic, Method};
use fieldcraft_model::*;
use fieldcraft_solver::{
    Limits,
    sensitivity::{StaticObservation, static_sensitivity},
    solve_static,
};

fn bar() -> Grid {
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(1000., 0.), Some(2.), full(), Procedure::Static {});
    g.face_set("bed", |p| near(p[0], 0.));
    g.face_set("tip", |p| near(p[0], 1.));
    g.model
        .add_foundation(
            "Bed",
            Foundation { set: "bed".into(), normal: 100., tangential: 100., normal_damping: 0., tangential_damping: 0., loss_factor: 0. },
        )
        .unwrap();
    g.load("pressure", Load::Pressure { set: "tip".into(), magnitude: -1. });
    g
}
#[test]
fn supported_bar_and_exact_design_derivative() {
    let g = bar();
    let node = g.nodes_where(|p| near(p[0], 1.))[0];
    let result = solve_static(&g.model, None, &Limits::default()).unwrap();
    let i = result.nodes.binary_search(&node).unwrap();
    // Uniform axial pressure: u_tip=pL/E+p/k, unit area and length.
    assert!((result.displacement[i][0] - 0.011).abs() < 1e-12);
    let p = Parameter::Foundation { name: "Bed".into(), component: FoundationComponent::Normal };
    let d = static_sensitivity(&g.model, None, &p, 0.01, StaticObservation { node, dof: 1, stress_element: None }, &Limits::default()).unwrap();
    assert!((d.direct + 0.0001).abs() < 1e-12);
    assert!((d.adjoint - d.direct).abs() < 1e-13);
    let mut plus = g.model.clone();
    let mut minus = g.model.clone();
    p.set(&mut plus, 100.01).unwrap();
    p.set(&mut minus, 99.99).unwrap();
    let fd = (solve_static(&plus, None, &Limits::default()).unwrap().displacement[i][0]
        - solve_static(&minus, None, &Limits::default()).unwrap().displacement[i][0])
        / 0.02;
    assert!((fd - d.direct).abs() < 2e-12);
    let mut large = g.model.clone();
    p.set(&mut large, 1e9).unwrap();
    assert!((solve_static(&large, None, &Limits::default()).unwrap().displacement[i][0] - 0.001).abs() < 2e-9);
}

#[test]
fn nonlinear_bar_keeps_reference_foundation_force_and_balance() {
    // Independent axial StVK patch: P=E lambda (lambda^2-1)/2,
    // u_tip=P/k + (lambda-1)L. Foundation work is P^2/(2k).
    for finite in [false, true] {
        let mut g = bar();
        let mut step = g.model.step("Step").unwrap().clone();
        step.procedure = Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: finite, ..NonlinearStatic::default() });
        g.model.replace_step(step).unwrap();
        let r = fieldcraft_solver::solve_nonlinear(&g.model, None, &Limits::default()).unwrap();
        let mut stretch = 1.001_f64;
        for _ in 0..8 {
            stretch -= (500. * stretch * (stretch * stretch - 1.) - 1.) / (500. * (3. * stretch * stretch - 1.));
        }
        let extension = if finite { stretch - 1. } else { 0.001 };
        for (label, p) in &g.param {
            let i = r.nodes.binary_search(label).unwrap();
            assert!((r.displacement[i][0] - (0.01 + p[0] * extension)).abs() < 1e-10);
        }
        assert!((r.applied_total[0] - 1.).abs() < 1e-12);
        let bed_u = r.displacement[r.nodes.binary_search(&g.nodes_where(|p| near(p[0], 0.))[0]).unwrap()][0];
        assert!((r.applied_total[0] + r.reaction_total[0] - 100. * bed_u).abs() < 1e-9);
        assert!(!r.warnings.iter().any(|w| w.contains("loads and reactions differ")), "{:?}", r.warnings);
        let history = r.nonlinear.as_ref().unwrap();
        assert!(history.increments.iter().all(|i| i.converged));
        assert_eq!(history.increments.last().unwrap().load_factor, 1.);
        println!("finite={finite} bed_u={bed_u:.17e} stretch={stretch:.17e} history={history:?}");
    }
}

#[test]
fn reference_foundation_terminal_tangent_has_the_exact_rigid_translation_modes() {
    // A rigid unit cube of density 2 on a unit-area bed has M=2, K=k.
    // Translation causes no solid strain/stress stiffness. The support contributes
    // exactly once in static equilibrium and in both perturbation procedures.
    for direction in 0..3 {
        let mut g = grid([1, 1, 1], ElementKind::Hex8, |p| p);
        g.define(
            iso(1000., 0.),
            Some(2.),
            full(),
            Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: true, ..NonlinearStatic::default() }),
        );
        g.face_set("bed", |p| near(p[0], 0.));
        g.model
            .add_foundation(
                "Bed",
                Foundation { set: "bed".into(), normal: 100., tangential: 30., normal_damping: 0., tangential_damping: 0., loss_factor: 0. },
            )
            .unwrap();
        g.model.add_reference_point("Motion", ReferencePoint { xyz: [0.5; 3], label: 100 }).unwrap();
        g.model
            .add_coupling("Rigid", Coupling { reference_point: "Motion".into(), set: "all".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() })
            .unwrap();
        let mut values = [Some(0.); 6];
        values[direction] = None;
        g.bc("Guide", BoundaryCondition::Displacement { target: Target::ReferencePoint("Motion".into()), values });
        let mut vector = [0.; 3];
        vector[direction] = 1.;
        g.load("Probe", Load::Force { target: Target::ReferencePoint("Motion".into()), vector, distribution: Default::default() });
        let mut modes = Step::new("Modes", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
        modes.base_state = Some("Step".into());
        modes.boundary_conditions.insert("Guide".into());
        g.model.add_step(modes, None).unwrap();
        let mut response = Step::new(
            "Response",
            Procedure::Harmonic(Harmonic {
                frequencies: Frequencies::List(vec![0., 0.5]),
                method: Method::Direct,
                damping: Damping::default(),
                load_factors: Default::default(),
            }),
        );
        response.base_state = Some("Step".into());
        response.boundary_conditions.insert("Guide".into());
        response.loads.insert("Probe".into());
        g.model.add_step(response, None).unwrap();
        let source = g.model.clone();
        let limits = Limits::default();
        let k = if direction == 0 { 100. } else { 30. };
        let r = fieldcraft_solver::solve_nonlinear(&g.model, Some("Step"), &limits).unwrap();
        for u in &r.displacement {
            assert!((u[direction] - 1. / k).abs() < 1e-11);
        }
        let f = fieldcraft_solver::solve_frequency(&g.model, Some("Modes"), &Default::default(), &limits).unwrap();
        assert_eq!(f.stats.equations, 1);
        assert!((f.total_mass - 2.).abs() < 1e-12);
        assert!((f.modes[0].eigenvalue - k / 2.).abs() < 1e-9);
        let h = fieldcraft_solver::solve_harmonic(&g.model, Some("Response"), &limits).unwrap();
        for frame in h.frames {
            let omega = 2. * std::f64::consts::PI * frame.frequency;
            let expected = 1. / (k - 2. * omega * omega);
            assert!((frame.real[0][direction] - expected).abs() < 1e-10);
            assert!(frame.relative_residual < 1e-10);
        }
        assert_eq!(g.model, source);
    }
}
#[test]
fn harmonic_zero_frequency_keeps_foundation_loss_and_direct_matches_modal() {
    let mut g = bar();
    let mut f = g.model.foundation("Bed").unwrap().clone();
    f.normal_damping = 0.2;
    f.tangential_damping = 0.3;
    f.loss_factor = 0.1;
    g.model.replace_foundation("Bed", f).unwrap();
    let mut step = Step::new(
        "H",
        Procedure::Harmonic(Harmonic {
            frequencies: Frequencies::List(vec![0., 0.5]),
            method: Method::Direct,
            damping: Damping::default(),
            load_factors: Default::default(),
        }),
    );
    step.loads.insert("pressure".into());
    g.model.add_step(step.clone(), None).unwrap();
    let a = fieldcraft_solver::harmonic::solve_harmonic(&g.model, Some("H"), &Limits::default()).unwrap();
    let node = g.nodes_where(|p| near(p[0], 1.))[0];
    let i = a.nodes.binary_search(&node).unwrap();
    assert!((a.frames[0].real[i][0] - (0.001 + 1. / 101.)).abs() < 1e-11);
    assert!((a.frames[0].imaginary[i][0] + 0.1 / 101.).abs() < 1e-11);
    if let Procedure::Harmonic(h) = &mut step.procedure {
        h.method = Method::Modal { modes: 36, damping_ratios: vec![], residual_vectors: 0 };
    }
    g.model.replace_step(step).unwrap();
    let b = fieldcraft_solver::harmonic::solve_harmonic(&g.model, Some("H"), &Limits::default()).unwrap();
    for (a, b) in a.frames.iter().zip(&b.frames) {
        assert!((a.real[i][0] - b.real[i][0]).abs() < 1e-9);
        assert!((a.imaginary[i][0] - b.imaginary[i][0]).abs() < 1e-9);
    }
}
#[test]
fn foundation_persists_and_protects_its_face_set() {
    let mut g = bar();
    assert!(g.model.remove_set("bed").is_err());
    let mut f = g.model.foundation("Bed").unwrap().clone();
    f.normal = -1.;
    assert!(g.model.replace_foundation("Bed", f).is_err());
    assert_eq!(g.model.foundation("Bed").unwrap().normal, 100.);
}

#[test]
fn joint_static_observations_identify_young_and_foundation() {
    use fieldcraft_model::design::MaterialConstant;
    use fieldcraft_solver::sensitivity::*;
    let g = bar();
    let parameters = vec![
        JointParameter {
            parameter: Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young },
            lower_si: 500.,
            upper_si: 2000.,
            scale_si: 1000.,
        },
        JointParameter {
            parameter: Parameter::Foundation { name: "Bed".into(), component: FoundationComponent::Normal },
            lower_si: 50.,
            upper_si: 200.,
            scale_si: 100.,
        },
    ];
    let mut truth = g.model.clone();
    parameters[0].parameter.set(&mut truth, 1200.).unwrap();
    parameters[1].parameter.set(&mut truth, 130.).unwrap();
    let target = solve_static(&truth, None, &Limits::default()).unwrap();
    let observations = [0., 1.]
        .into_iter()
        .map(|x| {
            let node = g.nodes_where(|p| near(p[0], x))[0];
            let i = target.nodes.binary_search(&node).unwrap();
            JointObservation {
                response: JointResponse::Static { step: "Step".into(), node, dof: 1 },
                value_si: target.displacement[i][0],
                standard_deviation_si: 1e-5,
            }
        })
        .collect();
    let fit = JointFit { parameters, observations, regularization: 0., max_iterations: 32, relative_increment: 1e-4 };
    let result = inverse_joint(&g.model, &fit, &Limits::default()).unwrap();
    for (x, t) in result.estimate_si.iter().zip([1200., 130.]) {
        assert!((x / t - 1.).abs() < 1e-6, "{result:?}");
    }
}

#[test]
fn modal_and_complex_foundation_derivatives_match_response_differences() {
    use fieldcraft_solver::sensitivity::{HarmonicObservation, frequency_sensitivity, harmonic_sensitivity};
    use fieldcraft_solver::{FrequencyOptions, solve_frequency};
    let mut g = bar();
    g.bc("guide", BoundaryCondition::Displacement { target: Target::Set("all".into()), values: [None, Some(0.), Some(0.), None, None, None] });
    let mut fstep = Step::new("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    fstep.boundary_conditions.insert("guide".into());
    g.model.add_step(fstep, None).unwrap();
    let p = Parameter::Foundation { name: "Bed".into(), component: FoundationComponent::Normal };
    let options = FrequencyOptions::default();
    let limits = Limits::default();
    let d = frequency_sensitivity(&g.model, Some("F"), &p, 0.01, &options, &limits).unwrap();
    let mut plus = g.model.clone();
    let mut minus = g.model.clone();
    p.set(&mut plus, 100.01).unwrap();
    p.set(&mut minus, 99.99).unwrap();
    let fd = (solve_frequency(&plus, Some("F"), &options, &limits).unwrap().modes[0].frequency
        - solve_frequency(&minus, Some("F"), &options, &limits).unwrap().modes[0].frequency)
        / 0.02;
    // The serialized modal derivative field is checked below through its typed result.
    assert!((fd - d.modes[0].frequency_derivative).abs() < 1e-8 * fd.abs().max(1.));
    let mut f = g.model.foundation("Bed").unwrap().clone();
    f.loss_factor = 0.1;
    f.normal_damping = 0.2;
    g.model.replace_foundation("Bed", f).unwrap();
    let mut h = Step::new(
        "H",
        Procedure::Harmonic(Harmonic {
            frequencies: Frequencies::List(vec![0.5]),
            method: Method::Direct,
            damping: Damping::default(),
            load_factors: Default::default(),
        }),
    );
    h.loads.insert("pressure".into());
    h.boundary_conditions.insert("guide".into());
    g.model.add_step(h, None).unwrap();
    let node = g.nodes_where(|x| near(x[0], 1.))[0];
    let derivative = harmonic_sensitivity(&g.model, Some("H"), &p, 0.01, HarmonicObservation { node, dof: 1 }, &limits).unwrap();
    let mut plus = g.model.clone();
    let mut minus = g.model.clone();
    p.set(&mut plus, 100.01).unwrap();
    p.set(&mut minus, 99.99).unwrap();
    let a = fieldcraft_solver::harmonic::solve_harmonic(&plus, Some("H"), &limits).unwrap();
    let b = fieldcraft_solver::harmonic::solve_harmonic(&minus, Some("H"), &limits).unwrap();
    let i = a.nodes.binary_search(&node).unwrap();
    for (fd, exact) in [
        ((a.frames[0].real[i][0] - b.frames[0].real[i][0]) / 0.02, derivative.frames[0].direct[0]),
        ((a.frames[0].imaginary[i][0] - b.frames[0].imaginary[i][0]) / 0.02, derivative.frames[0].direct[1]),
    ] {
        assert!((fd - exact).abs() < 1e-10, "{fd} {exact}");
    }
}
