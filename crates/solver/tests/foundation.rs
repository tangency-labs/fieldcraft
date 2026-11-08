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
