//! Whole-model checks against a consistent-mass axial bar's exact discrete equations.
use fieldcraft_model::harmonic::{Damping, Frequencies, Harmonic, LoadFactor, Method};
use fieldcraft_model::*;
use fieldcraft_solver::Limits;
use fieldcraft_solver::harmonic::solve_harmonic;
use std::collections::BTreeMap;
use std::f64::consts::PI;

fn rod(material: Damping, step_damping: Damping, method: Method, prescribed: f64) -> Model {
    let mut m = Model::default();
    m.set_node(1, [0.0; 3]).unwrap();
    m.set_node(2, [1.0, 0.0, 0.0]).unwrap();
    m.set_element(1, ElementKind::Rod2, vec![1, 2]).unwrap();
    for (name, nodes, elements) in [("all", vec![1, 2], vec![1]), ("root", vec![1], vec![]), ("tip", vec![2], vec![])] {
        m.add_set(name, Set { nodes, elements, ..Set::default() }).unwrap();
    }
    m.add_material(
        "bar",
        Material {
            elastic: Some(Elastic::Isotropic { young: 4.0, poisson: 0.25 }),
            density: Some(3.0),
            damping: Some(material),
            ..Material::default()
        },
    )
    .unwrap();
    m.add_section(
        "bar",
        Section {
            set: "all".into(),
            material: "bar".into(),
            kind: SectionKind::Rod { area: 1.0 },
            formulation: Default::default(),
            orientation: None,
        },
    )
    .unwrap();
    m.add_boundary_condition(
        "root",
        BoundaryCondition::Displacement { target: Target::Set("root".into()), values: [Some(prescribed), None, None, None, None, None] },
    )
    .unwrap();
    m.add_boundary_condition(
        "guide",
        BoundaryCondition::Displacement { target: Target::Set("all".into()), values: [None, Some(0.0), Some(0.0), None, None, None] },
    )
    .unwrap();
    m.add_load("force", Load::Force { target: Target::Set("tip".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    let h = Harmonic {
        frequencies: Frequencies::List(vec![0.0, 0.99 / PI, 1.0 / PI, 1.01 / PI]),
        method,
        damping: step_damping,
        load_factors: BTreeMap::from([("force".into(), LoadFactor { amplitude: 2.0, phase: 0.4 })]),
    };
    let mut step = Step::new("H", Procedure::Harmonic(h));
    step.boundary_conditions.extend(["root".into(), "guide".into()]);
    step.loads.insert("force".into());
    m.add_step(step, None).unwrap();
    m
}

use fieldcraft_model::design::{MaterialConstant, Parameter};
use fieldcraft_solver::sensitivity::*;
#[test]
fn complex_observations_identify_young_and_density() {
    let m = rod(Damping::default(), Damping { alpha: 0.2, beta: 0.01, loss_factor: 0.03 }, Method::Direct, 0.0);
    let parameters = vec![
        JointParameter {
            parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young },
            lower_si: 2.0,
            upper_si: 6.0,
            scale_si: 4.0,
        },
        JointParameter {
            parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Density },
            lower_si: 2.0,
            upper_si: 4.0,
            scale_si: 3.0,
        },
    ];
    let mut truth = m.clone();
    parameters[0].parameter.set(&mut truth, 4.2).unwrap();
    parameters[1].parameter.set(&mut truth, 3.1).unwrap();
    let target = solve_harmonic(&truth, Some("H"), &Limits::default()).unwrap();
    let mut observations = Vec::new();
    for (frame, v) in target.frames.iter().enumerate() {
        for (response, value) in [
            (JointResponse::HarmonicReal { step: "H".into(), node: 2, dof: 1, frame }, v.real[1][0]),
            (JointResponse::HarmonicImaginary { step: "H".into(), node: 2, dof: 1, frame }, v.imaginary[1][0]),
        ] {
            observations.push(JointObservation { response, value_si: value, standard_deviation_si: 0.01 });
        }
    }
    let fit = JointFit { parameters, observations, regularization: 0.0, max_iterations: 32, relative_increment: 1e-4 };
    let r = inverse_joint(&m, &fit, &Limits::default()).unwrap();
    for (a, b) in r.estimate_si.iter().zip([4.2, 3.1]) {
        assert!((a / b - 1.0).abs() < 1e-6, "{r:?}");
    }
}
#[test]
fn frequencies_only_cannot_separate_young_and_density_even_with_prior() {
    let p = |c| Parameter::Material { name: "bar".into(), constant: c };
    let fit = JointFit {
        parameters: vec![
            JointParameter { parameter: p(MaterialConstant::Young), lower_si: 1.0, upper_si: 10.0, scale_si: 4.0 },
            JointParameter { parameter: p(MaterialConstant::Density), lower_si: 1.0, upper_si: 10.0, scale_si: 3.0 },
        ],
        observations: vec![
            JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: 0 }, value_si: 1.0, standard_deviation_si: 0.01 },
            JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: 1 }, value_si: 2.0, standard_deviation_si: 0.01 },
        ],
        regularization: 100.0,
        max_iterations: 8,
        relative_increment: 1e-4,
    };
    let e = fit_joint(&[4.0, 3.0], &fit, |x| {
        let f = (x[0] / x[1]).sqrt();
        Ok((vec![f, 2.0 * f], vec![vec![f / (2.0 * x[0]), -f / (2.0 * x[1])], vec![f / x[0], -f / x[1]]]))
    })
    .unwrap_err()
    .to_string();
    assert!(e.contains("rank 1/2"), "{e}");
}

#[test]
fn homogeneous_joint_reuse_agrees_with_actual_finite_element_fit() {
    let mut m = rod(Damping::default(), Damping::default(), Method::Direct, 0.0);
    for (name, procedure) in [("S", Procedure::Static {}), ("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None })] {
        let mut step = Step::new(name, procedure);
        step.boundary_conditions.extend(["root".into(), "guide".into()]);
        if name == "S" {
            step.loads.insert("force".into());
        }
        m.add_step(step, None).unwrap();
    }
    let fit = JointFit {
        parameters: vec![
            JointParameter {
                parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young },
                lower_si: 2.0,
                upper_si: 6.0,
                scale_si: 4.0,
            },
            JointParameter {
                parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Density },
                lower_si: 2.0,
                upper_si: 4.0,
                scale_si: 3.0,
            },
        ],
        observations: vec![
            JointObservation {
                response: JointResponse::Static { step: "S".into(), node: 2, dof: 1 },
                value_si: 1.0 / 4.2,
                standard_deviation_si: 0.001,
            },
            JointObservation {
                response: JointResponse::Frequency { step: "F".into(), mode: 0 },
                value_si: (3.0_f64 * 4.2 / 3.1).sqrt() / (2.0 * PI),
                standard_deviation_si: 0.001,
            },
        ],
        regularization: 0.0,
        max_iterations: 32,
        relative_increment: 1e-4,
    };
    let a = inverse_joint(&m, &fit, &Limits::default()).unwrap();
    let b = inverse_joint_uniform_material(&m, &fit, &Limits::default()).unwrap();
    for ((x, y), truth) in a.estimate_si.iter().zip(&b.estimate_si).zip([4.2, 3.1]) {
        assert!((x / y - 1.0).abs() < 1e-6);
        assert!((y / truth - 1.0).abs() < 1e-6);
    }
}

#[test]
fn material_and_elastic_support_are_identified_together() {
    let mut m = rod(Damping::default(), Damping::default(), Method::Direct, 0.0);
    m.add_spring(
        "support",
        Spring::Spring { first: Target::Set("root".into()), second: None, direction: fieldcraft_model::discrete::Direction::Dof(1), stiffness: 4.0 },
    )
    .unwrap();
    for (name, procedure) in [("S", Procedure::Static {}), ("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None })] {
        let mut step = Step::new(name, procedure);
        step.boundary_conditions.insert("guide".into());
        if name == "S" {
            step.loads.insert("force".into());
        }
        m.add_step(step, None).unwrap();
    }
    let parameters = vec![
        JointParameter {
            parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young },
            lower_si: 2.0,
            upper_si: 8.0,
            scale_si: 4.0,
        },
        JointParameter { parameter: Parameter::SpringStiffness { name: "support".into() }, lower_si: 2.0, upper_si: 8.0, scale_si: 4.0 },
    ];
    let mut target = m.clone();
    parameters[0].parameter.set(&mut target, 4.2).unwrap();
    parameters[1].parameter.set(&mut target, 5.0).unwrap();
    let f = fieldcraft_solver::solve_frequency(&target, Some("F"), &Default::default(), &Limits::default()).unwrap().modes[0].frequency;
    let fit = JointFit {
        parameters,
        observations: vec![
            JointObservation { response: JointResponse::Static { step: "S".into(), node: 1, dof: 1 }, value_si: 0.2, standard_deviation_si: 0.0002 },
            JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: 0 }, value_si: f, standard_deviation_si: f * 0.001 },
        ],
        regularization: 0.0,
        max_iterations: 32,
        relative_increment: 1e-4,
    };
    let r = inverse_joint(&m, &fit, &Limits::default()).unwrap();
    for (x, truth) in r.estimate_si.iter().zip([4.2, 5.0]) {
        assert!((x / truth - 1.0).abs() < 1e-6, "{r:?}");
    }
    assert!(inverse_joint_uniform_material(&m, &fit, &Limits::default()).is_err());
}

#[test]
fn badly_conditioned_data_is_rejected_before_regularization() {
    let p = |name: &str| JointParameter { parameter: Parameter::SpringStiffness { name: name.into() }, lower_si: 0.1, upper_si: 10.0, scale_si: 1.0 };
    let fit = JointFit {
        parameters: vec![p("a"), p("b")],
        observations: (0..2)
            .map(|i| JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: i }, value_si: 0.0, standard_deviation_si: 1.0 })
            .collect(),
        regularization: 1000.0,
        max_iterations: 8,
        relative_increment: 1e-4,
    };
    let e = fit_joint(&[1.0, 1.0], &fit, |x| Ok((vec![x[0], 1e-9 * x[1]], vec![vec![1.0, 0.0], vec![0.0, 1e-9]]))).unwrap_err().to_string();
    assert!(e.contains("unidentifiable") && e.contains("weakest scaled parameter direction"), "{e}");
}

#[test]
fn step_damping_derivatives_and_joint_recovery() {
    use fieldcraft_model::design::DampingParameter as D;
    let m = rod(Damping::default(), Damping { alpha: 0.2, beta: 0.01, loss_factor: 0.03 }, Method::Direct, 0.0);
    for component in [D::Alpha, D::Beta, D::LossFactor] {
        let p = Parameter::Damping { step: "H".into(), component };
        let value = p.value(&m).unwrap();
        let h = value * 1e-4;
        let derivative = harmonic_sensitivity(&m, Some("H"), &p, h, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).unwrap();
        let mut plus = m.clone();
        let mut minus = m.clone();
        p.set(&mut plus, value + h).unwrap();
        p.set(&mut minus, value - h).unwrap();
        let rp = solve_harmonic(&plus, Some("H"), &Limits::default()).unwrap();
        let rm = solve_harmonic(&minus, Some("H"), &Limits::default()).unwrap();
        for ((d, p), m) in derivative.frames.iter().zip(&rp.frames).zip(&rm.frames) {
            for (i, fd) in [(0, (p.real[1][0] - m.real[1][0]) / (2.0 * h)), (1, (p.imaginary[1][0] - m.imaginary[1][0]) / (2.0 * h))] {
                assert!((d.direct[i] - fd).abs() < 1e-6 * fd.abs().max(1.0), "{:?}: {} {fd}", component, d.direct[i]);
            }
        }
    }
    let parameters = vec![
        JointParameter { parameter: Parameter::Damping { step: "H".into(), component: D::Alpha }, lower_si: 0.0, upper_si: 1.0, scale_si: 0.2 },
        JointParameter { parameter: Parameter::Damping { step: "H".into(), component: D::LossFactor }, lower_si: 0.0, upper_si: 0.2, scale_si: 0.03 },
    ];
    let mut target = m.clone();
    parameters[0].parameter.set(&mut target, 0.23).unwrap();
    parameters[1].parameter.set(&mut target, 0.04).unwrap();
    let truth = solve_harmonic(&target, Some("H"), &Limits::default()).unwrap();
    let mut observations = Vec::new();
    for (frame, v) in truth.frames.iter().enumerate() {
        observations.push(JointObservation {
            response: JointResponse::HarmonicReal { step: "H".into(), node: 2, dof: 1, frame },
            value_si: v.real[1][0],
            standard_deviation_si: 0.001,
        });
        observations.push(JointObservation {
            response: JointResponse::HarmonicImaginary { step: "H".into(), node: 2, dof: 1, frame },
            value_si: v.imaginary[1][0],
            standard_deviation_si: 0.001,
        });
    }
    let fit = JointFit { parameters, observations, regularization: 0.0, max_iterations: 32, relative_increment: 1e-5 };
    let r = inverse_joint(&m, &fit, &Limits::default()).unwrap();
    assert!((r.estimate_si[0] / 0.23 - 1.0).abs() < 1e-6, "{r:?}");
    assert!((r.estimate_si[1] / 0.04 - 1.0).abs() < 1e-6);
    let zero = Parameter::Damping { step: "H".into(), component: D::Alpha };
    let mut zero_model = m.clone();
    zero.set(&mut zero_model, 0.0).unwrap();
    assert!(harmonic_sensitivity(&zero_model, Some("H"), &zero, 1e-5, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).is_ok());
}

#[test]
fn full_modal_basis_derivatives_agree_with_direct_and_ratio_differences() {
    use fieldcraft_model::design::DampingParameter as D;
    let damping = Damping { alpha: 0.2, beta: 0.01, loss_factor: 0.03 };
    let direct = rod(Damping::default(), damping, Method::Direct, 0.0);
    let modal = rod(Damping::default(), damping, Method::Modal { modes: 1, damping_ratios: vec![], residual_vectors: 0 }, 0.0);
    for p in [
        Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young },
        Parameter::Material { name: "bar".into(), constant: MaterialConstant::Density },
        Parameter::Damping { step: "H".into(), component: D::Alpha },
        Parameter::Damping { step: "H".into(), component: D::Beta },
        Parameter::Damping { step: "H".into(), component: D::LossFactor },
    ] {
        let h = p.value(&direct).unwrap() * 1e-4;
        let a = harmonic_sensitivity(&direct, Some("H"), &p, h, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).unwrap();
        let b = harmonic_sensitivity(&modal, Some("H"), &p, h, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).unwrap();
        for (a, b) in a.frames.iter().zip(&b.frames) {
            for i in 0..2 {
                assert!((a.direct[i] - b.direct[i]).abs() < 1e-8 * a.direct[i].abs().max(1.0));
                assert!((b.direct[i] - b.adjoint[i]).abs() < 1e-8 * b.direct[i].abs().max(1.0));
            }
        }
    }
    let p = Parameter::Damping { step: "H".into(), component: D::ModalRatio { mode: 0 } };
    let mut modal = modal;
    p.set(&mut modal, 0.02).unwrap();
    let d = harmonic_sensitivity(&modal, Some("H"), &p, 1e-5, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).unwrap();
    let mut plus = modal.clone();
    let mut minus = modal.clone();
    p.set(&mut plus, 0.020001).unwrap();
    p.set(&mut minus, 0.019999).unwrap();
    let rp = solve_harmonic(&plus, Some("H"), &Limits::default()).unwrap();
    let rm = solve_harmonic(&minus, Some("H"), &Limits::default()).unwrap();
    for ((d, p), m) in d.frames.iter().zip(&rp.frames).zip(&rm.frames) {
        for (i, fd) in [(0, (p.real[1][0] - m.real[1][0]) / 2e-6), (1, (p.imaginary[1][0] - m.imaginary[1][0]) / 2e-6)] {
            assert!((d.direct[i] - fd).abs() < 1e-6 * fd.abs().max(1.0), "modal ratio component {i}: {} vs {fd}", d.direct[i]);
        }
    }
}

#[test]
fn modal_ratio_and_loss_are_identified_from_complex_frequency_data() {
    use fieldcraft_model::design::DampingParameter as D;
    let m = rod(
        Damping::default(),
        Damping { alpha: 0.1, beta: 0.0, loss_factor: 0.03 },
        Method::Modal { modes: 1, damping_ratios: vec![0.02], residual_vectors: 0 },
        0.0,
    );
    let parameters = vec![
        JointParameter {
            parameter: Parameter::Damping { step: "H".into(), component: D::ModalRatio { mode: 0 } },
            lower_si: 0.0,
            upper_si: 0.1,
            scale_si: 0.02,
        },
        JointParameter { parameter: Parameter::Damping { step: "H".into(), component: D::LossFactor }, lower_si: 0.0, upper_si: 0.1, scale_si: 0.03 },
    ];
    let truth = [0.025, 0.04];
    let mut target = m.clone();
    for (p, x) in parameters.iter().zip(truth) {
        p.parameter.set(&mut target, x).unwrap();
    }
    let response = solve_harmonic(&target, Some("H"), &Limits::default()).unwrap();
    let mut observations = Vec::new();
    for (frame, v) in response.frames.iter().enumerate() {
        observations.push(JointObservation {
            response: JointResponse::HarmonicReal { step: "H".into(), node: 2, dof: 1, frame },
            value_si: v.real[1][0],
            standard_deviation_si: 0.001,
        });
        observations.push(JointObservation {
            response: JointResponse::HarmonicImaginary { step: "H".into(), node: 2, dof: 1, frame },
            value_si: v.imaginary[1][0],
            standard_deviation_si: 0.001,
        });
    }
    let r = inverse_joint(
        &m,
        &JointFit { parameters, observations, regularization: 0.0, max_iterations: 32, relative_increment: 1e-5 },
        &Limits::default(),
    )
    .unwrap();
    for (x, t) in r.estimate_si.iter().zip(truth) {
        assert!((x / t - 1.0).abs() < 1e-6, "{r:?}");
    }
}
