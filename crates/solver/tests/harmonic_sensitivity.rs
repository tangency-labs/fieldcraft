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

#[test]
fn physical_material_derivatives_include_damping_mass_and_moving_supports() {
    use fieldcraft_model::design::{MaterialConstant, Parameter};
    use fieldcraft_solver::sensitivity::{HarmonicObservation, harmonic_sensitivity};
    let damping = Damping { alpha: 0.02, beta: 0.005, loss_factor: 0.03 };
    for prescribed in [0.0, 0.1] {
        let model = rod(damping, Damping::default(), Method::Direct, prescribed);
        for (constant, value) in [(MaterialConstant::Young, 4.0), (MaterialConstant::Density, 3.0)] {
            let parameter = Parameter::Material { name: "bar".into(), constant };
            let observation = HarmonicObservation { node: 2, dof: 1 };
            for h in [1e-4, 1e-5, 1e-6] {
                let result = harmonic_sensitivity(&model, None, &parameter, h, observation, &Limits::default()).unwrap();
                let (mut plus, mut minus) = (model.clone(), model.clone());
                parameter.set(&mut plus, value + h).unwrap();
                parameter.set(&mut minus, value - h).unwrap();
                let a = solve_harmonic(&plus, None, &Limits::default()).unwrap();
                let b = solve_harmonic(&minus, None, &Limits::default()).unwrap();
                for ((r, a), b) in result.frames.iter().zip(a.frames).zip(b.frames) {
                    let fd = [(a.real[1][0] - b.real[1][0]) / (2.0 * h), (a.imaginary[1][0] - b.imaginary[1][0]) / (2.0 * h)];
                    for j in 0..2 {
                        assert!((r.direct[j] - fd[j]).abs() < 1e-5 * r.direct[j].abs().max(1.0), "{r:?} {fd:?}");
                        assert!((r.direct[j] - r.adjoint[j]).abs() < 1e-10 * r.direct[j].abs().max(1.0));
                    }
                    assert!(r.primal_relative_residual < 1e-12);
                    assert!(r.tangent_relative_residual < 1e-12);
                }
            }
        }
    }
}

#[test]
fn undamped_pole_and_zero_response_do_not_report_derivatives() {
    use fieldcraft_model::design::{MaterialConstant, Parameter};
    use fieldcraft_solver::sensitivity::{HarmonicObservation, harmonic_sensitivity};
    let parameter = Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young };
    let model = rod(Damping::default(), Damping::default(), Method::Direct, 0.0);
    assert!(harmonic_sensitivity(&model, None, &parameter, 1e-4, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).is_err());
    let model = rod(Damping { loss_factor: 0.01, ..Damping::default() }, Damping::default(), Method::Direct, 0.0);
    assert!(harmonic_sensitivity(&model, None, &parameter, 1e-4, HarmonicObservation { node: 1, dof: 1 }, &Limits::default()).is_err());
}

#[test]
fn harmonic_inverse_recovers_material_with_complex_observations() {
    use fieldcraft_model::design::{MaterialConstant, Parameter};
    use fieldcraft_solver::sensitivity::{FrequencyFit, HarmonicObservation, inverse_harmonic};
    let model = rod(Damping { alpha: 0.02, beta: 0.005, loss_factor: 0.03 }, Damping::default(), Method::Direct, 0.0);
    let parameter = Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young };
    let mut truth = model.clone();
    parameter.set(&mut truth, 4.1).unwrap();
    let target = solve_harmonic(&truth, None, &Limits::default()).unwrap();
    let observed: Vec<f64> = target.frames.iter().flat_map(|f| [f.real[1][0], f.imaginary[1][0]]).collect();
    let fit = FrequencyFit {
        lower: 3.5,
        upper: 4.5,
        standard_deviation: Some(vec![0.01; observed.len()]),
        observed,
        mode_indices: None,
        regularization: 0.0,
        max_iterations: 32,
        relative_increment: 1e-5,
    };
    let r = inverse_harmonic(&model, None, &parameter, &fit, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).unwrap();
    assert!((r.estimate / 4.1 - 1.0).abs() < 1e-6, "{r:?}");
    let ci = r.confidence_interval.unwrap();
    assert!(ci[0] < 4.1 && ci[1] > 4.1);
    assert_eq!(parameter.value(&model).unwrap(), 4.0);
    let mut noisy = fit.clone();
    noisy.observed[0] += 0.002;
    let r = inverse_harmonic(&model, None, &parameter, &noisy, HarmonicObservation { node: 2, dof: 1 }, &Limits::default()).unwrap();
    let ci = r.confidence_interval.unwrap();
    assert!(ci[0] < 4.1 && ci[1] > 4.1, "{r:?}");
}

#[test]
fn harmonic_batch_matches_row_calls_bitwise() {
    use fieldcraft_model::design::{MaterialConstant, Parameter};
    use fieldcraft_solver::sensitivity::{HarmonicSensitivityRequest, harmonic_sensitivity_batch};
    let model = rod(Damping { alpha: 0.02, beta: 0.005, loss_factor: 0.03 }, Damping::default(), Method::Direct, 0.0);
    let parameter = Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young };
    let requests = vec![
        HarmonicSensitivityRequest {
            step: Some("H".into()),
            parameter: parameter.clone(),
            increment: 1e-4,
            observation: fieldcraft_solver::sensitivity::HarmonicObservation { node: 2, dof: 1 },
        },
        HarmonicSensitivityRequest {
            step: Some("H".into()),
            parameter: parameter.clone(),
            increment: 1e-4,
            observation: fieldcraft_solver::sensitivity::HarmonicObservation { node: 2, dof: 1 },
        },
    ];
    let batch = harmonic_sensitivity_batch(&model, &requests, &Limits::default()).unwrap();
    for (request, got) in requests.iter().zip(batch.results.iter()) {
        let row = fieldcraft_solver::sensitivity::harmonic_sensitivity(
            &model,
            request.step.as_deref(),
            &request.parameter,
            request.increment,
            request.observation,
            &Limits::default(),
        )
        .unwrap();
        assert_eq!(format!("{:?}", got.frames), format!("{:?}", row.frames));
    }
}
