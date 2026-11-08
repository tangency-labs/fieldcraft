use fieldcraft_model::{
    design::{MaterialConstant, Parameter},
    *,
};
use fieldcraft_solver::{Limits, sensitivity::*, solve_frequency, solve_static};
fn strip() -> Model {
    let mut m = Model::default();
    let cells = 24;
    for i in 0..=cells {
        for j in 0..2 {
            m.set_node(2 * i + j + 1, [2.0 * i as f64 / cells as f64, 0.1 * j as f64, 0.0]).unwrap();
        }
    }
    for i in 0..cells {
        m.set_element(i + 1, ElementKind::Quad4, vec![2 * i + 1, 2 * i + 3, 2 * i + 4, 2 * i + 2]).unwrap();
    }
    for (name, nodes, elements) in [
        ("all", (1..=2 * (cells + 1)).collect(), (1..=cells).collect()),
        ("root", vec![1, 2], vec![]),
        ("tip", vec![2 * cells + 1, 2 * cells + 2], vec![]),
    ] {
        m.add_set(name, Set { nodes, elements, ..Set::default() }).unwrap();
    }
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1000.0, poisson: 0.0 }), density: Some(2.0), ..Material::default() })
        .unwrap();
    m.add_section(
        "strip",
        Section {
            set: "all".into(),
            material: "Mat".into(),
            kind: SectionKind::Shell { thickness: 0.01, offset: 0.0, points: 3 },
            formulation: Default::default(),
            orientation: None,
        },
    )
    .unwrap();
    m.add_boundary_condition("root", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    m.add_load("force", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, 1e-9], distribution: Distribution::Total }).unwrap();
    for (name, procedure) in [("S", Procedure::Static {}), ("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None })] {
        let mut step = Step::new(name, procedure);
        step.boundary_conditions.insert("root".into());
        if name == "S" {
            step.loads.insert("force".into());
        }
        m.add_step(step, None).unwrap();
    }
    m
}

#[test]
fn actual_beam_young_and_thickness_joint_fit() {
    let m = strip();
    let parameters = vec![
        JointParameter {
            parameter: Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young },
            lower_si: 800.0,
            upper_si: 1400.0,
            scale_si: 1000.0,
        },
        JointParameter { parameter: Parameter::ShellThickness { name: "strip".into() }, lower_si: 0.008, upper_si: 0.014, scale_si: 0.01 },
    ];
    let truth = [1100.0, 0.011];
    let mut target = m.clone();
    for (p, x) in parameters.iter().zip(truth) {
        p.parameter.set(&mut target, x).unwrap();
    }
    let s = solve_static(&target, Some("S"), &Limits::default()).unwrap();
    let u = s.displacement[s.nodes.binary_search(&49).unwrap()][2];
    let f = solve_frequency(&target, Some("F"), &Default::default(), &Limits::default()).unwrap().modes[0].frequency;
    let mut fit = JointFit {
        parameters,
        observations: vec![
            JointObservation {
                response: JointResponse::Static { step: "S".into(), node: 49, dof: 3 },
                value_si: u,
                standard_deviation_si: u.abs() * 0.001,
            },
            JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: 0 }, value_si: f, standard_deviation_si: f * 0.001 },
        ],
        regularization: 0.0,
        max_iterations: 32,
        relative_increment: 1e-4,
    };
    for noisy in [false, true] {
        if noisy {
            fit.observations[0].value_si += 0.2 * fit.observations[0].standard_deviation_si;
            fit.observations[1].value_si -= 0.3 * fit.observations[1].standard_deviation_si;
        }
        let r = inverse_joint(&m, &fit, &Limits::default()).unwrap();
        let d = [r.estimate_si[0] - truth[0], r.estimate_si[1] - truth[1]];
        if !noisy {
            for (i, t) in truth.iter().enumerate() {
                assert!((r.estimate_si[i] / t - 1.0).abs() < 1e-6, "{r:?}");
            }
        }
        let c = r.covariance_si.unwrap();
        let det = c[0][0] * c[1][1] - c[0][1] * c[1][0];
        let mahal = (c[1][1] * d[0] * d[0] - 2.0 * c[0][1] * d[0] * d[1] + c[0][0] * d[1] * d[1]) / det;
        assert!(mahal < r.joint_chi_square, "{mahal}");
    }
}

#[test]
fn capped_ritz_values_are_not_used_as_sensitivity_or_joint_data() {
    let m = strip();
    let limits = Limits { max_lanczos_vectors: 6, ..Limits::default() };
    let young = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young };
    let error = frequency_sensitivity(&m, Some("F"), &young, 0.1, &Default::default(), &limits).unwrap_err().to_string();
    assert!(error.contains("Lanczos basis cap"), "{error}");
    let fit = JointFit {
        parameters: vec![
            JointParameter { parameter: young, lower_si: 800.0, upper_si: 1200.0, scale_si: 1000.0 },
            JointParameter {
                parameter: Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Density },
                lower_si: 1.0,
                upper_si: 3.0,
                scale_si: 2.0,
            },
        ],
        observations: vec![
            JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: 0 }, value_si: 1.0, standard_deviation_si: 0.01 },
            JointObservation {
                response: JointResponse::Static { step: "S".into(), node: 49, dof: 3 },
                value_si: 0.001,
                standard_deviation_si: 0.00001,
            },
        ],
        regularization: 0.0,
        max_iterations: 32,
        relative_increment: 1e-4,
    };
    let error = inverse_joint_uniform_material(&m, &fit, &limits).unwrap_err().to_string();
    assert!(error.contains("Lanczos basis cap"), "{error}");
}

#[test]
fn composite_sections_refuse_uniform_single_material_scaling() {
    let mut m = strip();
    let mut section = m.section("strip").unwrap().clone();
    section.material.clear();
    section.kind = SectionKind::CompositeShell {
        layup: Layup { offset: 0.0, plies: vec![Ply { material: "Mat".into(), thickness: 0.01, angle: 0.0, points: 3 }] },
    };
    m.replace_section("strip", section).unwrap();
    let young = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young };
    let scalar = FrequencyFit {
        lower: 800.0,
        upper: 1200.0,
        observed: vec![1.0],
        mode_indices: None,
        standard_deviation: None,
        regularization: 0.0,
        max_iterations: 2,
        relative_increment: 1e-4,
    };
    let err = inverse_frequency_uniform_young(&m, Some("F"), &young, &scalar, &Limits::default()).unwrap_err();
    assert!(matches!(err, fieldcraft_solver::Error::NotSupported(_)));
    let fit = JointFit {
        parameters: vec![
            JointParameter { parameter: young, lower_si: 800.0, upper_si: 1200.0, scale_si: 1000.0 },
            JointParameter {
                parameter: Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Density },
                lower_si: 1.0,
                upper_si: 3.0,
                scale_si: 2.0,
            },
        ],
        observations: vec![],
        regularization: 0.0,
        max_iterations: 2,
        relative_increment: 1e-4,
    };
    let err = inverse_joint_uniform_material(&m, &fit, &Limits::default()).unwrap_err();
    assert!(matches!(err, fieldcraft_solver::Error::NotSupported(_)));
}
