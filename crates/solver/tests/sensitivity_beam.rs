//! Euler–Bernoulli cantilever strip: u=FL³/(3EI), f₁=β₁² sqrt(EI/rho A)/(2πL²).
//! A thin shell strip approaches this theory; no external solver values are used.
use fieldcraft_model::{
    design::{MaterialConstant, Parameter},
    *,
};
use fieldcraft_solver::{
    Limits,
    sensitivity::{StaticObservation, frequency_sensitivity, static_sensitivity},
};

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
fn cantilever_material_and_thickness_derivatives_match_beam_theory() {
    let model = strip();
    let (length, width, thickness, young, density, force) = (2.0_f64, 0.1, 0.01, 1000.0, 2.0, 1e-9);
    let inertia = width * thickness * thickness * thickness / 12.0;
    let u = force * length * length * length / (3.0 * young * inertia);
    let beta = 1.875104068711961_f64;
    let frequency = beta * beta / (2.0 * std::f64::consts::PI * length * length) * (young * inertia / (density * width * thickness)).sqrt();
    for (parameter, value, du, df) in [
        (Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young }, young, -u / young, frequency / (2.0 * young)),
        (Parameter::ShellThickness { name: "strip".into() }, thickness, -3.0 * u / thickness, frequency / thickness),
    ] {
        for fraction in [1e-3, 1e-4, 1e-5] {
            let s = static_sensitivity(
                &model,
                Some("S"),
                &parameter,
                value * fraction,
                StaticObservation { node: 49, dof: 3, stress_element: None },
                &Limits::default(),
            )
            .unwrap();
            assert!((s.displacement / u - 1.0).abs() < 0.003, "{} {u}", s.displacement);
            assert!((s.direct / du - 1.0).abs() < 0.003, "{} {du}", s.direct);
            assert!((s.adjoint / s.direct - 1.0).abs() < 1e-8);
            let f = frequency_sensitivity(&model, Some("F"), &parameter, value * fraction, &Default::default(), &Limits::default()).unwrap();
            assert!((f.modes[0].frequency / frequency - 1.0).abs() < 0.005, "{} {frequency}", f.modes[0].frequency);
            assert!((f.modes[0].frequency_derivative / df - 1.0).abs() < 0.005, "{} {df}", f.modes[0].frequency_derivative);
        }
    }
}
