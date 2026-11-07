mod common;
use common::*;
use fieldcraft_model::design::{MaterialConstant as C, Parameter};
use fieldcraft_model::*;
use fieldcraft_solver::{
    Limits,
    sensitivity::{StaticObservation, frequency_sensitivity_with_vectors, static_sensitivity},
    solve_static,
};

#[test]
fn all_nine_orthotropic_constants_differentiate_actual_response() {
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| [p[0] * 2.0, p[1] * 0.3, p[2] * 0.2]);
    g.define(
        Elastic::EngineeringConstants { young: [1000.0, 800.0, 600.0], poisson: [0.2, 0.15, 0.1], shear: [300.0, 250.0, 200.0] },
        Some(2.0),
        full(),
        Procedure::Static {},
    );
    g.fixed("fixed", |p| near(p[0], 0.0));
    let node = g.nodes_where(|p| near(p[0], 1.0))[0];
    g.model.add_set("tip", Set { nodes: vec![node], ..Set::default() }).unwrap();
    g.load("force", Load::Force { target: Target::Set("tip".into()), vector: [1.0, 0.3, -0.2], distribution: Distribution::Total });
    for constant in [C::E1, C::E2, C::E3, C::Nu12, C::Nu13, C::Nu23, C::G12, C::G13, C::G23] {
        let parameter = Parameter::Material { name: "Mat".into(), constant };
        let value = parameter.value(&g.model).unwrap();
        let h = value * 1e-4;
        let result =
            static_sensitivity(&g.model, None, &parameter, h, StaticObservation { node, dof: 3, stress_element: Some(1) }, &Limits::default())
                .unwrap();
        let (mut plus, mut minus) = (g.model.clone(), g.model.clone());
        parameter.set(&mut plus, value + h).unwrap();
        parameter.set(&mut minus, value - h).unwrap();
        let a = solve_static(&plus, None, &Limits::default()).unwrap();
        let b = solve_static(&minus, None, &Limits::default()).unwrap();
        let i = a.nodes.binary_search(&node).unwrap();
        let fd = (a.displacement[i][2] - b.displacement[i][2]) / (2.0 * h);
        assert!((result.direct - fd).abs() < 1e-5 * fd.abs().max(1e-8), "{constant:?} {} {fd}", result.direct);
        assert!((result.direct - result.adjoint).abs() < 1e-10 * fd.abs().max(1e-8));
    }
}

fn membrane() -> Model {
    let mut m = Model::default();
    for (i, p) in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]].into_iter().enumerate() {
        m.set_node(i as u32 + 1, p).unwrap();
    }
    m.set_element(1, ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
    for (name, nodes, elements) in [("all", vec![1, 2, 3, 4], vec![1]), ("root", vec![1, 4], vec![]), ("tip", vec![2, 3], vec![])] {
        m.add_set(name, Set { nodes, elements, ..Set::default() }).unwrap();
    }
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1000.0, poisson: 0.0 }), density: Some(2.0), ..Material::default() })
        .unwrap();
    m.add_section(
        "shell",
        Section {
            set: "all".into(),
            material: "Mat".into(),
            kind: SectionKind::Shell { thickness: 0.1, offset: 0.0, points: 3 },
            formulation: Default::default(),
            orientation: None,
        },
    )
    .unwrap();
    m.add_boundary_condition("fixed", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    m.add_boundary_condition(
        "guide",
        BoundaryCondition::Displacement { target: Target::Set("all".into()), values: [None, Some(0.0), Some(0.0), Some(0.0), Some(0.0), Some(0.0)] },
    )
    .unwrap();
    m.add_load("force", Load::Force { target: Target::Set("tip".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    for (name, procedure) in [("S", Procedure::Static {}), ("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None })] {
        let mut step = Step::new(name, procedure);
        step.boundary_conditions.extend(["fixed".into(), "guide".into()]);
        if name == "S" {
            step.loads.insert("force".into());
        }
        m.add_step(step, None).unwrap();
    }
    m
}

#[test]
fn shell_thickness_changes_membrane_stiffness_mass_and_mode_normalization() {
    let m = membrane();
    let p = Parameter::ShellThickness { name: "shell".into() };
    for h in [1e-3, 1e-4, 1e-5] {
        let s = static_sensitivity(&m, Some("S"), &p, h, StaticObservation { node: 2, dof: 1, stress_element: Some(1) }, &Limits::default()).unwrap();
        assert!((s.displacement - 0.01).abs() < 1e-12);
        assert!((s.direct + 0.1).abs() < 1e-10);
        assert_eq!(s.stress.len(), 12);
        for (stress, derivative) in s.stress.iter().zip(&s.stress_derivative) {
            assert!((stress[0] - 10.0).abs() < 1e-9);
            assert!((derivative[0] + 100.0).abs() < 1e-7);
        }

        let f = frequency_sensitivity_with_vectors(&m, Some("F"), &p, h, &Default::default(), true, &Limits::default()).unwrap();
        let mode = &f.modes[0];
        assert!(mode.frequency_derivative.abs() < 1e-8);
        let v = mode.vector.as_ref().unwrap();
        for (shape, derivative) in v.shape_si.iter().flatten().zip(v.derivative_si.iter().flatten()) {
            assert!((derivative + 5.0 * shape).abs() < 1e-8);
        }
    }
}
