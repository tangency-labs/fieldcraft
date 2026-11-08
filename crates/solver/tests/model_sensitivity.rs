mod common;
use common::*;
use fieldcraft_model::{
    ElementKind, Load, Procedure,
    design::{MaterialConstant, Parameter},
};
use fieldcraft_solver::{
    Limits,
    sensitivity::{StaticObservation, static_sensitivity},
    solve_static,
};

fn block() -> Grid {
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(1000.0, 0.0), Some(2.0), full(), Procedure::Static {});
    g.fixed("fixed", |p| near(p[0], 0.0));
    g.face_set("right", |p| near(p[0], 1.0));
    g.load("traction", Load::Pressure { set: "right".into(), magnitude: -1.0 });
    g
}
fn young() -> Parameter {
    Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young }
}

#[test]
fn traction_bar_direct_adjoint_stress_and_work_match_theory() {
    let g = block();
    let node = g.nodes_where(|p| near(p[0], 1.0))[0];
    let query = StaticObservation { node, dof: 1, stress_element: Some(1) };
    for h in [1.0, 0.1, 0.01] {
        let r = static_sensitivity(&g.model, None, &young(), h, query, &Limits::default()).unwrap();
        assert!((r.displacement - 0.001).abs() < 1e-13);
        assert!((r.direct + 1e-6).abs() < 1e-13);
        assert!((r.adjoint - r.direct).abs() < 1e-16);
        assert!((r.load_work_derivative + 1e-6).abs() < 1e-13);
        assert!(r.stress_derivative.iter().flatten().all(|v| v.abs() < 1e-10));
        assert!(r.primal_relative_residual < 1e-12);
        assert!(r.tangent_relative_residual < 1e-12);
    }
}

#[test]
fn poisson_derivative_at_zero_and_density_load_are_not_lost() {
    let mut g = block();
    let node = g.nodes_where(|p| near(p[0], 1.0))[0];
    let p = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Poisson };
    let observation = StaticObservation { node, dof: 2, stress_element: None };
    let r = static_sensitivity(&g.model, None, &p, 1e-5, observation, &Limits::default()).unwrap();
    let mut plus = g.model.clone();
    let mut minus = g.model.clone();
    p.set(&mut plus, 1e-5).unwrap();
    p.set(&mut minus, -1e-5).unwrap();
    let a = solve_static(&plus, None, &Limits::default()).unwrap();
    let b = solve_static(&minus, None, &Limits::default()).unwrap();
    let i = a.nodes.binary_search(&node).unwrap();
    assert!((r.direct - (a.displacement[i][1] - b.displacement[i][1]) / 2e-5).abs() < 1e-10);
    g.load("gravity", Load::Gravity { set: Some("all".into()), acceleration: [1.0, 0.0, 0.0] });
    let rho = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Density };
    let r = static_sensitivity(&g.model, None, &rho, 1e-3, StaticObservation { dof: 1, ..observation }, &Limits::default()).unwrap();
    assert!((r.direct - 0.0005).abs() < 1e-11);
    assert!((r.direct - r.adjoint).abs() < 1e-15);
}

#[test]
fn parameters_refuse_bad_values_without_mutation() {
    let mut g = block();
    let original = g.model.clone();
    assert!(young().set(&mut g.model, -1.0).is_err());
    assert_eq!(g.model, original);
    let p = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Poisson };
    assert!(p.set(&mut g.model, 0.5).is_err());
    assert_eq!(g.model, original);
    assert!(Parameter::Material { name: "Mat".into(), constant: MaterialConstant::E1 }.value(&g.model).is_err());
    assert!(
        static_sensitivity(&g.model, None, &young(), 0.0, StaticObservation { node: 1, dof: 1, stress_element: None }, &Limits::default()).is_err()
    );
}

#[test]
fn actual_fe_modal_derivatives_follow_young_and_density_scaling() {
    let mut g = grid([3, 1, 1], ElementKind::Hex8, |p| [p[0] * 2.0, p[1] * 0.3, p[2] * 0.2]);
    g.define(iso(1000.0, 0.25), Some(2.0), full(), Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
    g.fixed("fixed", |p| near(p[0], 0.0));
    for (p, h, sign, value) in
        [(young(), 1.0, 1.0, 1000.0), (Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Density }, 0.01, -1.0, 2.0)]
    {
        let result = fieldcraft_solver::sensitivity::frequency_sensitivity(&g.model, None, &p, h, &Default::default(), &Limits::default()).unwrap();
        assert_eq!(result.modes.len(), 3);
        for mode in result.modes {
            assert!((mode.frequency_derivative / (sign * mode.frequency / (2.0 * value)) - 1.0).abs() < 1e-8, "{p:?} {mode:?}");
        }
    }
}

#[test]
fn prescribed_motion_is_differentiated_in_the_reduced_right_hand_side() {
    let mut g = block();
    g.model
        .replace_boundary_condition(
            "fixed",
            fieldcraft_model::BoundaryCondition::Displacement {
                target: fieldcraft_model::Target::Set("fixed".into()),
                values: [Some(0.1), Some(0.0), Some(0.0), None, None, None],
            },
        )
        .unwrap();
    let node = g.nodes_where(|p| near(p[0], 1.0))[0];
    let r =
        static_sensitivity(&g.model, None, &young(), 0.01, StaticObservation { node, dof: 1, stress_element: Some(1) }, &Limits::default()).unwrap();
    assert!((r.displacement - 0.101).abs() < 1e-12);
    assert!((r.direct + 1e-6).abs() < 1e-11);
    assert!(r.stress_derivative.iter().flatten().all(|v| v.abs() < 1e-8));
}

#[test]
fn oriented_solid_sensitivity_recovers_global_stress() {
    use fieldcraft_model::{Elastic, MaterialOrientation, SectionOrientation};
    let mut g = block();
    let mut mat = g.model.material("Mat").unwrap().clone();
    mat.elastic = Some(Elastic::EngineeringConstants { young: [1000.0, 100.0, 200.0], poisson: [0.0; 3], shear: [50.0, 40.0, 30.0] });
    g.model.replace_material("Mat", mat).unwrap();
    let mut section = g.model.section("solid").unwrap().clone();
    section.orientation = Some(SectionOrientation {
        frame: MaterialOrientation::Rectangular { axis1: [0.0, 1.0, 0.0], axis2: [-1.0, 0.0, 0.0] },
        shell_reference_axis: 1,
    });
    g.model.replace_section("solid", section).unwrap();
    let p = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::E2 };
    let node = g.nodes_where(|p| near(p[0], 1.0))[0];
    let r = static_sensitivity(&g.model, None, &p, 0.01, StaticObservation { node, dof: 1, stress_element: Some(1) }, &Limits::default()).unwrap();
    assert!((r.displacement - 0.01).abs() < 1e-12);
    assert!((r.direct + 0.0001).abs() < 1e-12);
    assert!((r.direct - r.adjoint).abs() < 1e-14);
    for point in r.stress {
        assert!((point[0] - 1.0).abs() < 1e-12);
    }
    assert!(r.stress_derivative.iter().flatten().all(|v| v.abs() < 1e-9));
}

#[test]
fn nelson_mass_normalization_derivatives_for_uniform_material_changes() {
    use fieldcraft_solver::sensitivity::frequency_sensitivity_with_vectors;
    let mut g = grid([3, 1, 1], ElementKind::Hex8, |p| [p[0] * 2.0, p[1] * 0.3, p[2] * 0.2]);
    g.define(iso(1000.0, 0.25), Some(2.0), full(), Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
    g.fixed("fixed", |p| near(p[0], 0.0));
    for (parameter, h, factor) in
        [(young(), 0.1, 0.0), (Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Density }, 1e-3, -0.25)]
    {
        let result = frequency_sensitivity_with_vectors(&g.model, None, &parameter, h, &Default::default(), true, &Limits::default()).unwrap();
        for mode in result.modes {
            let v = mode.vector.unwrap();
            let scale = v.shape_si.iter().flatten().map(|x| x.abs()).fold(0.0_f64, f64::max);
            for (shape, derivative) in v.shape_si.iter().flatten().zip(v.derivative_si.iter().flatten()) {
                assert!((derivative - factor * shape).abs() < 1e-7 * scale, "{derivative} {shape}");
            }
            assert!(v.equation_relative_residual < 1e-7);
            assert!(v.normalization_derivative_residual.abs() < 1e-10);
        }
    }
}

#[test]
fn exact_uniform_young_fit_agrees_with_recomputed_finite_elements() {
    use fieldcraft_solver::sensitivity::{FrequencyFit, frequency_sensitivity, inverse_frequency, inverse_frequency_uniform_young};
    let mut g = grid([3, 1, 1], ElementKind::Hex8, |p| [p[0] * 2.0, p[1] * 0.3, p[2] * 0.2]);
    g.define(iso(1000.0, 0.25), Some(2.0), full(), Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
    g.fixed("fixed", |p| near(p[0], 0.0));
    let mut truth = g.model.clone();
    young().set(&mut truth, 1100.0).unwrap();
    let target = frequency_sensitivity(&truth, None, &young(), 0.01, &Default::default(), &Limits::default()).unwrap();
    let fit = FrequencyFit {
        lower: 800.0,
        upper: 1300.0,
        observed: target.modes.iter().map(|m| m.frequency).collect(),
        mode_indices: None,
        standard_deviation: None,
        regularization: 0.0,
        max_iterations: 16,
        relative_increment: 1e-4,
    };
    let fast = inverse_frequency_uniform_young(&g.model, None, &young(), &fit, &Limits::default()).unwrap();
    let general = inverse_frequency(&g.model, None, &young(), &fit, &Limits::default()).unwrap();
    assert!((fast.estimate / 1100.0 - 1.0).abs() < 1e-6);
    assert!((fast.estimate / general.estimate - 1.0).abs() < 1e-7);
    assert_eq!(fast.model_evaluations, 1);
    assert!(fast.baseline.is_some());
    assert_eq!(general.model_evaluations, general.evaluations);
    let rho = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Density };
    assert!(inverse_frequency_uniform_young(&g.model, None, &rho, &fit, &Limits::default()).is_err());
    let mut material = g.model.material("Mat").unwrap().clone();
    material.elastic = Some(fieldcraft_model::Elastic::EngineeringConstants { young: [1000.0, 900.0, 800.0], poisson: [0.2; 3], shear: [300.0; 3] });
    g.model.replace_material("Mat", material).unwrap();
    assert!(inverse_frequency_uniform_young(&g.model, None, &young(), &fit, &Limits::default()).is_err());
}
