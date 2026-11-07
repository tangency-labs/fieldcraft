//! Whole-model checks against a consistent-mass axial bar's exact discrete equations.
use fieldcraft_model::harmonic::{Damping, Frequencies, Harmonic, LoadFactor, Method};
use fieldcraft_model::*;
use fieldcraft_solver::harmonic::solve_harmonic;
use fieldcraft_solver::{Limits, solve_static};
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
            orientation: None,
            set: "all".into(),
            material: "bar".into(),
            kind: SectionKind::Rod { area: 1.0 },
            formulation: Default::default(),
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
fn material_and_step_damping_phase_stress_and_prescribed_motion() {
    let d = Damping { alpha: 0.02, beta: 0.005, loss_factor: 0.03 };
    for prescribed in [0.0, 0.1] {
        let material = rod(d, Damping::default(), Method::Direct, prescribed);
        let step = rod(Damping::default(), d, Method::Direct, prescribed);
        let a = solve_harmonic(&material, None, &Limits::default()).unwrap();
        let b = solve_harmonic(&step, None, &Limits::default()).unwrap();
        for (a, b) in a.frames.iter().zip(&b.frames) {
            let w = 2.0 * PI * a.frequency;
            let (x, y) = (4.0 - w * w, w * (d.alpha + 4.0 * d.beta) + 4.0 * d.loss_factor);
            let fr = 2.0 * 0.4_f64.cos() - prescribed * (-4.0 - 0.5 * w * w);
            let fi = 2.0 * 0.4_f64.sin() - prescribed * (w * (0.5 * d.alpha - 4.0 * d.beta) - 4.0 * d.loss_factor);
            let re = (fr * x + fi * y) / (x * x + y * y);
            let im = (fi * x - fr * y) / (x * x + y * y);
            assert!((a.real[1][0] - re).abs() < 1e-11);
            assert!((a.imaginary[1][0] - im).abs() < 1e-11);
            assert!((a.real[1][0] - b.real[1][0]).abs() < 1e-11);
            assert!((a.imaginary[1][0] - b.imaginary[1][0]).abs() < 1e-11);
            assert!((a.stress_real.as_ref().unwrap()[1][0] - 4.0 * (re - prescribed)).abs() < 1e-10);
            assert!((a.stress_imaginary.as_ref().unwrap()[1][0] - 4.0 * im).abs() < 1e-10);
            assert!((a.input_power - a.dissipated_power).abs() < 1e-10);
        }
    }
}

#[test]
fn residual_vectors_refuse_nonzero_prescribed_motion() {
    let d = Damping { alpha: 0.02, beta: 0.005, loss_factor: 0.0 };
    let method = Method::Modal { modes: 1, damping_ratios: vec![], residual_vectors: 1 };
    let m = rod(d, Damping::default(), method, 0.1);
    let error = solve_harmonic(&m, None, &Limits::default()).unwrap_err().to_string();
    assert!(error.contains("residual vectors") && error.contains("prescribed motion"), "{error}");
}

#[test]
fn zero_frequency_static_limit_and_full_modal_agreement() {
    let d = Damping { alpha: 0.02, beta: 0.005, loss_factor: 0.0 };
    let direct = rod(d, Damping::default(), Method::Direct, 0.0);
    let modal = rod(d, Damping::default(), Method::Modal { modes: 1, damping_ratios: vec![], residual_vectors: 0 }, 0.0);
    let a = solve_harmonic(&direct, None, &Limits::default()).unwrap();
    let b = solve_harmonic(&modal, None, &Limits::default()).unwrap();
    for (a, b) in a.frames.iter().zip(&b.frames) {
        assert!((a.real[1][0] - b.real[1][0]).abs() < 1e-10);
        assert!((a.imaginary[1][0] - b.imaginary[1][0]).abs() < 1e-10);
        assert!(b.relative_residual < 1e-10);
    }
    let mut m = direct.clone();
    let mut s = Step::new("S", Procedure::Static {});
    s.boundary_conditions.extend(["root".into(), "guide".into()]);
    s.loads.insert("force".into());
    m.add_step(s, None).unwrap();
    let s = solve_static(&m, Some("S"), &Limits::default()).unwrap();
    assert!((a.frames[0].real[1][0] - 2.0 * 0.4_f64.cos() * s.displacement[1][0]).abs() < 1e-13);
}
