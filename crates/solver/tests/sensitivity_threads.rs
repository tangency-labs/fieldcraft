//! Same-platform 1/8-thread bit identity for the new physical design paths.
mod common;
use common::*;
use fieldcraft_model::{
    ElementKind, Procedure, Step,
    design::{MaterialConstant, Parameter},
    harmonic::{Damping, Frequencies, Harmonic, Method},
};
use fieldcraft_solver::{Limits, sensitivity::*};

#[test]
fn physical_design_outputs_are_identical_across_thread_counts() {
    let mut g = grid([4, 2, 2], ElementKind::Hex8, |p| [p[0] * 2.0, p[1] * 0.3, p[2] * 0.2]);
    g.define(iso(1000.0, 0.25), Some(2.0), full(), Procedure::Static {});
    g.fixed("fixed", |p| near(p[0], 0.0));
    g.face_set("tip", |p| near(p[0], 1.0));
    g.load("traction", fieldcraft_model::Load::Pressure { set: "tip".into(), magnitude: -1.0 });
    let mut f = Step::new("F", Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
    f.boundary_conditions = g.model.step("Step").unwrap().boundary_conditions.clone();
    g.model.add_step(f, None).unwrap();
    let mut h = Step::new(
        "H",
        Procedure::Harmonic(Harmonic {
            frequencies: Frequencies::List(vec![0.1, 1.0]),
            method: Method::Direct,
            damping: Damping { beta: 0.01, ..Damping::default() },
            load_factors: Default::default(),
        }),
    );
    h.boundary_conditions = g.model.step("Step").unwrap().boundary_conditions.clone();
    h.loads.insert("traction".into());
    g.model.add_step(h, None).unwrap();
    let node = g.nodes_where(|p| near(p[0], 1.0))[0];
    let parameter = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::Young };
    let run = |threads| {
        fieldcraft_parallel::with_threads(threads, || {
            let mut values = Vec::new();
            let s = static_sensitivity(
                &g.model,
                Some("Step"),
                &parameter,
                0.01,
                StaticObservation { node, dof: 1, stress_element: Some(1) },
                &Limits::default(),
            )
            .unwrap();
            values.extend([s.displacement, s.direct, s.adjoint, s.load_work_derivative, s.primal_relative_residual, s.tangent_relative_residual]);
            values.extend(s.stress.into_iter().flatten());
            values.extend(s.stress_derivative.into_iter().flatten());
            let f = frequency_sensitivity_with_vectors(&g.model, Some("F"), &parameter, 0.01, &Default::default(), true, &Limits::default()).unwrap();
            let observed = f.modes.iter().map(|m| m.frequency * 1.1_f64.sqrt()).collect();
            for m in f.modes {
                values.extend([m.frequency, m.eigenvalue, m.frequency_derivative, m.eigenvalue_derivative, m.modal_mass]);
                let v = m.vector.unwrap();
                values.extend(v.shape_si.into_iter().flatten());
                values.extend(v.derivative_si.into_iter().flatten());
            }
            let h = harmonic_sensitivity(&g.model, Some("H"), &parameter, 0.01, HarmonicObservation { node, dof: 1 }, &Limits::default()).unwrap();
            let harmonic_observed = h.frames.iter().flat_map(|f| f.displacement).collect();
            for f in h.frames {
                values.extend(f.displacement);
                values.extend(f.direct);
                values.extend(f.adjoint);
                values.extend([f.amplitude, f.phase, f.amplitude_derivative, f.phase_derivative]);
            }
            let fit = FrequencyFit {
                lower: 800.0,
                upper: 1300.0,
                observed,
                standard_deviation: None,
                regularization: 0.0,
                max_iterations: 16,
                relative_increment: 1e-5,
            };
            let estimate = inverse_frequency_uniform_young(&g.model, Some("F"), &parameter, &fit, &Limits::default()).unwrap();
            values.push(estimate.estimate);
            values.extend(estimate.predicted);
            values.extend(estimate.residuals);
            let harmonic_fit = FrequencyFit { observed: harmonic_observed, ..fit };
            let estimate =
                inverse_harmonic(&g.model, Some("H"), &parameter, &harmonic_fit, HarmonicObservation { node, dof: 1 }, &Limits::default()).unwrap();
            values.push(estimate.estimate);
            values.extend(estimate.predicted);
            values.into_iter().map(f64::to_bits).collect::<Vec<_>>()
        })
    };
    let one = run(1);
    let eight = run(8);
    assert_eq!(one, eight);
    assert!(one.len() > 500);
}
