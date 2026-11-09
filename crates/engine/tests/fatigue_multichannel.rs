//! Independent tensor superposition checks; actual static command wiring is tested in fatigue.rs.
use fieldcraft_solver::{
    StaticResult,
    fatigue::{Basquin, history::evaluate, life::StressMean, stress::Criterion},
};
use serde_json::json;
fn reference(stress: [f64; 6]) -> StaticResult {
    serde_json::from_value(json!({"step":"Reference","nodes":[1],"displacement":[[0,0,0]],"rotation":[],"reaction":[[0,0,0]],"reaction_moment":[],
        "stress":[stress],"von_mises":[0],"elements":[],"shells":[],"beams":[],"reference_points":[],
        "applied_total":[0,0,0],"reaction_total":[0,0,0],"stats":fieldcraft_solver::Stats::default(),"warnings":[],"temperature":null,"nonlinear":null,"contacts":[]})).unwrap()
}
#[test]
fn tensor_channels_sum_before_nonlinear_scalarization_and_cancel_exactly() {
    let a = reference([3., 0., 0., 0., 0., 0.]);
    let b = reference([0., 0., 0., 2., 0., 0.]);
    let curve = Basquin { amplitude_ref: 1., cycles_ref: 100., exponent: 2., endurance_limit: 0. };
    let coeff = vec![vec![-1., 1., -1.], vec![-1., 1., -1.]];
    // Plane-stress eigenvalues of [[3,2],[2,0]] are 4 and -1; VM^2=21.
    // Summing individual equivalent stresses would produce (3+sqrt(12))^2 instead.
    for (criterion, expected) in [(Criterion::SignedVonMises, 0.21), (Criterion::MaxPrincipal, 0.0225)] {
        let r = evaluate(&[&a, &b], &coeff, &[0., 1., 2.], &[1], criterion, &curve, StressMean::None).unwrap();
        assert!((r[0].damage - expected).abs() < 1e-13);
        assert_eq!(r[0].cycle_count, 1.);
    }
    let b = reference([-3., 0., 0., 0., 0., 0.]);
    let r = evaluate(&[&a, &b], &coeff, &[0., 1., 2.], &[1], Criterion::SignedVonMises, &curve, StressMean::None).unwrap();
    assert_eq!(r[0].damage, 0.);
}
#[test]
fn combined_nonproportional_history_and_sign_ties_are_named_refusals() {
    let a = reference([3., 0., 0., 0., 0., 0.]);
    let b = reference([0., 0., 0., 2., 0., 0.]);
    let curve = Basquin { amplitude_ref: 1., cycles_ref: 100., exponent: 2., endurance_limit: 0. };
    for criterion in [Criterion::SignedVonMises, Criterion::MaxPrincipal] {
        let e = evaluate(&[&a, &b], &[vec![-1., 1., -1.], vec![0., 1., 0.]], &[0., 1., 2.], &[1], criterion, &curve, StressMean::None)
            .unwrap_err()
            .to_string();
        assert!(e.contains("nonproportional"), "{e}");
    }
    let e = evaluate(&[&a, &b], &[vec![0.; 3], vec![-1., 1., -1.]], &[0., 1., 2.], &[1], Criterion::SignedVonMises, &curve, StressMean::None)
        .unwrap_err()
        .to_string();
    assert!(e.contains("ambiguous"), "{e}");
}
