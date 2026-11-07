//! Harmonic matrix systems against independent oscillator equations and full solves.
use fieldcraft_model::harmonic::Damping;
use fieldcraft_solver::Limits;
use fieldcraft_solver::complex::Vector;
use fieldcraft_solver::harmonic::{Matrices, Modal};
use fieldcraft_solver::sparse::Csr;

fn dense(values: &[f64], n: usize) -> Csr {
    Csr { n, ptr: (0..=n).map(|i| i * n).collect(), idx: (0..n).flat_map(|_| 0..n).collect(), val: values.to_vec() }
}

fn assert_close(a: &Vector, b: &Vector) {
    let scale = b.real.iter().chain(&b.imaginary).map(|v| v.abs()).fold(1.0_f64, f64::max);
    for (a, b) in a.real.iter().chain(&a.imaginary).zip(b.real.iter().chain(&b.imaginary)) {
        assert!((a - b).abs() < 2e-10 * scale, "{a} != {b}");
    }
}

#[test]
fn all_modes_match_direct_for_nonproportional_viscous_and_loss_damping() {
    let k = dense(&[2.0, -1.0, -1.0, 2.0], 2);
    let m = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let c = dense(&[0.15, -0.02, -0.02, 0.06], 2);
    let g = dense(&[0.03, 0.01, 0.01, 0.05], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: Some(&c), loss: Some(&g) };
    let limits = Limits::default();
    let modal = Modal::new(&matrices, 2, &limits).unwrap();
    let force = Vector { real: vec![1.0, 0.2], imaginary: vec![-0.3, 0.8] };
    let damping = Damping { alpha: 0.02, beta: 0.01, loss_factor: 0.03 };
    for omega in [0.0, 0.99, 1.0, 1.01, 3.0_f64.sqrt(), 2.0] {
        let direct = matrices.direct(omega, damping, &force, &limits).unwrap();
        let reduced = modal.solve(omega, damping, &[], &force, &limits).unwrap();
        assert_close(&reduced.displacement, &direct.displacement);
        assert!(!reduced.diagonal);
        assert!(reduced.full_relative_residual < 2e-12);
        assert!(reduced.projected_residual.relative < 2e-12);
    }
}

#[test]
fn proportional_and_modal_damping_use_the_diagonal_and_agree() {
    let k = dense(&[4.0], 1);
    let m = dense(&[1.0], 1);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let modal = Modal::new(&matrices, 1, &limits).unwrap();
    let force = Vector { real: vec![1.0], imaginary: vec![0.0] };
    // alpha=0.04 + beta=0.02 * omega_n^2=0.08 => c=0.12 => zeta=0.03.
    for omega in [0.0, 1.9, 2.0, 2.1, 5.0] {
        let rayleigh = matrices.direct(omega, Damping { alpha: 0.04, beta: 0.02, loss_factor: 0.0 }, &force, &limits).unwrap();
        let ratio = modal.solve(omega, Damping::default(), &[0.03], &force, &limits).unwrap();
        assert!(ratio.diagonal);
        assert_close(&ratio.displacement, &rayleigh.displacement);
        assert!(ratio.full_relative_residual < 1e-13);
    }
}

#[test]
fn truncated_modes_report_the_unprojected_residual() {
    let k = dense(&[2.0, -1.0, -1.0, 2.0], 2);
    let m = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let modal = Modal::new(&matrices, 1, &limits).unwrap();
    let force = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.0] };
    let out = modal.solve(0.5, Damping { alpha: 0.04, ..Damping::default() }, &[], &force, &limits).unwrap();
    assert!(out.projected_residual.relative < 1e-13);
    assert!((out.full_relative_residual - 0.5).abs() < 1e-12);
}

#[test]
fn modal_basis_and_response_are_identical_across_thread_counts() {
    let k = dense(&[3.0, -1.0, -1.0, 2.0], 2);
    let m = dense(&[2.0, 0.2, 0.2, 1.0], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let force = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.5] };
    let run = || {
        let modes = Modal::new(&matrices, 2, &limits).unwrap();
        modes.solve(1.0, Damping { alpha: 0.05, beta: 0.01, loss_factor: 0.03 }, &[], &force, &limits).unwrap().displacement
    };
    let a = fieldcraft_parallel::with_threads(1, run);
    for threads in [2, 4, 6] {
        let b = fieldcraft_parallel::with_threads(threads, run);
        for (a, b) in a.real.iter().chain(&a.imaginary).zip(b.real.iter().chain(&b.imaginary)) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}

#[test]
fn free_body_zero_frequency_is_not_regularized_by_the_modal_extraction_shift() {
    // Two unit masses joined by a unit spring: [1,1] is an exact rigid mode.
    // A uniform nonzero static force cannot be balanced; at positive omega the
    // same mode has the exact inertial response -F / omega².
    let k = dense(&[1.0, -1.0, -1.0, 1.0], 2);
    let m = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let force = Vector { real: vec![1.0, 1.0], imaginary: vec![0.0, 0.0] };
    let modal = Modal::new(&matrices, 2, &limits).unwrap();
    assert!(matrices.direct(0.0, Damping::default(), &force, &limits).is_err());
    for damping in [Damping::default(), Damping { alpha: 0.1, beta: 0.02, loss_factor: 0.03 }] {
        assert!(modal.solve(0.0, damping, &[0.02], &force, &limits).is_err());
    }
    let response = modal.solve(0.5, Damping::default(), &[], &force, &limits).unwrap();
    assert_close(&response.displacement, &Vector { real: vec![-4.0, -4.0], imaginary: vec![0.0, 0.0] });
    assert!(response.full_relative_residual < 1e-12);
}

#[test]
fn undamped_modal_poles_are_refused_without_shifting_nearby_frequencies() {
    // Eigenvalues are exactly 1 and 3; the extracted values may differ by a ULP.
    let k = dense(&[2.0, -1.0, -1.0, 2.0], 2);
    let m = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let force = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.0] };
    let modal = Modal::new(&matrices, 2, &limits).unwrap();
    for omega in [1.0, 3.0_f64.sqrt()] {
        assert!(matrices.direct(omega, Damping::default(), &force, &limits).is_err());
        assert!(modal.solve(omega, Damping::default(), &[], &force, &limits).is_err());
    }
    for omega in [1.0 - 1e-5, 1.0 + 1e-5] {
        let direct = matrices.direct(omega, Damping::default(), &force, &limits).unwrap();
        let reduced = modal.solve(omega, Damping::default(), &[], &force, &limits).unwrap();
        assert_close(&reduced.displacement, &direct.displacement);
        assert!(reduced.full_relative_residual < 1e-8);
    }
}
