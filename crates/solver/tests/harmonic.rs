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
    // The missing modal force is [1/2,-1/2], independent of frequency for
    // proportional damping. This separates truncation from linear solve error.
    for omega in [0.0, 0.5, 0.99, 1.01, 2.0] {
        let out = modal.solve(omega, Damping { alpha: 0.04, beta: 0.01, loss_factor: 0.02 }, &[], &force, &limits).unwrap();
        assert!(out.projected_residual.relative < 1e-13);
        assert!((out.full_relative_residual - 0.5).abs() < 1e-12);
    }
}

#[test]
fn residual_vectors_recover_omitted_response_with_full_damping_projection() {
    let k = dense(&[2.0, -1.0, -1.0, 2.0], 2);
    let m = dense(&[2.0, 0.2, 0.2, 1.0], 2);
    let c = dense(&[0.15, -0.02, -0.02, 0.06], 2);
    let g = dense(&[0.03, 0.01, 0.01, 0.05], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: Some(&c), loss: Some(&g) };
    let limits = Limits::default();
    let force = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.5] };
    let run = || {
        let modal = Modal::with_residual_vectors(&matrices, 1, &force, 1, &limits).unwrap();
        assert_eq!(modal.eigenvalues().len(), 1, "Ritz vectors are not reported as natural modes");
        assert_eq!(modal.residual_vectors(), 1);
        let damping = Damping { alpha: 0.04, beta: 0.01, loss_factor: 0.02 };
        let mut response = Vec::new();
        for omega in [0.0, 0.5, 0.99, 1.01, 2.0] {
            let full = matrices.direct(omega, damping, &force, &limits).unwrap();
            let corrected = modal.solve(omega, damping, &[], &force, &limits).unwrap();
            assert_close(&corrected.displacement, &full.displacement);
            assert!(corrected.full_relative_residual < 1e-12);
            // A modal ratio acts only on the retained physical mode, not the Ritz vector.
            let ratio = 0.03;
            let phi = &modal.shapes()[0];
            let mut mphi = vec![0.0; 2];
            m.mul(phi, &mut mphi);
            let mut with_ratio = c.clone();
            for i in 0..2 {
                for j in 0..2 {
                    with_ratio.val[i * 2 + j] += 2.0 * ratio * modal.eigenvalues()[0].sqrt() * mphi[i] * mphi[j];
                }
            }
            let full = Matrices { viscous: Some(&with_ratio), ..matrices }.direct(omega, damping, &force, &limits).unwrap();
            let corrected = modal.solve(omega, damping, &[ratio], &force, &limits).unwrap();
            assert_close(&corrected.displacement, &full.displacement);
            assert!(corrected.full_relative_residual < 1e-12);
            response.extend(corrected.displacement.real.into_iter().chain(corrected.displacement.imaginary).map(f64::to_bits));
        }
        response
    };
    let reference = fieldcraft_parallel::with_threads(1, run);
    for threads in [2, 4, 6] {
        assert_eq!(reference, fieldcraft_parallel::with_threads(threads, run));
    }
}

#[test]
fn residual_krylov_chain_converges_toward_independent_oscillator_solutions() {
    let n = 12;
    let mut k = dense(&vec![0.0; n * n], n);
    let mut m = k.clone();
    for j in 0..n {
        m.val[j * n + j] = 1.0 + j as f64 / 10.0;
        k.val[j * n + j] = m.val[j * n + j] * (j + 1).pow(2) as f64;
    }
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let force = Vector { real: vec![1.0; n], imaginary: vec![0.0; n] };
    let damping = Damping { alpha: 0.01, beta: 0.02, loss_factor: 0.03 };
    let plain = Modal::new(&matrices, 1, &limits).unwrap().solve(0.2, damping, &[], &force, &limits).unwrap();
    assert!(plain.full_relative_residual > 0.9);
    let corrected = Modal::with_residual_vectors(&matrices, 1, &force, 6, &limits).unwrap();
    assert!(corrected.residual_vectors() <= 6);
    assert!(1 + corrected.residual_vectors() < n, "test must retain actual truncation");
    for omega in [0.0, 0.1, 0.2] {
        let out = corrected.solve(omega, damping, &[], &force, &limits).unwrap();
        assert!(out.full_relative_residual < 1e-8, "{}", out.full_relative_residual);
        for j in 0..n {
            let kj = k.val[j * n + j];
            let mj = m.val[j * n + j];
            let a = kj - omega * omega * mj;
            let b = omega * (damping.alpha * mj + damping.beta * kj) + damping.loss_factor * kj;
            let re = a / (a * a + b * b);
            let im = -b / (a * a + b * b);
            assert!((out.displacement.real[j] - re).abs() < 1e-10);
            assert!((out.displacement.imaginary[j] - im).abs() < 1e-10);
        }
    }
    let bad = Vector { real: vec![f64::NAN; n], imaginary: vec![0.0; n] };
    assert!(Modal::with_residual_vectors(&matrices, 1, &bad, 1, &limits).is_err());
    assert!(Modal::with_residual_vectors(&matrices, 1, &force, 17, &limits).is_err());
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
    assert!(Modal::with_residual_vectors(&matrices, 1, &force, 1, &limits).is_err());
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

#[test]
fn reusable_modal_basis_retains_per_call_limits_and_dynamic_checks() {
    let k = dense(&[4.0, 0.0, 0.0, 9.0], 2);
    let m = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let limits = Limits::default();
    let modal = Modal::new(&matrices, 1, &limits).unwrap();
    let force = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.0] };
    // The projected one-mode system fits; the two physical equations do not.
    let reduced = Limits { max_dofs: 1, ..Limits::default() };
    assert!(modal.solve(0.5, Damping::default(), &[], &force, &reduced).is_err());
    assert!(modal.solve(f64::NAN, Damping::default(), &[], &force, &limits).is_err());
    assert!(modal.solve(0.5, Damping { alpha: -1.0, ..Damping::default() }, &[], &force, &limits).is_err());
    assert!(modal.solve(f64::MAX, Damping::default(), &[], &force, &limits).is_err());
    let direct = matrices.direct(0.5, Damping::default(), &force, &limits).unwrap();
    let reused = modal.solve(0.5, Damping::default(), &[], &force, &limits).unwrap();
    assert_close(&reused.displacement, &direct.displacement);
    assert!(reused.full_relative_residual < 1e-12);
}

#[test]
fn enriched_modal_damping_keeps_exact_transpose_symmetry() {
    let n = 12;
    let mut k = dense(&vec![0.; n * n], n);
    let mut m = k.clone();
    for i in 0..n {
        m.val[i * n + i] = 1. + i as f64 / 10.;
        k.val[i * n + i] = m.val[i * n + i] * (i + 1).pow(2) as f64;
    }
    let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
    let force = Vector { real: (0..n).map(|i| 1. / (i + 1) as f64).collect(), imaginary: (0..n).map(|i| 0.1 * i as f64).collect() };
    let limits = Limits::default();
    let modal = Modal::with_residual_vectors(&matrices, 3, &force, 1, &limits).unwrap();
    assert!(modal.residual_vectors() > 0);
    let ratios = [0.01, 0.02, 0.03];
    // Independent full-space damping operator, then a Galerkin projection of
    // K, M and C together. Build each symmetric pair once, without tolerances.
    let mut c = dense(&vec![0.; n * n], n);
    for (j, phi) in modal.shapes().iter().enumerate() {
        let mut mp = vec![0.; n];
        m.mul(phi, &mut mp);
        let factor = 2. * ratios[j] * modal.eigenvalues()[j].sqrt();
        for row in 0..n {
            for col in 0..=row {
                let v = c.val[row * n + col] + factor * (mp[row] * mp[col]);
                c.val[row * n + col] = v;
                c.val[col * n + row] = v;
            }
        }
    }
    let basis = modal.basis_vectors();
    let count = basis.len();
    let project = |a: &Csr| {
        let mut out = dense(&vec![0.; count * count], count);
        let mut v = vec![0.; n];
        for j in 0..count {
            a.mul(&basis[j], &mut v);
            for (i, phi) in basis.iter().enumerate().take(j + 1) {
                let x = phi.iter().zip(&v).map(|(a, b)| a * b).sum();
                out.val[i * count + j] = x;
                out.val[j * count + i] = x;
            }
        }
        out
    };
    let (pk, pm, pc) = (project(&k), project(&m), project(&c));
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(a, b)| a * b).sum::<f64>();
    let rhs =
        Vector { real: basis.iter().map(|v| dot(v, &force.real)).collect(), imaginary: basis.iter().map(|v| dot(v, &force.imaginary)).collect() };
    for omega in [0.1, 1., 3., 10.] {
        let reference =
            Matrices { stiffness: &pk, mass: &pm, viscous: Some(&pc), loss: None }.direct(omega, Damping::default(), &rhs, &limits).unwrap();
        let mut expected = Vector::zeros(n);
        for (j, phi) in basis.iter().enumerate() {
            for (i, p) in phi.iter().enumerate() {
                expected.real[i] += p * reference.displacement.real[j];
                expected.imaginary[i] += p * reference.displacement.imaginary[j];
            }
        }
        let out = modal.solve(omega, Damping::default(), &ratios, &force, &limits).unwrap();
        assert_close(&out.displacement, &expected);
        assert!(out.projected_residual.relative < 1e-12);
    }
}
