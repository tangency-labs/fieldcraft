//! Analytic checks for the complex symmetric direct kernel. These matrices are
//! written here from oscillator equations; no external solver data is used.

use fieldcraft_solver::complex::{Vector, residual, solve};
use fieldcraft_solver::sparse::Csr;
use fieldcraft_solver::{Error, Limits};

fn dense(values: &[f64], n: usize) -> Csr {
    Csr { n, ptr: (0..=n).map(|i| i * n).collect(), idx: (0..n).flat_map(|_| 0..n).collect(), val: values.to_vec() }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn oscillator_closed_form_across_resonance_and_at_zero() {
    // m = 1, k = 4, c = 0.04; e^(i ωt) gives a negative displacement phase.
    for omega in [0.0_f64, 0.01, 1.98, 2.0, 2.02, 3.0] {
        let (a, b) = (4.0 - omega * omega, omega * 0.04);
        let out = solve(&dense(&[a], 1), &dense(&[b], 1), &Vector { real: vec![1.0], imaginary: vec![0.0] }, &Limits::default()).unwrap();
        let denominator = a * a + b * b;
        close(out.displacement.real[0], a / denominator, 2e-14);
        close(out.displacement.imaginary[0], -b / denominator, 2e-14);
        close(out.displacement.amplitudes()[0], 1.0 / a.hypot(b), 2e-14);
        close(out.displacement.phases()[0], -b.atan2(a), 2e-14);
        assert!(out.residual.relative < 1e-14);
        assert_eq!(out.factor_equations, if omega == 0.0 { 1 } else { 2 });
        if omega == 2.0 {
            close(out.displacement.imaginary[0], -12.5, 2e-14);
        }
    }
}

#[test]
fn structural_loss_is_not_silently_switched_off_at_zero() {
    // Constant loss stiffness G = ηK remains in Z(0). Its limit is NOT K^-1.
    let out = solve(&dense(&[4.0], 1), &dense(&[0.4], 1), &Vector { real: vec![1.0], imaginary: vec![0.0] }, &Limits::default()).unwrap();
    close(out.displacement.real[0], 4.0 / 16.16, 2e-14);
    close(out.displacement.imaginary[0], -0.4 / 16.16, 2e-14);
}

#[test]
fn phase_and_transpose_reciprocity_for_nonproportional_damping() {
    let a = dense(&[2.36, -1.0, -1.0, 1.36], 2);
    let b = dense(&[0.212, -0.062, -0.062, 0.106], 2);
    let rhs = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.0] };
    let first = solve(&a, &b, &rhs, &Limits::default()).unwrap();
    let second = solve(&a, &b, &Vector { real: vec![0.0, 1.0], imaginary: vec![0.0, 0.0] }, &Limits::default()).unwrap();
    close(first.displacement.real[1], second.displacement.real[0], 2e-14);
    close(first.displacement.imaginary[1], second.displacement.imaginary[0], 2e-14);
    assert!(first.displacement.imaginary[1].abs() > 0.01); // reciprocity is not conjugate symmetry.
    let rotated = solve(&a, &b, &Vector { real: vec![0.0, 0.0], imaginary: vec![1.0, 0.0] }, &Limits::default()).unwrap();
    for i in 0..2 {
        close(rotated.displacement.real[i], -first.displacement.imaginary[i], 2e-14);
        close(rotated.displacement.imaginary[i], first.displacement.real[i], 2e-14);
    }
    // Average input power = ω/2 Im(f^T conj(u)); dissipation = ω/2 u^H B u.
    let omega = 0.8;
    let u = &first.displacement;
    let input = -0.5 * omega * u.imaginary[0];
    let mut br = vec![0.0; 2];
    let mut bi = vec![0.0; 2];
    b.mul(&u.real, &mut br);
    b.mul(&u.imaginary, &mut bi);
    let dissipated = 0.5 * omega * (0..2).map(|i| u.real[i] * br[i] + u.imaginary[i] * bi[i]).sum::<f64>();
    assert!(input > 0.0);
    close(input, dissipated, 1e-14);
}

#[test]
fn real_indefinite_system_solves_both_right_hand_sides() {
    let a = dense(&[0.0, 2.0, 2.0, 0.0], 2);
    let b = dense(&[0.0; 4], 2);
    let out = solve(&a, &b, &Vector { real: vec![4.0, 6.0], imaginary: vec![8.0, 10.0] }, &Limits::default()).unwrap();
    assert_eq!(out.displacement, Vector { real: vec![3.0, 2.0], imaginary: vec![5.0, 4.0] });
    assert_eq!(out.factor_equations, 2);
    assert_eq!(out.residual.relative, 0.0);
}

#[test]
fn exact_undamped_resonance_is_a_located_refusal() {
    let a = dense(&[0.0], 1);
    let rhs = Vector { real: vec![1.0], imaginary: vec![0.0] };
    let e = solve(&a, &a, &rhs, &Limits::default()).unwrap_err();
    assert!(matches!(e, Error::Singular(_)));
    assert!(e.to_string().contains("equation 0"));
}

#[test]
fn residual_is_against_original_complex_equations() {
    let a = dense(&[3.0, -1.0, -1.0, 2.0], 2);
    let b = dense(&[0.2, 0.05, 0.05, 0.1], 2);
    let rhs = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.0] };
    let bad = Vector { real: vec![0.6, 0.3], imaginary: vec![0.0, 0.0] };
    assert!(residual(&a, &b, &rhs, &bad, &Limits::default()).unwrap().relative > 0.5);
}

#[test]
fn residual_norm_includes_every_complex_equation_and_checks_zero_force() {
    let a = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let b = dense(&[0.0; 4], 2);
    let zero = Vector::zeros(2);
    let rhs = Vector { real: vec![3.0, 0.0], imaginary: vec![0.0, 4.0] };
    let r = residual(&a, &b, &rhs, &zero, &Limits::default()).unwrap();
    assert_eq!(r.norm, 5.0);
    assert_eq!(r.relative, 1.0);
    assert_eq!(residual(&a, &b, &zero, &zero, &Limits::default()).unwrap().norm, 0.0);
    assert!(residual(&a, &b, &zero, &rhs, &Limits::default()).is_err());
}

#[test]
fn an_unrepresentable_residual_ratio_is_an_error_not_an_infinite_certificate() {
    let a = dense(&[1.0], 1);
    let b = dense(&[0.0], 1);
    let tiny_force = Vector { real: vec![1e-310], imaginary: vec![0.0] };
    let bad = Vector { real: vec![1.0], imaginary: vec![0.0] };
    let error = residual(&a, &b, &tiny_force, &bad, &Limits::default()).unwrap_err();
    assert!(matches!(error, Error::Invalid(_)));
    assert!(error.to_string().contains("residual ratios"));
    let out = solve(&a, &b, &tiny_force, &Limits::default()).unwrap();
    assert_eq!(out.displacement, tiny_force);
    assert_eq!(out.residual.relative, 0.0);
}

#[test]
fn equation_residual_preserves_small_terms_under_cancellation() {
    let a = dense(&[1e16, 1.0, -1e16, 1.0, 2.0, 0.0, -1e16, 0.0, 1e16], 3);
    let b = dense(&[0.0; 9], 3);
    let rhs = Vector { real: vec![1.0, 3.0, 0.0], imaginary: vec![0.0; 3] };
    let u = Vector { real: vec![1.0; 3], imaginary: vec![0.0; 3] };
    assert_eq!(residual(&a, &b, &rhs, &u, &Limits::default()).unwrap().norm, 0.0);
}

#[test]
fn malformed_csr_nonfinite_and_resource_limits_are_refused() {
    let a = dense(&[2.0, -1.0, -1.0, 2.0], 2);
    let b = dense(&[0.1, 0.0, 0.0, 0.1], 2);
    let rhs = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0; 2] };
    let mut malformed = vec![
        Csr { n: usize::MAX, ptr: vec![], idx: vec![], val: vec![] },
        Csr { n: 2, ptr: vec![0, usize::MAX, 4], ..a.clone() },
        Csr { n: 2, ptr: vec![0, 4, 3], ..a.clone() },
    ];
    for (at, value) in [(0, f64::NAN), (1, f64::INFINITY), (1, 0.5)] {
        let mut m = a.clone();
        m.val[at] = value;
        malformed.push(m);
    }
    for cols in [vec![0, 2, 0, 1], vec![0, 0, 0, 1], vec![1, 0, 0, 1]] {
        malformed.push(Csr { idx: cols, ..a.clone() });
    }
    for matrix in malformed {
        assert!(solve(&matrix, &b, &rhs, &Limits::default()).is_err());
    }
    assert!(solve(&a, &b, &Vector { real: vec![], imaginary: vec![] }, &Limits::default()).is_err());
    assert!(solve(&a, &b, &rhs, &Limits { max_factor_entries: 3, ..Limits::default() }).is_err());
    assert!(solve(&a, &b, &rhs, &Limits { max_dofs: 1, ..Limits::default() }).is_err());
}

#[test]
fn interleaving_preserves_sparse_pairs_and_thread_bit_identity() {
    let n = 80;
    let mut a = Csr { n, ptr: vec![0], idx: vec![], val: vec![] };
    for r in 0..n {
        for c in r.saturating_sub(1)..=(r + 1).min(n - 1) {
            a.idx.push(c);
            a.val.push(if c == r { 0.0 } else { -1.0 });
        }
        a.ptr.push(a.idx.len());
    }
    let b = Csr { val: a.idx.iter().map(|_| 0.2).collect(), ..a.clone() };
    let rhs = Vector { real: (0..n).map(|r| (r % 7) as f64).collect(), imaginary: vec![0.0; n] };
    let run = || solve(&a, &b, &rhs, &Limits::default()).unwrap();
    let reference = fieldcraft_parallel::with_threads(1, run);
    for threads in [2, 4, 6] {
        let actual = fieldcraft_parallel::with_threads(threads, run);
        assert_eq!(actual.factor_entries, reference.factor_entries);
        for (x, y) in actual
            .displacement
            .real
            .iter()
            .chain(&actual.displacement.imaginary)
            .zip(reference.displacement.real.iter().chain(&reference.displacement.imaginary))
        {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }
    assert!(reference.factor_entries < n * 12);
    assert!(reference.residual.relative < 1e-13);
}
