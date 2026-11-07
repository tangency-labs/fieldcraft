use fieldcraft_solver::{Limits, complex::Vector, sensitivity, sparse::Csr};

fn matrix(v: &[f64], n: usize) -> Csr {
    Csr { n, ptr: (0..=n).map(|i| i * n).collect(), idx: (0..n).flat_map(|_| 0..n).collect(), val: v.to_vec() }
}

#[test]
fn oscillator_derivatives_include_resonance_and_static_limit() {
    let limits = Limits::default();
    let f = Vector { real: vec![1.0], imaginary: vec![0.0] };
    let dk = matrix(&[1.0], 1);
    let zero = matrix(&[0.0], 1);
    for omega in [0.0_f64, 1.9, 2.0, 2.1] {
        let a = matrix(&[4.0 - omega * omega], 1);
        let b = matrix(&[0.12 * omega], 1);
        let u = sensitivity::direct(&a, &b, &f, &limits).unwrap().displacement;
        let load = sensitivity::tangent_load(&dk, &zero, &u, &Vector::zeros(1), &limits).unwrap();
        let du = sensitivity::direct(&a, &b, &load, &limits).unwrap();
        let expected = [u.imaginary[0] * u.imaginary[0] - u.real[0] * u.real[0], -2.0 * u.real[0] * u.imaginary[0]];
        assert!((du.displacement.real[0] - expected[0]).abs() < 1e-12);
        assert!((du.displacement.imaginary[0] - expected[1]).abs() < 1e-12);
        let ap = sensitivity::amplitude_phase([u.real[0], u.imaginary[0]], expected).unwrap();
        let norm = a.val[0].hypot(b.val[0]);
        assert!((ap[0] + a.val[0] / norm.powi(3)).abs() < 1e-11);
        assert!((ap[1] - b.val[0] / norm.powi(2)).abs() < 1e-11);
        for h in [1e-3, 1e-4, 1e-5] {
            let plus = sensitivity::direct(&matrix(&[a.val[0] + h], 1), &b, &f, &limits).unwrap().displacement;
            let minus = sensitivity::direct(&matrix(&[a.val[0] - h], 1), &b, &f, &limits).unwrap().displacement;
            for (fd, exact) in
                [((plus.real[0] - minus.real[0]) / (2.0 * h), expected[0]), ((plus.imaginary[0] - minus.imaginary[0]) / (2.0 * h), expected[1])]
            {
                assert!((fd - exact).abs() < 3e-5 * exact.abs().max(1.0));
            }
        }
    }
}

#[test]
fn complex_transpose_adjoint_matches_direct_observation() {
    let a = matrix(&[3.0, -1.0, -1.0, 2.0], 2);
    let b = matrix(&[0.2, 0.03, 0.03, 0.4], 2);
    let rhs = Vector { real: vec![0.4, -0.7], imaginary: vec![0.8, 0.3] };
    let weights = Vector { real: vec![0.3, 0.9], imaginary: vec![-0.4, 0.2] };
    let limits = Limits::default();
    let direct = sensitivity::direct(&a, &b, &rhs, &limits).unwrap().displacement;
    let adjoint = sensitivity::adjoint(&a, &b, &weights, &rhs, &limits).unwrap();
    let real: f64 = (0..2).map(|i| weights.real[i] * direct.real[i] - weights.imaginary[i] * direct.imaginary[i]).sum();
    let imag: f64 = (0..2).map(|i| weights.real[i] * direct.imaginary[i] + weights.imaginary[i] * direct.real[i]).sum();
    assert!((real - adjoint[0]).abs() < 1e-14);
    assert!((imag - adjoint[1]).abs() < 1e-14);
}

#[test]
fn invalid_inputs_are_errors_and_zero_phase_is_not_differentiable() {
    assert!(sensitivity::amplitude_phase([0.0, 0.0], [1.0, 0.0]).is_err());
    assert!(sensitivity::amplitude_phase([f64::NAN, 0.0], [1.0, 0.0]).is_err());
    let k = matrix(&[1.0], 1);
    assert!(sensitivity::tangent_load(&k, &k, &Vector::zeros(0), &Vector::zeros(1), &Limits::default()).is_err());
    let bad = Csr { ptr: vec![0, 2], ..k.clone() };
    assert!(sensitivity::tangent_load(&bad, &k, &Vector::zeros(1), &Vector::zeros(1), &Limits::default()).is_err());
}

#[test]
fn tangent_load_differentiates_both_operator_parts_and_force() {
    let da = matrix(&[2.0, -1.0, -1.0, 3.0], 2);
    let db = matrix(&[0.5, 0.2, 0.2, -0.1], 2);
    let u = Vector { real: vec![1.0, 2.0], imaginary: vec![-2.0, 3.0] };
    let df = Vector { real: vec![0.7, -0.3], imaginary: vec![0.1, 0.6] };
    let rhs = sensitivity::tangent_load(&da, &db, &u, &df, &Limits::default()).unwrap();
    for (a, b) in rhs.real.iter().chain(&rhs.imaginary).zip([0.3, -6.0, 6.2, -10.4]) {
        assert!((a - b).abs() < 2e-15, "{a} != {b}");
    }
}
