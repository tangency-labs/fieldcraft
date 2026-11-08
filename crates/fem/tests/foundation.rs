use fieldcraft_fem::foundation::matrix;

#[test]
fn triangle_exact_consistent_matrix_and_local_directions() {
    // Right triangle, area 1, in the yz plane: x is normal.
    let xyz = [[0., 0., 0.], [0., 2., 0.], [0., 0., 1.]];
    let k = matrix(&xyz, 12., 6.).unwrap();
    for a in 0..3 {
        for b in 0..3 {
            for r in 0..3 {
                for c in 0..3 {
                    let expected = if r != c { 0. } else { (if a == b { 2. } else { 1. }) * (if r == 0 { 1. } else { 0.5 }) };
                    assert!((k[(3 * a + r) * 9 + 3 * b + c] - expected).abs() < 2e-14);
                }
            }
        }
    }
}

#[test]
fn all_face_orders_preserve_rigid_translation_energy() {
    let tri = vec![[0., 0., 0.], [2., 0., 0.], [0., 1., 0.]];
    let quad = vec![[0., 0., 0.], [2., 0., 0.], [2., 1., 0.], [0., 1., 0.]];
    for corners in [tri, quad] {
        for quadratic in [false, true] {
            let mut xyz = corners.clone();
            if quadratic {
                for i in 0..corners.len() {
                    let a = corners[i];
                    let b = corners[(i + 1) % corners.len()];
                    xyz.push([(a[0] + b[0]) / 2., (a[1] + b[1]) / 2., 0.]);
                }
            }
            let k = matrix(&xyz, 7., 3.).unwrap();
            let n = xyz.len();
            let area = if corners.len() == 3 { 1. } else { 2. };
            for c in 0..3 {
                let sum: f64 = (0..n).flat_map(|a| (0..n).map(move |b| (a, b))).map(|(a, b)| k[(3 * a + c) * 3 * n + 3 * b + c]).sum();
                assert!((sum - area * if c == 2 { 7. } else { 3. }).abs() < 1e-12);
            }
            for a in 0..3 * n {
                for b in 0..3 * n {
                    assert!((k[a * 3 * n + b] - k[b * 3 * n + a]).abs() < 1e-14);
                }
            }
        }
    }
}

#[test]
fn rejects_bad_geometry_and_coefficients_without_panicking() {
    assert!(matrix(&[], 1., 0.).is_err());
    assert!(matrix(&[[0.; 3]; 3], 1., 0.).is_err());
    let xyz = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    assert!(matrix(&xyz, f64::NAN, 0.).is_err());
    assert!(matrix(&xyz, 1., -1.).is_err());
    assert!(matrix(&xyz, 0., 0.).unwrap().iter().all(|v| *v == 0.));
}

#[test]
fn quadratic_triangle_integrates_shape_products_not_just_area() {
    // ∫L1²(2L1-1)² dA = A/30; ∫(4L1L2)² dA = 8A/45.
    // These moments detect an under-integrated quadratic face even if its area is right.
    let xyz = [[0., 0., 0.], [2., 0., 0.], [0., 1., 0.], [1., 0., 0.], [1., 0.5, 0.], [0., 0.5, 0.]];
    let k = matrix(&xyz, 1., 0.).unwrap();
    for a in 0..6 {
        let expected = if a < 3 { 1. / 30. } else { 8. / 45. };
        assert!((k[(3 * a + 2) * 18 + 3 * a + 2] - expected).abs() < 1e-14);
    }
}

#[test]
fn coefficients_are_linear_at_zero_and_normals_are_unoriented() {
    let xyz = [[0., 0., 0.], [1., 1., 0.], [0., 0., 1.]];
    let kn = matrix(&xyz, 1., 0.).unwrap();
    let kt = matrix(&xyz, 0., 1.).unwrap();
    let combined = matrix(&xyz, 2., 3.).unwrap();
    for ((n, t), k) in kn.iter().zip(&kt).zip(&combined) {
        assert!((2. * n + 3. * t - k).abs() < 1e-14);
    }
    let reversed = matrix(&[xyz[0], xyz[2], xyz[1]], 2., 3.).unwrap();
    let permutation = [0, 2, 1];
    for a in 0..3 {
        for b in 0..3 {
            for r in 0..3 {
                for c in 0..3 {
                    assert!((combined[(3 * a + r) * 9 + 3 * b + c] - reversed[(3 * permutation[a] + r) * 9 + 3 * permutation[b] + c]).abs() < 1e-14);
                }
            }
        }
    }
}
