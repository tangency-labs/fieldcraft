//! Independent closed forms and published stress-state checks. These are kernel tests;
//! they do not stand in for the command oracle or a finite-element benchmark solve.
use fieldcraft_fem::failure::plane_stress_indices;
use fieldcraft_fem::laminate::{Lamina, Laminate, PlyProperties};
use fieldcraft_model::{Elastic, PlyStrength};

const ORTHO: Elastic = Elastic::EngineeringConstants { young: [100.0, 10.0, 10.0], poisson: [0.0; 3], shear: [5.0, 4.0, 3.0] };

fn ply(angle: f64, thickness: f64) -> PlyProperties {
    PlyProperties { elastic: ORTHO, thickness, angle, points: 3, density: Some(2.0), strength: None }
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() < 3e-14 * scale, "{a:.16e} != {b:.16e}");
}

#[test]
fn three_stacking_sequences_have_the_closed_form_abd() {
    // Q0=diag(100,10,5); no Poisson coupling. h=1. Outer two quarter plies
    // have ∫z²dz=7/96, inner two 1/96. These constants are not computed by the kernel.
    let half_pi = std::f64::consts::FRAC_PI_2;
    let qtr_pi = std::f64::consts::FRAC_PI_4;
    let cases = [
        (
            vec![ply(0.0, 0.25), ply(half_pi, 0.25), ply(half_pi, 0.25), ply(0.0, 0.25)],
            [55.0, 55.0, 0.0, 5.0],
            [710.0 / 96.0, 170.0 / 96.0, 0.0, 5.0 / 12.0, 0.0],
            [0.0, 0.0],
        ),
        (
            vec![ply(qtr_pi, 0.25), ply(-qtr_pi, 0.25), ply(-qtr_pi, 0.25), ply(qtr_pi, 0.25)],
            [32.5, 32.5, 22.5, 27.5],
            [32.5 / 12.0, 32.5 / 12.0, 22.5 / 12.0, 27.5 / 12.0, 22.5 / 16.0],
            [0.0, 0.0],
        ),
        (vec![ply(0.0, 0.5), ply(half_pi, 0.5)], [55.0, 55.0, 0.0, 5.0], [55.0 / 12.0, 55.0 / 12.0, 0.0, 5.0 / 12.0, 0.0], [-11.25, 11.25]),
    ];
    for (plies, a, d, b) in cases {
        let section = Laminate::new(&plies, 0.0).unwrap();
        let k = section.abd();
        let mut expected = [[0.0; 6]; 6];
        expected[0][0] = a[0];
        expected[1][1] = a[1];
        expected[0][1] = a[2];
        expected[1][0] = a[2];
        expected[2][2] = a[3];
        expected[3][3] = d[0];
        expected[4][4] = d[1];
        expected[3][4] = d[2];
        expected[4][3] = d[2];
        expected[5][5] = d[3];
        expected[3][5] = d[4];
        expected[5][3] = d[4];
        expected[4][5] = d[4];
        expected[5][4] = d[4];
        expected[0][3] = b[0];
        expected[3][0] = b[0];
        expected[1][4] = b[1];
        expected[4][1] = b[1];
        for i in 0..6 {
            for j in 0..6 {
                close(k[i][j], expected[i][j], 100.0);
                assert_eq!(k[i][j], k[j][i]);
            }
        }
        println!(
            "angles={:?}, A11={:.15}, B11={:.15}, D11={:.15}, D16={:.15}",
            plies.iter().map(|p| p.angle).collect::<Vec<_>>(),
            k[0][0],
            k[0][3],
            k[3][3],
            k[3][5]
        );
    }
}

#[test]
fn eccentricity_shifts_stiffness_and_density_moments_together() {
    let plies = [ply(0.2, 0.25), PlyProperties { density: Some(7.0), ..ply(-0.7, 0.75) }];
    let center = Laminate::new(&plies, 0.0).unwrap();
    let offset = Laminate::new(&plies, 0.3).unwrap();
    let (c, o, e) = (center.abd(), offset.abd(), -0.3);
    for i in 0..3 {
        for j in 0..3 {
            close(o[i][j], c[i][j], 100.0);
            close(o[i][j + 3], c[i][j + 3] + e * c[i][j], 100.0);
            close(o[i + 3][j + 3], c[i + 3][j + 3] + 2.0 * e * c[i][j + 3] + e * e * c[i][j], 100.0);
        }
    }
    let m = center.mass_moments().unwrap();
    let n = offset.mass_moments().unwrap();
    close(m[0], 5.75, 10.0);
    close(m[1], 0.46875, 1.0); // asymmetric density despite a centred reference plane
    close(n[0], m[0], 10.0);
    close(n[1], m[1] + e * m[0], 10.0);
    close(n[2], m[2] + 2.0 * e * m[1] + e * e * m[0], 10.0);
    let single = Laminate::new(&[ply(0.0, 1.0)], 0.5).unwrap();
    assert_eq!(single.geometric_moments(), [1.0, -0.5, 1.0 / 3.0]);
}

#[test]
fn off_axis_stress_and_strain_are_energy_conjugate() {
    let p = ply(std::f64::consts::FRAC_PI_4, 1.0);
    let section = Laminate::new(&[p], 0.0).unwrap();
    let g = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let fields = section.recover(g).unwrap();
    let at = &fields[0].middle;
    close(at.strain[0], 0.5, 1.0);
    close(at.strain[1], 0.5, 1.0);
    close(at.strain[3], -1.0, 1.0);
    close(at.stress[0], 50.0, 100.0);
    close(at.stress[1], 5.0, 10.0);
    close(at.stress[3], -5.0, 10.0);
    let local_work: f64 = at.strain.iter().zip(at.stress).map(|(e, s)| e * s).sum();
    let generalized_work: f64 = g.iter().zip(section.resultants(g).unwrap()).map(|(e, s)| e * s).sum();
    close(local_work, generalized_work, 100.0);
    let q = Lamina::new(&ORTHO).unwrap().rotated(p.angle).unwrap();
    close(q[0][0], 32.5, 100.0);
    close(q[0][1], 22.5, 100.0);
    close(q[0][2], 22.5, 100.0);
    let ks = section.shear();
    close(ks[0][1], 5.0 / 12.0, 1.0);
}

#[test]
fn one_sided_interface_fields_and_single_point_integration_preserve_bending() {
    let p0 = ply(0.0, 0.5);
    let p1 = ply(std::f64::consts::FRAC_PI_2, 0.5);
    let section = Laminate::new(&[p0, p1], 0.0).unwrap();
    let g = [0.1, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0];
    let f = section.recover(g).unwrap();
    assert_eq!(f[0].bottom.z, -0.5);
    assert_eq!(f[0].top.z, f[1].bottom.z);
    assert_eq!(f[1].top.z, 0.5);
    close(f[0].bottom.stress[0], 0.0, 1.0);
    close(f[0].top.stress[0], 10.0, 100.0);
    close(f[1].bottom.stress[1], 1.0, 10.0);
    close(f[1].top.stress[1], 2.0, 10.0);
    let one = Laminate::new(&[PlyProperties { points: 1, ..p0 }, PlyProperties { points: 1, ..p1 }], 0.0).unwrap();
    assert_eq!(one.abd(), section.abd());
    assert_eq!(one.mass_moments().unwrap(), section.mass_moments().unwrap());
    assert_eq!(one.recover(g).unwrap()[0].points.len(), 1);
    let sampled = section.recover(g).unwrap();
    close(sampled[0].points[0].z, -0.25 - 0.25 * (3.0_f64 / 5.0).sqrt(), 1.0);
    for points in 1..=15 {
        let layer = Laminate::new(&[PlyProperties { points, ..p0 }], 0.0).unwrap().recover(g).unwrap();
        assert_eq!(layer[0].points.len(), points as usize);
        assert!(layer[0].points.windows(2).all(|w| w[0].z < w[1].z));
        for p in &layer[0].points {
            assert!(p.z > layer[0].bottom.z && p.z < layer[0].top.z);
        }
    }
    // Integrate linearly varying in-plane stress with Simpson's rule independently of ABD.
    let flat = Laminate::new(&[ply(0.0, 1.0)], 0.0).unwrap();
    let f = flat.recover(g).unwrap();
    let p = &f[0];
    let force = (p.bottom.stress[0] + 4.0 * p.middle.stress[0] + p.top.stress[0]) / 6.0;
    let moment = (p.bottom.z * p.bottom.stress[0] + 4.0 * p.middle.z * p.middle.stress[0] + p.top.z * p.top.stress[0]) / 6.0;
    let resultant = flat.resultants(g).unwrap();
    close(force, resultant[0], 100.0);
    close(moment, resultant[3], 100.0);
}

fn strengths() -> PlyStrength {
    PlyStrength { tension: [500.0, 5.0], compression: [350.0, 75.0], shear12: 35.0, shear23: 35.0, tsai_wu_interaction: -0.5 }
}

#[test]
fn failure_thresholds_and_sign_branches_are_not_reserve_factors() {
    let s = strengths();
    for stress in [[500.0, 0.0, 0.0], [-350.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, -75.0, 0.0], [0.0, 0.0, 35.0]] {
        let f = plane_stress_indices(stress, &s).unwrap();
        close(f.maximum_stress, 1.0, 1.0);
        close(f.tsai_hill, 1.0, 1.0);
        close(f.tsai_wu, 1.0, 1.0);
        close(f.hashin.maximum, 1.0, 1.0);
    }
    let f = plane_stress_indices([-175.0, 0.0, 0.0], &s).unwrap();
    assert_eq!(f.maximum_stress, 0.5);
    assert_eq!(f.tsai_hill, 0.25);
    assert_eq!(f.hashin.fiber_compression, Some(0.25));
    assert_eq!(f.hashin.fiber_tension, None);
    let different_s23 = PlyStrength { shear23: 15.0, ..s };
    let a = plane_stress_indices([0.0, -10.0, 3.0], &s).unwrap();
    let b = plane_stress_indices([0.0, -10.0, 3.0], &different_s23).unwrap();
    assert_ne!(a.hashin.matrix_compression, b.hashin.matrix_compression);
    assert_eq!(a.tsai_wu, b.tsai_wu);
    let extreme = PlyStrength { tension: [1e18, 5.0], compression: [1.0, 1e9], shear23: 1.0, ..s };
    close(plane_stress_indices([1e18, 0.0, 0.0], &extreme).unwrap().tsai_wu, 1.0, 1.0);
    close(plane_stress_indices([0.0, -1e9, 0.0], &extreme).unwrap().hashin.matrix_compression.unwrap(), 1.0, 1.0);
}

#[test]
fn nafems_r0093_published_stress_states_keep_the_negative_tsai_wu_value() {
    // Hopkins, NAFEMS R0093 (2005), problem 1; public reproduction: COMSOL 6.3
    // "Failure Prediction in a Laminated Composite Shell", Tables 8–12:
    // https://doc.comsol.com/6.3/doc/com.comsol.help.models.compmat.failure_prediction_in_a_laminated_composite_shell/failure_prediction_in_a_laminated_composite_shell.html
    // Stresses below are benchmark (not computed-column) values in MPa. The index is
    // dimensionless, so strengths use the same units. 2e-4 absolute allows rounding in
    // the published stresses and FI (not a finite-element discretization tolerance).
    let stresses = [[-5.128, 4.407, -1.663], [12.59, 1.983, 2.572], [8.520, 0.1256, -2.051], [9.357, -1.859, -0.5557]];
    for ((stress, tw), th) in stresses.into_iter().zip([0.8840, 0.3730, 0.0199, -0.3430]).zip([0.7795, 0.1632, 0.0043, 0.00140]) {
        let f = plane_stress_indices(stress, &strengths()).unwrap();
        assert!((f.tsai_wu - tw).abs() < 2e-4, "TW: {} != {tw}", f.tsai_wu);
        assert!((f.tsai_hill - th).abs() < 2e-4, "TH: {} != {th}", f.tsai_hill);
        println!("NAFEMS stress={stress:?}, TW={:.9}, TH={:.9}", f.tsai_wu, f.tsai_hill);
    }
}

#[test]
fn hostile_laminates_and_missing_density_are_refused() {
    let good = ply(0.0, 1.0);
    for bad in [
        PlyProperties { thickness: 0.0, ..good },
        PlyProperties { thickness: f64::MAX, ..good },
        PlyProperties { thickness: f64::MIN_POSITIVE, ..good },
        PlyProperties { angle: f64::NAN, ..good },
        PlyProperties { points: 16, ..good },
        PlyProperties { density: Some(f64::NAN), ..good },
    ] {
        assert!(Laminate::new(&[bad], 0.0).is_err(), "{bad:?}");
    }
    assert!(Laminate::new(&[], 0.0).is_err());
    assert!(Laminate::new(&vec![good; 257], 0.0).is_err());
    for bad in [f64::NAN, f64::INFINITY, 0.51] {
        assert!(Laminate::new(&[good], bad).is_err());
    }
    let missing = Laminate::new(&[PlyProperties { density: None, ..good }], 0.0).unwrap();
    assert!(missing.mass_moments().is_err());
    assert!(missing.recover([f64::MAX; 8]).is_err());
    assert!(missing.resultants([f64::NAN; 8]).is_err());
    assert!(plane_stress_indices([f64::NAN; 3], &strengths()).is_err());
    assert!(plane_stress_indices([f64::MAX; 3], &strengths()).is_err());
    let tiny = f64::from_bits(1);
    assert!(Laminate::new(&[PlyProperties { density: Some(tiny), thickness: 0.1, ..good }], 0.0).is_err());
    let small_shear = Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0, tiny, 1.0] };
    assert!(Laminate::new(&[PlyProperties { elastic: small_shear, thickness: 0.1, ..good }], 0.0).is_err());
}
