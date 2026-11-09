//! Independent continuum fields and energy integrals for Lee–Bathe MITC6-a.
use fieldcraft_fem::{
    basis::Basis,
    dense::{Cholesky, Matrix},
    frame::cross,
    mitc8::Mitc6,
    shell::ShellSection,
};
use fieldcraft_model::Elastic;
fn section(h: f64) -> ShellSection {
    ShellSection::new(h, 0.0, &Elastic::Isotropic { young: 1200.0, poisson: 0.0 }).unwrap()
}
fn nodes(curve: bool) -> Vec<[f64; 3]> {
    Basis::Tri6.nodes().iter().map(|p| [p[0], p[1], if curve { 0.15 * p[0] * p[0] } else { 0.0 }]).collect()
}
fn field(x: &[[f64; 3]], f: impl Fn([f64; 3]) -> [f64; 6]) -> Vec<f64> {
    x.iter().flat_map(|p| f(*p)).collect()
}
fn energy(k: &Matrix, u: &[f64]) -> f64 {
    0.5 * k.mul_vec(u).iter().zip(u).map(|(a, b)| a * b).sum::<f64>()
}
#[test]
fn rigid_and_numbering_invariance_on_curved_geometry() {
    for curved in [false, true] {
        let x = nodes(curved);
        let shell = Mitc6::new(&x, section(0.1)).unwrap();
        let k = shell.stiffness().unwrap();
        for dof in 0..6 {
            let mut v = [0.0; 3];
            v[dof % 3] = 1.0;
            let u = field(&x, |p| {
                if dof < 3 {
                    [v[0], v[1], v[2], 0.0, 0.0, 0.0]
                } else {
                    let t = cross(v, p);
                    [t[0], t[1], t[2], v[0], v[1], v[2]]
                }
            });
            for q in [[0.2, 0.3, -0.8], [0.6, 0.1, 0.7], [1.0 / 3.0, 1.0 / 3.0, 0.0]] {
                for row in shell.strain_operator(q, true).unwrap() {
                    assert!(row.iter().zip(&u).map(|(a, b)| a * b).sum::<f64>().abs() < 2e-12);
                }
            }
            assert!(k.mul_vec(&u).iter().all(|v| v.abs() < 2e-10));
        }
        for perm in [[1, 2, 0, 4, 5, 3], [0, 2, 1, 5, 4, 3]] {
            let xp: Vec<_> = perm.iter().map(|i| x[*i]).collect();
            let kp = Mitc6::new(&xp, section(0.1)).unwrap().stiffness().unwrap();
            let scale = k.data.iter().map(|v| v.abs()).fold(0.0, f64::max);
            for i in 0..36 {
                for j in 0..36 {
                    assert!(
                        (kp.get(i, j) - k.get(6 * perm[i / 6] + i % 6, 6 * perm[j / 6] + j % 6)).abs() < 2e-11 * scale,
                        "numbering {curved} {i} {j}"
                    );
                }
            }
        }
        let mut constrained = Matrix::zeros(30);
        for i in 0..30 {
            for j in 0..30 {
                constrained.add(i, j, k.get(i + 6, j + 6));
            }
        }
        Cholesky::new(&constrained).unwrap();
    }
}
#[test]
fn membrane_bending_shear_and_mass_integrals() {
    let x = nodes(false);
    for h in [0.2, 0.02, 0.002, 0.0002] {
        let shell = Mitc6::new(&x, section(h)).unwrap();
        let k = shell.stiffness().unwrap();
        for (u, want, tol) in [
            (
                field(&x, |p| [0.001 * p[0], 0.002 * p[1], 0.0, 0.0, 0.0, 0.0]),
                0.5 * 0.5 * 1200.0 * h * (0.001_f64.powi(2) + 0.002_f64.powi(2)),
                2e-12,
            ),
            (field(&x, |p| [0.0, 0.0, -0.5 * p[0] * p[0], 0.0, p[0], 0.0]), 0.5 * 0.5 * 1200.0 * h.powi(3) / 12.0, 2e-7),
            (field(&x, |p| [0.0, 0.0, p[0], 0.0, 0.0, 0.0]), 0.5 * 0.5 * (5.0 / 6.0) * 600.0 * h, 2e-12),
        ] {
            let err = (energy(&k, &u) / want - 1.0).abs();
            println!("h={h} energy relative error={err}");
            assert!(err < tol);
        }
        for m in [shell.mass(3.0).unwrap(), shell.lumped_mass(3.0).unwrap()] {
            let u = field(&x, |_| [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
            assert!((2.0 * energy(&m, &u) / (0.5 * 3.0 * h) - 1.0).abs() < 2e-12);
            Cholesky::new(&m).unwrap();
        }
    }
}

#[test]
fn laminated_mitc6a_single_ply_coupling_and_material_directions() {
    use fieldcraft_fem::laminate::PlyProperties;
    use fieldcraft_model::MaterialOrientation;
    let ply = |young, thickness, angle| PlyProperties {
        elastic: Elastic::Isotropic { young, poisson: 0.0 },
        thickness,
        angle,
        points: 3,
        density: None,
        strength: None,
    };
    for curved in [false, true] {
        let xyz = nodes(curved);
        let sec = ShellSection::new(0.1, 0.2, &Elastic::Isotropic { young: 1200.0, poisson: 0.0 }).unwrap();
        let h = Mitc6::new(&xyz, sec).unwrap();
        let l = Mitc6::laminated(&xyz, &[ply(1200.0, 0.1, 0.37)], 0.2, None, 1).unwrap();
        let a = h.stiffness().unwrap();
        let b = l.stiffness().unwrap();
        assert!(a.data.iter().zip(&b.data).all(|(a, b)| (a - b).abs() < 2e-12));
    }
    let xyz = nodes(false);
    let plies = [ply(1200.0, 0.04, 0.0), ply(3000.0, 0.06, 0.0)];
    let shell = Mitc6::laminated(&xyz, &plies, 0.0, None, 1).unwrap();
    let (e, kappa) = (0.002, 0.003);
    let u = field(&xyz, |p| [e * p[0], 0.0, -0.5 * kappa * p[0] * p[0], 0.0, kappa * p[0], 0.0]);
    let want = 0.5
        * 0.5
        * [(1200.0, -0.05f64, -0.01f64), (3000.0, -0.01, 0.05)]
            .iter()
            .map(|&(c, a, b)| c * (e * e * (b - a) + e * kappa * (b * b - a * a) + kappa * kappa * (b.powi(3) - a.powi(3)) / 3.0))
            .sum::<f64>();
    let got = energy(&shell.stiffness().unwrap(), &u);
    assert!((got / want - 1.0).abs() < 2e-12, "coupling energy {got} != {want}");
    let r = shell.stresses(&u).unwrap();
    for n in &r.plies {
        assert!((n[0].top.stress[0] - 1200.0 * (e - 0.01 * kappa)).abs() < 1e-12);
        assert!((n[1].bottom.stress[0] - 3000.0 * (e - 0.01 * kappa)).abs() < 1e-12);
        assert_eq!(n[0].top.z, n[1].bottom.z);
    }
    // Orthotropic fibres fixed in physical space, independent of local r/s order.
    let mut orth = ply(1.0, 0.1, 0.31);
    orth.elastic = Elastic::EngineeringConstants { young: [1000.0, 100.0, 100.0], poisson: [0.0; 3], shear: [50.0, 35.0, 20.0] };
    let axes = MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [0.0, 1.0, 0.0] };
    let base = Mitc6::laminated(&xyz, &[orth], 0.0, Some(axes), 1).unwrap().stiffness().unwrap();
    for perm in [[1, 2, 0, 4, 5, 3], [0, 2, 1, 5, 4, 3]] {
        let reflected = perm == [0, 2, 1, 5, 4, 3];
        let x: Vec<_> = perm.iter().map(|&i| xyz[i]).collect();
        let mut p = orth;
        if reflected {
            p.angle = -p.angle;
        }
        let m = Mitc6::laminated(&x, &[p], 0.0, Some(axes), 1).unwrap().stiffness().unwrap();
        for i in 0..36 {
            for j in 0..36 {
                assert!((m.get(i, j) - base.get(6 * perm[i / 6] + i % 6, 6 * perm[j / 6] + j % 6)).abs() < 2e-11);
            }
        }
    }
    // Reversing the signed normal also reverses physical layer order and offset.
    let perm = [0, 2, 1, 5, 4, 3];
    let x: Vec<_> = perm.iter().map(|&i| xyz[i]).collect();
    let a = Mitc6::laminated(&xyz, &plies, 0.2, None, 1).unwrap().stiffness().unwrap();
    let b = Mitc6::laminated(&x, &[plies[1], plies[0]], -0.2, None, 1).unwrap().stiffness().unwrap();
    for i in 0..36 {
        for j in 0..36 {
            assert!((b.get(i, j) - a.get(6 * perm[i / 6] + i % 6, 6 * perm[j / 6] + j % 6)).abs() < 2e-11);
        }
    }
    println!("MITC6-a laminate coupling energy: computed={got:.16e}, independent={want:.16e}");
}

#[test]
fn invalid_topology_and_natural_points_are_rejected() {
    let x = nodes(false);
    let shell = Mitc6::new(&x, section(0.1)).unwrap();
    assert!(Mitc6::new(&x[..5], section(0.1)).is_err());
    assert!(fieldcraft_fem::mitc8::QuadraticShell::<6, 48>::new(&x, section(0.1)).is_err());
    for p in [[0.7, 0.7, 0.0], [-0.1, 0.0, 0.0], [0.2, 0.3, f64::NAN]] {
        assert!(shell.strain_operator(p, true).is_err());
    }
    assert!(shell.stresses(&[0.0; 48]).is_err());
    // Exercise layered coordinate evaluation at all six nodes, guarding against
    // accidentally retaining an eight-node loop in the shared implementation.
    let ply = fieldcraft_fem::laminate::PlyProperties {
        elastic: Elastic::Isotropic { young: 1200.0, poisson: 0.0 },
        thickness: 0.1,
        angle: 0.0,
        points: 3,
        density: Some(3.0),
        strength: None,
    };
    let layer = Mitc6::laminated(&x, &[ply], 0.0, None, 1).unwrap();
    assert_eq!(layer.stresses(&[0.0; 36]).unwrap().plies.len(), 6);
}
