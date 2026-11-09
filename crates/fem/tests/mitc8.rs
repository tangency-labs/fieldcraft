//! Independent polynomial fields and energy integrals; no reference solver.
use fieldcraft_fem::{
    basis::Basis,
    dense::{Cholesky, Matrix},
    frame::cross,
    mitc8::Mitc8,
    shell::{SHEAR_CORRECTION, ShellSection},
};
use fieldcraft_model::Elastic;

fn section(h: f64, offset: f64) -> ShellSection {
    ShellSection::new(h, offset, &Elastic::Isotropic { young: 1200.0, poisson: 0.25 }).unwrap()
}
fn nodes(curved: bool) -> Vec<[f64; 3]> {
    Basis::Quad8.nodes().iter().map(|p| [p[0], p[1], if curved { 0.15 * p[0] * p[0] } else { 0.0 }]).collect()
}
fn field(x: &[[f64; 3]], f: impl Fn([f64; 3]) -> [f64; 6]) -> Vec<f64> {
    x.iter().flat_map(|p| f(*p)).collect()
}
fn energy(k: &Matrix, u: &[f64]) -> f64 {
    0.5 * k.mul_vec(u).iter().zip(u).map(|(a, b)| a * b).sum::<f64>()
}

#[test]
fn six_rigid_motions_vanish_on_curved_and_offset_surfaces() {
    let xyz = nodes(true);
    for offset in [0.0, 0.3] {
        let shell = Mitc8::new(&xyz, section(0.2, offset)).unwrap();
        for dof in 0..6 {
            let mut v = [0.0; 3];
            v[dof % 3] = 1.0;
            let u = field(&xyz, |p| {
                if dof < 3 {
                    [v[0], v[1], v[2], 0.0, 0.0, 0.0]
                } else {
                    let d = cross(v, p);
                    [d[0], d[1], d[2], v[0], v[1], v[2]]
                }
            });
            for p in [[0.0, 0.0, 0.0], [-0.6, 0.3, -1.0], [0.7, -0.5, 1.0]] {
                for assumed in [false, true] {
                    let b = shell.strain_operator(p, assumed).unwrap();
                    for row in b {
                        let e: f64 = row.iter().zip(&u).map(|(a, b)| a * b).sum();
                        assert!(e.abs() < 2e-12, "rigid {dof}: {p:?} {assumed} {e}");
                    }
                }
            }
            let k = shell.stiffness().unwrap();
            let residual = k.mul_vec(&u).iter().map(|v| v.abs()).fold(0.0, f64::max);
            assert!(residual < 2e-10, "rigid stiffness residual {residual}");
        }
    }
}

#[test]
fn constant_membrane_strain_and_curved_axial_extension() {
    for curved in [false, true] {
        let xyz = nodes(curved);
        let shell = Mitc8::new(&xyz, section(0.1, 0.0)).unwrap();
        let u = field(&xyz, |p| [0.0, 0.003 * p[1], 0.0, 0.0, 0.0, 0.0]);
        for p in [[-0.4, 0.7, -0.8], [0.3, -0.2, 0.6], [0.0, 0.0, 0.0]] {
            for assumed in [false, true] {
                let b = shell.strain_operator(p, assumed).unwrap();
                let got: [f64; 5] = std::array::from_fn(|i| b[i].iter().zip(&u).map(|(a, b)| a * b).sum());
                for (a, b) in got.into_iter().zip([0.0, 0.003, 0.0, 0.0, 0.0]) {
                    assert!((a - b).abs() < 2e-13, "{curved}: {got:?}");
                }
            }
        }
    }
    // Distorted straight-sided quadrilateral: affine physical membrane field,
    // checked by total energy (independent polygon area), including in-plane spin.
    let corners = [[0.0, 0.0, 0.0], [2.0, 0.2, 0.0], [1.6, 1.4, 0.0], [-0.3, 1.0, 0.0]];
    let mut xyz = corners.to_vec();
    for i in 0..4 {
        xyz.push(std::array::from_fn(|k| 0.5 * (corners[i][k] + corners[(i + 1) % 4][k])));
    }
    let shell = Mitc8::new(&xyz, section(0.1, 0.0)).unwrap();
    let u = field(&xyz, |p| [0.001 * p[0] + 0.002 * p[1], -0.003 * p[0] + 0.004 * p[1], 0.0, 0.0, 0.0, -0.0025]);
    let area: f64 = (0..4).map(|i| 0.5 * (corners[i][0] * corners[(i + 1) % 4][1] - corners[(i + 1) % 4][0] * corners[i][1])).sum();
    let want = 0.5
        * 0.1
        * area
        * (1200.0 / (1.0 - 0.25f64.powi(2)) * (0.001f64.powi(2) + 0.004f64.powi(2) + 2.0 * 0.25 * 0.001 * 0.004) + 480.0 * 0.001f64.powi(2));
    let got = energy(&shell.stiffness().unwrap(), &u);
    assert!((got / want - 1.0).abs() < 2e-12, "membrane {got} {want}");
}

#[test]
fn bending_and_shear_energy_remain_correct_in_the_thin_limit() {
    let xyz = nodes(false);
    for h in [0.2, 0.02, 0.002, 0.0002] {
        let shell = Mitc8::new(&xyz, section(h, 0.0)).unwrap();
        let k = shell.stiffness().unwrap();
        let u = field(&xyz, |p| [0.0, 0.0, -0.5 * p[0] * p[0], 0.0, p[0], 0.0]);
        let want = 0.5 * 4.0 * 1200.0 / (1.0 - 0.25f64.powi(2)) * h.powi(3) / 12.0;
        let got = energy(&k, &u);
        let bending_error = (got / want - 1.0).abs();
        assert!((got / want - 1.0).abs() < 2e-7, "bending h={h} {got} {want}");
        for row in shell.strain_operator([0.37, -0.25, 0.6], true).unwrap().iter().skip(3) {
            assert!(row.iter().zip(&u).map(|(a, b)| a * b).sum::<f64>().abs() < 1e-12);
        }
        let u = field(&xyz, |p| [0.0, 0.0, p[0], 0.0, 0.0, 0.0]);
        let want = 0.5 * 4.0 * SHEAR_CORRECTION * 480.0 * h;
        assert!((energy(&k, &u) / want - 1.0).abs() < 2e-12);
        println!("MITC8 pure bending h={h}: relative energy error={:.6e}", bending_error);
    }
}

#[test]
fn no_extra_zero_energy_modes_after_six_rigid_constraints() {
    for curved in [false, true] {
        let shell = Mitc8::new(&nodes(curved), section(0.2, 0.0)).unwrap();
        let k = shell.stiffness().unwrap();
        let mut fixed = Matrix::zeros(42);
        for i in 0..42 {
            for j in 0..42 {
                assert_eq!(k.get(i + 6, j + 6), k.get(j + 6, i + 6));
                fixed.add(i, j, k.get(i + 6, j + 6));
            }
        }
        Cholesky::new(&fixed).unwrap();
    }
}

#[test]
fn body_force_and_invalid_geometry() {
    let xyz = nodes(false);
    let shell = Mitc8::new(&xyz, section(0.2, 0.3)).unwrap();
    let f = shell.body_force([1.0, 0.0, 0.0]).unwrap();
    assert!((f.iter().step_by(6).sum::<f64>() - 0.8).abs() < 1e-13);
    assert!((f.iter().skip(4).step_by(6).sum::<f64>() + 0.048).abs() < 1e-13);
    assert!(Mitc8::new(&xyz[..7], section(0.2, 0.0)).is_err());
    let mut folded = xyz.clone();
    folded[4] = [0.0, 3.0, 0.0];
    assert!(Mitc8::new(&folded, section(0.2, 0.0)).is_err());
    assert!(shell.strain_operator([f64::NAN, 0.0, 0.0], true).is_err());
    assert!(shell.stresses(&[0.0; 47]).is_err());
}

#[test]
fn laminated_mitc8_single_ply_coupling_and_material_directions() {
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
        let h = Mitc8::new(&xyz, sec).unwrap();
        let l = Mitc8::laminated(&xyz, &[ply(1200.0, 0.1, 0.37)], 0.2, None, 1).unwrap();
        let a = h.stiffness().unwrap();
        let b = l.stiffness().unwrap();
        assert!(a.data.iter().zip(&b.data).all(|(a, b)| (a - b).abs() < 2e-12));
    }
    let xyz = nodes(false);
    let plies = [ply(1200.0, 0.04, 0.0), ply(3000.0, 0.06, 0.0)];
    let shell = Mitc8::laminated(&xyz, &plies, 0.0, None, 1).unwrap();
    let (e, kappa) = (0.002, 0.003);
    let u = field(&xyz, |p| [e * p[0], 0.0, -0.5 * kappa * p[0] * p[0], 0.0, kappa * p[0], 0.0]);
    let want = 0.5
        * 4.0
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
    let base = Mitc8::laminated(&xyz, &[orth], 0.0, Some(axes), 1).unwrap().stiffness().unwrap();
    for perm in [[1, 2, 3, 0, 5, 6, 7, 4], [0, 3, 2, 1, 7, 6, 5, 4]] {
        let reflected = perm[1] == 3;
        let x: Vec<_> = perm.iter().map(|&i| xyz[i]).collect();
        let mut p = orth;
        if reflected {
            p.angle = -p.angle;
        }
        let m = Mitc8::laminated(&x, &[p], 0.0, Some(axes), 1).unwrap().stiffness().unwrap();
        for i in 0..48 {
            for j in 0..48 {
                assert!((m.get(i, j) - base.get(6 * perm[i / 6] + i % 6, 6 * perm[j / 6] + j % 6)).abs() < 2e-11);
            }
        }
    }
    // Reversing the signed normal also reverses physical layer order and offset.
    let perm = [0, 3, 2, 1, 7, 6, 5, 4];
    let x: Vec<_> = perm.iter().map(|&i| xyz[i]).collect();
    let a = Mitc8::laminated(&xyz, &plies, 0.2, None, 1).unwrap().stiffness().unwrap();
    let b = Mitc8::laminated(&x, &[plies[1], plies[0]], -0.2, None, 1).unwrap().stiffness().unwrap();
    for i in 0..48 {
        for j in 0..48 {
            assert!((b.get(i, j) - a.get(6 * perm[i / 6] + i % 6, 6 * perm[j / 6] + j % 6)).abs() < 2e-11);
        }
    }
    println!("MITC8 laminate coupling energy: computed={got:.16e}, independent={want:.16e}");
}

#[test]
fn layered_weight_keeps_the_first_mass_moment() {
    use fieldcraft_fem::laminate::PlyProperties;
    let ply = |density| PlyProperties {
        elastic: Elastic::Isotropic { young: 1200.0, poisson: 0.0 },
        thickness: 0.05,
        angle: 0.0,
        points: 3,
        density: Some(density),
        strength: None,
    };
    let shell = Mitc8::laminated(&nodes(false), &[ply(2.0), ply(4.0)], 0.0, None, 1).unwrap();
    let f = shell.gravity(3.0, [1.0, 0.0, 0.0]).unwrap();
    let force: f64 = f.as_chunks::<6>().0.iter().map(|v| v[0]).sum();
    let moment: f64 = f.as_chunks::<6>().0.iter().map(|v| v[4]).sum();
    assert!((force - 1.2).abs() < 1e-13);
    assert!((moment - 0.01).abs() < 1e-13);
}
