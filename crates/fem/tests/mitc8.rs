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

#[test]
fn inertia_integrates_rigid_motion_and_density_moments() {
    use fieldcraft_fem::{frame::dot, quadrature::Rule};
    for curved in [false, true] {
        let xyz = nodes(curved);
        let s = Mitc8::new(&xyz, section(0.2, 0.3)).unwrap();
        let m = s.mass(3.0).unwrap();
        Cholesky::new(&m).unwrap();
        let omega = [0.4, -0.3, 0.7];
        let velocity = [0.2, 0.5, -0.1];
        let u = field(&xyz, |x| {
            let v = cross(omega, x);
            [v[0] + velocity[0], v[1] + velocity[1], v[2] + velocity[2], omega[0], omega[1], omega[2]]
        });
        // Independently form x(r,s,t), its three tangents and omega×x+v.
        // Same declared integration rule: this checks the discrete kinetic energy,
        // not an exact quadrature claim for arbitrary curved Jacobians.
        let mut reference = 0.0;
        let mut total = 0.0;
        let mut rotary = 0.0;
        for q in Rule::Quad(3).points() {
            for (t, w) in fieldcraft_fem::quadrature::gauss_legendre(3) {
                let sh = Basis::Quad8.eval(q.xi);
                let z = -0.06 + 0.1 * t;
                let mut x = [0.0; 3];
                let mut g = [[0.0; 3]; 3];
                let mut drill = 0.0;
                for (i, node) in xyz.iter().enumerate() {
                    let n = s.directors()[i];
                    drill += sh.n[i] * dot(omega, n);
                    for a in 0..3 {
                        x[a] += sh.n[i] * (node[a] + z * n[a]);
                        for (b, tangent) in g.iter_mut().enumerate().take(2) {
                            tangent[a] += sh.dn[i][b] * (node[a] + z * n[a]);
                        }
                        g[2][a] += 0.1 * sh.n[i] * n[a];
                    }
                }
                let v = cross(omega, x);
                let speed: f64 = (0..3).map(|a| (v[a] + velocity[a]).powi(2)).sum();
                let dv = 3.0 * q.weight * w * dot(cross(g[0], g[1]), g[2]);
                reference += 0.5 * dv * (speed + z * z * drill * drill);
                total += dv;
                rotary += dv * z * z;
            }
        }
        assert!((energy(&m, &u) / reference - 1.0).abs() < 2e-13);
        let lump = s.lumped_mass(3.0).unwrap();
        for a in 0..3 {
            let trans = field(&xyz, |_| {
                let mut v = [0.0; 6];
                v[a] = 1.0;
                v
            });
            assert!((2.0 * energy(&m, &trans) / total - 1.0).abs() < 2e-13);
            assert!((2.0 * energy(&lump, &trans) / total - 1.0).abs() < 2e-13);
            let rot = field(&xyz, |_| {
                let mut v = [0.0; 6];
                v[a + 3] = 1.0;
                v
            });
            assert!((2.0 * energy(&lump, &rot) / rotary - 1.0).abs() < 2e-13);
            let f = s.gravity(3.0, std::array::from_fn(|c| if c == a { 1.0 } else { 0.0 })).unwrap();
            assert!(m.mul_vec(&trans).iter().zip(f).all(|(x, y)| (x - y).abs() < 2e-14));
        }
        assert!(lump.data.iter().enumerate().all(|(i, x)| if i / 48 == i % 48 { *x > 0.0 } else { *x == 0.0 }));
    }
    use fieldcraft_fem::laminate::PlyProperties;
    let plies: Vec<_> = [(0.08, 2.0), (0.12, 5.0)]
        .into_iter()
        .map(|(thickness, density)| PlyProperties {
            thickness,
            density: Some(density),
            elastic: Elastic::Isotropic { young: 1200.0, poisson: 0.25 },
            angle: 0.0,
            points: 3,
            strength: None,
        })
        .collect();
    let s = Mitc8::laminated(&nodes(false), &plies, 0.0, None, 1).unwrap();
    let m = s.mass(1.0).unwrap();
    let t = field(&nodes(false), |_| [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let r = field(&nodes(false), |_| [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    // area=4; I0=.76; I1=.0144; I2=2*(.1³-.02³)/3+5*(.1³+.02³)/3.
    assert!((2.0 * energy(&m, &t) - 4.0 * 0.76).abs() < 1e-13);
    assert!((t.iter().zip(m.mul_vec(&r)).map(|(a, b)| a * b).sum::<f64>() - 4.0 * 0.0144).abs() < 1e-13);
    let i2 = 2.0 * (0.001 - 0.000008) / 3.0 + 5.0 * (0.001 + 0.000008) / 3.0;
    assert!((2.0 * energy(&m, &r) - 4.0 * i2).abs() < 1e-13);
    assert!(s.lumped_mass(1.0).unwrap().mul_vec(&r).iter().zip(&t).map(|(a, b)| a * b).sum::<f64>() == 0.0);
}

#[test]
fn plate_initial_stress_has_independent_membrane_energy() {
    let xyz = nodes(false);
    let s = Mitc8::new(&xyz, section(0.2, 0.0)).unwrap();
    // epsilon_x=-.001, epsilon_y=nu*.001 -> sigma_x=-1.2, sigma_y=0.
    let preload = field(&xyz, |x| [-0.001 * x[0], 0.00025 * x[1], 0.0, 0.0, 0.0, 0.0]);
    let kg = s.plate_geometric_stiffness(&preload).unwrap();
    let w = field(&xyz, |x| [0.0, 0.0, x[0] * x[1], 0.0, 0.0, 0.0]);
    assert!((energy(&kg, &w) - (-0.5 * 1.2 * 0.2 * 4.0 / 3.0)).abs() < 2e-13);
    assert!(
        Mitc8::new(&nodes(true), section(0.2, 0.0)).unwrap().plate_geometric_stiffness(&preload).unwrap_err().to_string().contains("curved-shell")
    );
}
