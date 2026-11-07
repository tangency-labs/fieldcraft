use fieldcraft_fem::dense::{Cholesky, Matrix};
use fieldcraft_fem::frame::cross;
use fieldcraft_fem::laminate::PlyProperties;
use fieldcraft_fem::shell::{Shell, ShellSection};
use fieldcraft_model::{Elastic, ElementKind, MaterialOrientation, PlyStrength};

const QUAD: [[f64; 3]; 4] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
const ORTHO: Elastic = Elastic::EngineeringConstants { young: [100.0, 10.0, 10.0], poisson: [0.0; 3], shear: [5.0, 4.0, 3.0] };

fn ply(elastic: Elastic, angle: f64, thickness: f64, density: f64) -> PlyProperties {
    PlyProperties { elastic, angle, thickness, density: Some(density), strength: None, points: 3 }
}

fn close_matrices(a: &Matrix, b: &Matrix) {
    let scale = a.data.iter().map(|v| v.abs()).fold(0.0, f64::max);
    assert_eq!(a.n, b.n);
    for (x, y) in a.data.iter().zip(&b.data) {
        assert!((x - y).abs() < 2e-13 * scale, "{x} != {y}");
    }
}

#[test]
fn isotropic_plies_reproduce_the_homogeneous_mitc4_and_mitc3() {
    let elastic = Elastic::Isotropic { young: 1000.0, poisson: 0.25 };
    let plies = [ply(elastic, 0.2, 0.05, 3.0), ply(elastic, -0.3, 0.15, 3.0)];
    for (kind, xyz) in [(ElementKind::Quad4, &QUAD[..]), (ElementKind::Tri3, &QUAD[..3])] {
        for offset in [0.0, -0.5, 0.25] {
            let a = Shell::new(kind, xyz, ShellSection::new(0.2, offset, &elastic).unwrap()).unwrap();
            let b = Shell::laminated(kind, xyz, &plies, offset, None, 1).unwrap();
            close_matrices(&a.stiffness().unwrap(), &b.stiffness().unwrap());
            close_matrices(&a.mass(3.0).unwrap(), &b.ply_mass().unwrap());
            let al = a.lumped_mass(3.0).unwrap();
            let bl = b.ply_lumped_mass().unwrap();
            for (x, y) in al.iter().zip(bl) {
                assert!((x - y).abs() < 1e-13);
            }
            let force = a.body_force([0.0, 0.0, -2.0]).unwrap();
            let other = b.body_force([0.0, 0.0, -2.0]).unwrap();
            for (x, y) in force.iter().zip(other) {
                assert!((x - y).abs() < 1e-13);
            }
        }
    }
}

#[test]
fn rotated_laminates_have_six_rigid_motions_and_preserve_density_eccentricity() {
    let plies = [ply(ORTHO, 0.2, 0.25, 2.0), ply(ORTHO, -0.7, 0.75, 7.0)];
    let tilted = QUAD.map(|p| [p[0], 0.6 * p[1], 0.8 * p[1]]);
    for (kind, xyz) in [(ElementKind::Quad4, &tilted[..]), (ElementKind::Tri3, &tilted[..3])] {
        let shell = Shell::laminated(kind, xyz, &plies, 0.3, None, 1).unwrap();
        let k = shell.stiffness().unwrap();
        let scale = k.data.iter().map(|v| v.abs()).fold(0.0, f64::max);
        for i in 0..k.n {
            for j in 0..k.n {
                assert_eq!(k.get(i, j), k.get(j, i));
            }
        }
        for dof in 0..6 {
            let rigid: Vec<_> = xyz
                .iter()
                .flat_map(|p| {
                    let mut u = [0.0; 6];
                    if dof < 3 {
                        u[dof] = 1.0;
                    } else {
                        let mut w = [0.0; 3];
                        w[dof - 3] = 1.0;
                        u[..3].copy_from_slice(&cross(w, *p));
                        u[3..].copy_from_slice(&w);
                    }
                    u
                })
                .collect();
            assert!(k.mul_vec(&rigid).iter().all(|f| f.abs() < 1e-13 * scale));
        }
        Cholesky::new(&shell.ply_mass().unwrap()).unwrap();
        assert!(shell.mass(1.0).is_err(), "cannot silently replace per-ply densities by one density");
    }
    let flat = Shell::laminated(ElementKind::Quad4, &QUAD, &plies, 0.0, None, 1).unwrap();
    let gravity = flat.ply_gravity([1.0, 0.0, 0.0]).unwrap();
    assert!((gravity.iter().step_by(6).sum::<f64>() - 5.75).abs() < 1e-13);
    assert!((gravity.iter().skip(4).step_by(6).sum::<f64>() - 0.46875).abs() < 1e-13);
}

#[test]
fn a_fixed_material_frame_is_independent_of_the_first_element_edge() {
    let orientation = MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [0.0, 1.0, 0.0] };
    let plies = [ply(ORTHO, 0.3, 0.1, 2.0), ply(ORTHO, -0.6, 0.2, 3.0)];
    let a = Shell::laminated(ElementKind::Quad4, &QUAD, &plies, 0.2, Some(&orientation), 1).unwrap();
    let rotated = [QUAD[1], QUAD[2], QUAD[3], QUAD[0]];
    let b = Shell::laminated(ElementKind::Quad4, &rotated, &plies, 0.2, Some(&orientation), 1).unwrap();
    let ka = a.stiffness().unwrap();
    let kb = b.stiffness().unwrap();
    for i in 0..24 {
        for j in 0..24 {
            let ri = ((i / 6 + 3) % 4) * 6 + i % 6;
            let rj = ((j / 6 + 3) % 4) * 6 + j % 6;
            assert!((ka.get(i, j) - kb.get(ri, rj)).abs() < 1e-12, "{i} {j}");
        }
    }
}

#[test]
fn a_large_finite_ply_angle_does_not_erase_the_material_frame() {
    let frame = MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [0.0, 1.0, 0.0] };
    let plies = [ply(ORTHO, 1e20, 0.1, 2.0)];
    let a = Shell::laminated(ElementKind::Quad4, &QUAD, &plies, 0.0, Some(&frame), 1).unwrap();
    let shifted = [QUAD[1], QUAD[2], QUAD[3], QUAD[0]];
    let b = Shell::laminated(ElementKind::Quad4, &shifted, &plies, 0.0, Some(&frame), 1).unwrap();
    let ka = a.stiffness().unwrap();
    let kb = b.stiffness().unwrap();
    for i in 0..24 {
        for j in 0..24 {
            let ri = ((i / 6 + 3) % 4) * 6 + i % 6;
            let rj = ((j / 6 + 3) % 4) * 6 + j % 6;
            assert!((ka.get(i, j) - kb.get(ri, rj)).abs() < 1e-12);
        }
    }
}

#[test]
fn nafems_r0093_one_element_laminate_stress() {
    // Original model of NAFEMS R0093/1, Hopkins (2005), from physical description and
    // benchmark columns in the COMSOL 6.3 public verification page (Tables 1–3,8):
    // https://doc.comsol.com/6.3/doc/com.comsol.help.models.compmat.failure_prediction_in_a_laminated_composite_shell/failure_prediction_in_a_laminated_composite_shell.html
    // SI internally. No vendor input file. Fixed acceptance before this run: 1% per
    // component plus 0.0005 MPa rounding, allowing differing linear shell formulations.
    let elastic = Elastic::EngineeringConstants { young: [207e9, 7.6e9, 7.6e9], poisson: [0.3, 0.3, 0.0], shear: [5e9, 5e9, 3.8e9] };
    let strength = PlyStrength { tension: [500e6, 5e6], compression: [350e6, 75e6], shear12: 35e6, shear23: 35e6, tsai_wu_interaction: -0.5 };
    let plies =
        [90.0_f64, -45.0, 45.0, 0.0].map(|angle| PlyProperties { strength: Some(strength), ..ply(elastic, angle.to_radians(), 0.00005, 1500.0) });
    let xyz = QUAD.map(|p| p.map(|v| v * 0.01));
    let shell = Shell::laminated(ElementKind::Quad4, &xyz, &plies, 0.0, None, 1).unwrap();
    let held = [0, 1, 2, 3, 4, 5, 11, 17, 18, 23];
    let free: Vec<_> = (0..24).filter(|i| !held.contains(i)).collect();
    let k = shell.stiffness().unwrap();
    let mut reduced = Matrix::zeros(free.len());
    for (i, &a) in free.iter().enumerate() {
        for (j, &b) in free.iter().enumerate() {
            reduced.add(i, j, k.get(a, b));
        }
    }
    let mut rhs: Vec<_> = free.iter().map(|i| if [6, 12].contains(i) { 7.5 } else { 0.0 }).collect();
    Cholesky::new(&reduced).unwrap().solve(&mut rhs);
    let mut u = vec![0.0; 24];
    for (&i, &value) in free.iter().zip(&rhs) {
        u[i] = value;
    }
    let fields = shell.stresses(&u).unwrap();
    let expected = [[-5.128, 4.407, -1.663], [12.59, 1.983, 2.572], [8.520, 0.1256, -2.051], [9.357, -1.859, -0.5557]];
    for (ply, reference) in expected.iter().enumerate() {
        for (j, &component) in [0, 1, 3].iter().enumerate() {
            let actual = fields.plies.iter().map(|p| p[ply].middle.stress[component] / 1e6 / 4.0).sum::<f64>();
            println!("NAFEMS R0093/1 ply={} component={component}, stress_MPa={actual:.9}, reference={}", ply + 1, reference[j]);
            assert!((actual - reference[j]).abs() <= 0.01 * reference[j].abs() + 0.0005, "ply {ply}: {actual} != {}", reference[j]);
        }
    }
}

#[test]
fn flat_plate_initial_stress_stiffness_has_the_membrane_energy_and_compression_sign() {
    let plies = [ply(ORTHO, 0.0, 0.1, 2.0)];
    for (kind, xyz) in [(ElementKind::Quad4, &QUAD[..]), (ElementKind::Tri3, &QUAD[..3])] {
        let shell = Shell::laminated(kind, xyz, &plies, 0.0, None, 1).unwrap();
        let u: Vec<_> = xyz.iter().flat_map(|p| [-0.01 * p[0], 0.0, 0.0, 0.0, 0.0, 0.0]).collect();
        let w: Vec<_> = xyz.iter().flat_map(|p| [0.0, 0.0, p[0], 0.0, 0.0, 0.0]).collect();
        let kg = shell.plate_geometric_stiffness(&u).unwrap();
        let work: f64 = w.iter().zip(kg.mul_vec(&w)).map(|(x, y)| x * y).sum();
        assert!((work - (-0.01 * 100.0 * 0.1 * shell.area())).abs() < 1e-14);
        assert!(shell.plate_geometric_stiffness(&[]).is_err());
    }
}

fn navier_amplitude(h: f64) -> (f64, f64) {
    // Independent cross-ply closed forms: Q0=diag(100,10,5) GPa and outer plies at 0°.
    // FSDT modal amplitudes (φx,φy,w); hard simply supported edges. α=β=π.
    let a = std::f64::consts::PI;
    let d11 = 710.0 / 96.0 * 1e9 * h.powi(3);
    let d22 = 170.0 / 96.0 * 1e9 * h.powi(3);
    let d66 = 5e9 * h.powi(3) / 12.0;
    let s = (5.0 / 6.0) * 3.5e9 * h;
    let (k11, k22, k12) = ((d11 + d66) * a * a + s, (d22 + d66) * a * a + s, d66 * a * a);
    let det = k11 * k22 - k12 * k12;
    let mindlin = 1.0 / (2.0 * s * a * a - (s * a).powi(2) * (k11 + k22 - 2.0 * k12) / det);
    let kirchhoff = 1.0 / ((d11 + 4.0 * d66 + d22) * a.powi(4));
    (mindlin, kirchhoff)
}

fn sine_plate(n: usize, h: f64) -> f64 {
    let elastic = Elastic::EngineeringConstants { young: [100e9, 10e9, 10e9], poisson: [0.0; 3], shear: [5e9, 4e9, 3e9] };
    let plies = [0.0, std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2, 0.0].map(|angle| ply(elastic, angle, h / 4.0, 1600.0));
    let nodes = (n + 1) * (n + 1);
    let mut free = vec![None; 6 * nodes];
    let mut count = 0;
    for y in 0..=n {
        for x in 0..=n {
            let node = y * (n + 1) + x;
            for (dof, held) in [(2, x == 0 || x == n || y == 0 || y == n), (3, x == 0 || x == n), (4, y == 0 || y == n)] {
                if !held {
                    free[6 * node + dof] = Some(count);
                    count += 1;
                }
            }
        }
    }
    let mut k = Matrix::zeros(count);
    let mut f = vec![0.0; count];
    for y in 0..n {
        for x in 0..n {
            let labels = [y * (n + 1) + x, y * (n + 1) + x + 1, (y + 1) * (n + 1) + x + 1, (y + 1) * (n + 1) + x];
            let xyz = labels.map(|node| [(node % (n + 1)) as f64 / n as f64, (node / (n + 1)) as f64 / n as f64, 0.0]);
            let shell = Shell::laminated(ElementKind::Quad4, &xyz, &plies, 0.0, None, 1).unwrap();
            let local = shell.stiffness().unwrap();
            for i in 0..24 {
                for j in 0..24 {
                    if let (Some(gi), Some(gj)) = (free[6 * labels[i / 6] + i % 6], free[6 * labels[j / 6] + j % 6]) {
                        k.add(gi, gj, local.get(i, j));
                    }
                }
            }
            // Consistent integration of our sinusoidal load, 4x4 Gauss, independently of shell.
            for &(r, wr) in fieldcraft_fem::quadrature::gauss_legendre(4) {
                for &(s, ws) in fieldcraft_fem::quadrature::gauss_legendre(4) {
                    let px = (x as f64 + (r + 1.0) / 2.0) / n as f64;
                    let py = (y as f64 + (s + 1.0) / 2.0) / n as f64;
                    let q = (std::f64::consts::PI * px).sin() * (std::f64::consts::PI * py).sin();
                    let shape = [(1.0 - r) * (1.0 - s), (1.0 + r) * (1.0 - s), (1.0 + r) * (1.0 + s), (1.0 - r) * (1.0 + s)];
                    for a in 0..4 {
                        if let Some(i) = free[6 * labels[a] + 2] {
                            f[i] += q * shape[a] * wr * ws / (16.0 * (n * n) as f64);
                        }
                    }
                }
            }
        }
    }
    Cholesky::new(&k).unwrap().solve(&mut f);
    f[free[6 * ((n / 2) * (n + 1) + n / 2) + 2].unwrap()]
}

#[test]
fn thin_and_thick_cross_ply_plates_converge_to_navier() {
    for h in [0.01, 0.1] {
        let (mindlin, kirchhoff) = navier_amplitude(h);
        let coarse = sine_plate(4, h);
        let fine = sine_plate(8, h);
        let coarse_error = (coarse / mindlin - 1.0).abs();
        let fine_error = (fine / mindlin - 1.0).abs();
        println!("Navier h={h}: W4={coarse:.12e}, W8={fine:.12e}, FSDT={mindlin:.12e}, CLT={kirchhoff:.12e}, error8={fine_error:.6}");
        assert!(fine_error < 0.03, "{h}: {fine_error}");
        assert!(fine_error < 0.35 * coarse_error, "second-order refinement: {coarse_error} -> {fine_error}");
        if h == 0.01 {
            assert!((mindlin / kirchhoff - 1.0).abs() < 0.005);
        } else {
            assert!(mindlin / kirchhoff > 1.1, "thick-plate shear cannot be silently omitted");
        }
    }
}

#[test]
fn laminated_shell_recovery_is_objective_and_rejects_nonfinite_outputs() {
    let plies = [ply(ORTHO, 0.4, 0.25, 2.0), ply(ORTHO, -0.7, 0.75, 7.0)];
    let flat = Shell::laminated(ElementKind::Quad4, &QUAD, &plies, 0.0, None, 1).unwrap();
    let axes = [[1.0, 0.0, 0.0], [0.0, 0.6, 0.8], [0.0, -0.8, 0.6]];
    let rotated = QUAD.map(|v| fieldcraft_fem::frame::to_global(&axes, v));
    let tilt = Shell::laminated(ElementKind::Quad4, &rotated, &plies, 0.0, None, 1).unwrap();
    let u: Vec<_> = QUAD.iter().flat_map(|p| [0.001 * p[0], 0.002 * p[1], 0.0, 0.0, 0.0, 0.0]).collect();
    let ru = fieldcraft_fem::frame::rotate_to_global(&u, &axes);
    let a = flat.stresses(&u).unwrap();
    let b = tilt.stresses(&ru).unwrap();
    for (a, b) in a.plies.iter().zip(b.plies) {
        for (a, b) in a.iter().zip(b) {
            for (x, y) in a.middle.stress.iter().zip(b.middle.stress) {
                assert!((x - y).abs() < 1e-14);
            }
        }
    }
    let a = flat.ply_gravity([0.0, 0.0, -1.0]).unwrap();
    let b = tilt.ply_gravity(fieldcraft_fem::frame::to_global(&axes, [0.0, 0.0, -1.0])).unwrap();
    for (x, y) in fieldcraft_fem::frame::rotate_to_global(&a, &axes).iter().zip(b) {
        assert!((x - y).abs() < 1e-13);
    }
    assert!(flat.ply_gravity([f64::NAN, 0.0, 0.0]).is_err());
    assert!(flat.body_force([f64::INFINITY, 0.0, 0.0]).is_err());
    assert!(flat.stresses(&[f64::NAN; 24]).is_err());
    let elastic = Elastic::Isotropic { young: 1e300, poisson: 0.0 };
    let huge = QUAD.map(|v| v.map(|x| x * 1e8));
    let s = Shell::laminated(ElementKind::Quad4, &huge, &[ply(elastic, 0.0, 1.0, 1e300)], 0.0, None, 1).unwrap();
    assert!(s.ply_mass().is_err());
    assert!(s.ply_lumped_mass().is_err());
}
