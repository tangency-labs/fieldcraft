//! Thirteen-variable enhanced-gradient Hex8, infinitesimal linear elasticity.
//!
//! Simo & Rifai (1990), IJNME 29, 1595–1638, doi:10.1002/nme.1620290802,
//! sections 2.3 and 3.1: orthogonal enhanced strains and element-local elimination.
//! The explicit 3D basis is documented in the public Abaqus Theory Guide,
//! "Continuum elements with incompatible modes", geometrically linear formulation:
//! https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-incompatible.htm
//! Nine columns are sym(e_c ⊗ (j0/j * xi_k * J0^{-T} e_k)); four are
//! j0/j * {xi*eta, xi*zeta, eta*zeta, xi*eta*zeta} I. Engineering shear is used.
//! Full 2x2x2 integration and a Schur complement eliminate all 13 amplitudes.
//! This is an independently selectable native formulation, not a claim of deck equivalence.
use crate::{
    Error, Result,
    dense::{Cholesky, Matrix},
    element::{PointStress, Solid, b_block},
    material::{Elasticity, stress},
};

const N: usize = 24;
const A: usize = 13;
type Columns = [[f64; 6]; A];
struct Blocks {
    uu: Matrix,
    ua: [[f64; A]; N],
    aa: Matrix,
}

impl Solid<'_> {
    fn enhanced13_columns(&self, xi: [f64; 3], det: f64, centre: &([[f64; 3]; 3], f64)) -> Columns {
        let (inverse, det0) = centre;
        let ratio = det0 / det;
        let mut columns = [[0.; 6]; A];
        for k in 0..3 {
            let gradient = [0, 1, 2].map(|b| ratio * xi[k] * inverse[b][k]);
            let b = b_block(&gradient);
            for c in 0..3 {
                for r in 0..6 {
                    columns[3 * k + c][r] = b[r][c];
                }
            }
        }
        let [x, y, z] = xi;
        for (i, value) in [x * y, x * z, y * z, x * y * z].into_iter().enumerate() {
            for r in 0..3 {
                columns[9 + i][r] = ratio * value;
            }
        }
        columns
    }

    fn enhanced13_blocks(&self, d: &Elasticity) -> Result<Blocks> {
        let mut blocks = Blocks { uu: Matrix::zeros(N), ua: [[0.; A]; N], aa: Matrix::zeros(A) };
        let centre = self.centre_jacobian()?;
        for p in self.integration_points() {
            let g = self.gradients(p.xi)?;
            let w = g.det * p.weight;
            let a = self.enhanced13_columns(p.xi, g.det, &centre);
            let mut b = [[0.; 6]; N];
            for (i, gradient) in g.dx.iter().take(8).enumerate() {
                let block = b_block(gradient);
                for c in 0..3 {
                    for r in 0..6 {
                        b[3 * i + c][r] = block[r][c];
                    }
                }
            }
            let da = a.map(|col| stress(d, &col));
            let db = b.map(|col| stress(d, &col));
            for i in 0..N {
                for j in 0..N {
                    blocks.uu.add(i, j, dot(&b[i], &db[j]) * w);
                }
                for j in 0..A {
                    blocks.ua[i][j] += dot(&b[i], &da[j]) * w;
                }
            }
            for i in 0..A {
                for j in 0..A {
                    blocks.aa.add(i, j, dot(&a[i], &da[j]) * w);
                }
            }
        }
        Ok(blocks)
    }

    pub(crate) fn enhanced13_stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        let Blocks { mut uu, ua, aa } = self.enhanced13_blocks(d)?;
        let chol = Cholesky::new(&aa)?;
        for i in 0..N {
            let mut solution = ua[i];
            chol.solve(&mut solution);
            for j in 0..=i {
                let correction = ua[j].iter().zip(solution).map(|(a, b)| a * b).sum::<f64>();
                uu.add(j, i, -correction);
                if i != j {
                    uu.add(i, j, -correction);
                }
            }
        }
        Ok(uu)
    }

    pub(crate) fn enhanced13_stresses(&self, d: &Elasticity, u: &[[f64; 3]]) -> Result<Vec<PointStress>> {
        if u.len() != 8 || u.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("enhanced_strain13 requires eight finite displacement vectors".into()));
        }
        let Blocks { ua, aa, .. } = self.enhanced13_blocks(d)?;
        let mut alpha = [0.; A];
        for (row, value) in ua.iter().zip(u.iter().flatten()) {
            for (a, coefficient) in alpha.iter_mut().zip(row) {
                *a -= coefficient * value;
            }
        }
        Cholesky::new(&aa)?.solve(&mut alpha);
        let centre = self.centre_jacobian()?;
        let mut result = Vec::new();
        for p in self.integration_points() {
            let g = self.gradients(p.xi)?;
            let mut strain = [0.; 6];
            for (gradient, displacement) in g.dx.iter().zip(u) {
                let b = b_block(gradient);
                for r in 0..6 {
                    strain[r] += b[r].iter().zip(displacement).map(|(a, v)| a * v).sum::<f64>();
                }
            }
            for (column, value) in self.enhanced13_columns(p.xi, g.det, &centre).iter().zip(alpha) {
                for r in 0..6 {
                    strain[r] += column[r] * value;
                }
            }
            let mut xyz = [0.; 3];
            for (n, x) in g.shape.values().iter().zip(self.xyz) {
                for r in 0..3 {
                    xyz[r] += n * x[r];
                }
            }
            result.push(PointStress { xyz, strain, stress: stress(d, &strain) });
        }
        Ok(result)
    }
}
fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{element::Formulation, material::isotropic};
    use fieldcraft_model::ElementKind;

    fn geometry(warp: bool) -> Vec<[f64; 3]> {
        let mut x = vec![[0., 0., 0.], [2., 0., 0.], [2., 1., 0.], [0., 1., 0.], [0., 0., 3.], [2., 0., 3.], [2., 1., 3.], [0., 1., 3.]];
        for p in &mut x {
            p[0] += 0.2 * p[2];
            p[1] += 0.1 * p[0];
        }
        if warp {
            x[6][0] += 0.25;
            x[4][1] -= 0.13;
        }
        x
    }
    #[test]
    fn enhanced13_zero_mean_rigid_modes_and_affine_patch() {
        for warped in [false, true] {
            let xyz = geometry(warped);
            let s = Solid::new(ElementKind::Hex8, &xyz, Formulation::EnhancedStrain13).unwrap();
            assert_eq!(s.dofs().internal.incompatible_modes, 13);
            let centre = s.centre_jacobian().unwrap();
            let mut integral = [[0.; 6]; A];
            for p in s.integration_points() {
                let g = s.gradients(p.xi).unwrap();
                for (row, col) in integral.iter_mut().zip(s.enhanced13_columns(p.xi, g.det, &centre)) {
                    for r in 0..6 {
                        row[r] += col[r] * g.det * p.weight;
                    }
                }
            }
            assert!(integral.iter().flatten().all(|v| v.abs() < 2e-14));
            for nu in [0., 0.3, 0.45, 0.499] {
                let d = isotropic(2e5, nu).unwrap();
                let k = s.stiffness(&d).unwrap();
                let full = Solid::new(ElementKind::Hex8, &xyz, Formulation::Full).unwrap().stiffness(&d).unwrap();
                for mode in 0..7 {
                    let u: Vec<_> = xyz
                        .iter()
                        .map(|&[x, y, z]| match mode {
                            0 => [1., 0., 0.],
                            1 => [0., 1., 0.],
                            2 => [0., 0., 1.],
                            3 => [0., -z, y],
                            4 => [z, 0., -x],
                            5 => [-y, x, 0.],
                            _ => [0.01 * x + 0.02 * y, -0.03 * y + 0.01 * z, 0.04 * z - 0.01 * x],
                        })
                        .collect();
                    let flat: Vec<_> = u.iter().flatten().copied().collect();
                    for i in 0..N {
                        let force = (0..N).map(|j| k.get(i, j) * flat[j]).sum::<f64>();
                        let expected = if mode < 6 { 0. } else { (0..N).map(|j| full.get(i, j) * flat[j]).sum() };
                        assert!((force - expected).abs() < 2e-7, "warp={warped} nu={nu} mode={mode} force={force} expected={expected}");
                    }
                    let expected = if mode < 6 { [0.; 6] } else { [0.01, -0.03, 0.04, 0.02, -0.01, 0.01] };
                    for p in s.stresses(&d, &u).unwrap() {
                        for (a, b) in p.strain.iter().zip(expected) {
                            assert!((a - b).abs() < 2e-12);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn enhanced13_condensation_energy_and_recovered_internal_equilibrium() {
        let xyz = geometry(true);
        let s = Solid::new(ElementKind::Hex8, &xyz, Formulation::EnhancedStrain13).unwrap();
        let d = isotropic(70000., 0.3).unwrap();
        let u: Vec<_> = (0..8).map(|i| [((i * 7 + 1) as f64).sin() * 0.01, ((i * 3 + 2) as f64).cos() * 0.02, (i as f64 - 4.) * 0.003]).collect();
        let flat: Vec<_> = u.iter().flatten().copied().collect();
        let k = s.stiffness(&d).unwrap();
        let energy: f64 = (0..N).map(|i| (0..N).map(|j| flat[i] * k.get(i, j) * flat[j]).sum::<f64>()).sum();
        let mut integrated = 0.;
        let mut residual = [0.; A];
        let centre = s.centre_jacobian().unwrap();
        for (p, recovered) in s.integration_points().iter().zip(s.stresses(&d, &u).unwrap()) {
            let g = s.gradients(p.xi).unwrap();
            let w = g.det * p.weight;
            integrated += dot(&recovered.strain, &recovered.stress) * w;
            for (r, col) in residual.iter_mut().zip(s.enhanced13_columns(p.xi, g.det, &centre)) {
                *r += dot(&col, &recovered.stress) * w;
            }
        }
        assert!((energy - integrated).abs() < energy.abs() * 2e-12);
        assert!(residual.iter().all(|r| r.abs() < 2e-10));
        assert!(energy > 0.);
    }
    #[test]
    fn enhanced13_rectangular_pure_bending_matches_continuum_energy() {
        let xyz = vec![[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.], [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]];
        for nu in [0., 0.3, 0.45] {
            let e = 2e5;
            let curvature = 0.02;
            let s = Solid::new(ElementKind::Hex8, &xyz, Formulation::EnhancedStrain13).unwrap();
            let d = isotropic(e, nu).unwrap();
            let u: Vec<_> = xyz
                .iter()
                .map(|&[x, y, z]| [-curvature * x * z, nu * curvature * y * z, 0.5 * curvature * (x * x + nu * z * z - nu * y * y)])
                .collect();
            let flat: Vec<_> = u.iter().flatten().copied().collect();
            let k = s.stiffness(&d).unwrap();
            let energy: f64 = flat.iter().zip(k.mul_vec(&flat)).map(|(a, b)| a * b).sum();
            assert!((energy - e * curvature * curvature * 8. / 3.).abs() < 1e-9);
            for point in s.stresses(&d, &u).unwrap() {
                assert!((point.stress[0] + e * curvature * point.xyz[2]).abs() < 1e-9);
                assert!(point.stress[1..].iter().all(|v| v.abs() < 1e-9));
            }
        }
    }

    #[test]
    fn enhanced13_no_hourglass_and_rotation_covariance_on_distorted_hex() {
        let x = geometry(true);
        let d = isotropic(2e5, 0.3).unwrap();
        let angle = 0.63_f64;
        let (sn, cs) = angle.sin_cos();
        let rotate = |[x, y, z]: [f64; 3]| [cs * x - sn * y, sn * x + cs * y, z];
        let xr: Vec<_> = x.iter().copied().map(rotate).collect();
        let s = Solid::new(ElementKind::Hex8, &x, Formulation::EnhancedStrain13).unwrap();
        let sr = Solid::new(ElementKind::Hex8, &xr, Formulation::EnhancedStrain13).unwrap();
        let k = s.stiffness(&d).unwrap();
        let kr = sr.stiffness(&d).unwrap();
        let free: Vec<_> = (0..24).filter(|i| ![0, 1, 2, 4, 5, 11].contains(i)).collect();
        let mut constrained = Matrix::zeros(18);
        for (i, &a) in free.iter().enumerate() {
            for (j, &b) in free.iter().enumerate() {
                constrained.add(i, j, k.get(a, b));
            }
        }
        Cholesky::new(&constrained).unwrap();
        for column in 0..24 {
            let mut u = [[0.; 3]; 8];
            u[column / 3][column % 3] = 1.;
            let ur: Vec<_> = u.iter().copied().map(rotate).flat_map(|p| p.into_iter()).collect();
            let f = k.mul_vec(&u.iter().flatten().copied().collect::<Vec<_>>());
            let fr = kr.mul_vec(&ur);
            for (a, b) in f.as_chunks::<3>().0.iter().zip(fr.as_chunks::<3>().0) {
                for (v, w) in rotate([a[0], a[1], a[2]]).into_iter().zip(b) {
                    assert!((v - w).abs() < 2e-9);
                }
            }
        }
    }
    #[test]
    fn enhanced13_rejects_other_element_and_combined_formulations() {
        use fieldcraft_model::{Formulation as ModelFormulation, Integration, SectionKind};
        let base = ModelFormulation { enhanced_strain13: true, ..Default::default() };
        assert!(Formulation::from_section(ElementKind::Tet4, &base).is_err());
        assert!(base.validate_for(&SectionKind::PlaneStrain { thickness: 1. }).is_err());
        for f in [
            ModelFormulation { b_bar: true, ..base },
            ModelFormulation { hybrid: true, ..base },
            ModelFormulation { incompatible_modes: true, ..base },
            ModelFormulation { integration: Integration::Reduced, ..base },
        ] {
            assert!(f.validate_for(&SectionKind::Solid {}).is_err());
            assert!(Formulation::from_section(ElementKind::Hex8, &f).is_err());
        }
        let x = geometry(false);
        let solid = Solid::new(ElementKind::Hex8, &x, Formulation::EnhancedStrain13).unwrap();
        let d = isotropic(1000., 0.3).unwrap();
        assert!(solid.initial_stress_load(&d, &[[0.; 6]; 8]).is_err());
        assert!(solid.stresses_from(&d, &[[0.; 3]; 8], None, Some(&[[0.; 6]; 8])).is_err());
    }
}
