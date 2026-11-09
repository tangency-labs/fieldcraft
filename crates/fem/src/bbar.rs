//! Volume-averaged dilatation (B-bar), small-strain Hex8.
//!
//! T. J. R. Hughes, Generalization of selective integration procedures to
//! anisotropic and nonlinear media, IJNME 15 (1980), 1413–1418,
//! doi:10.1002/nme.1620150914. This implementation is restricted at the model
//! boundary to linear isotropic elasticity, and uses volume-weighted 2x2x2
//! quadrature, not the gradient evaluated at the element centre.
use crate::{
    Result,
    dense::Matrix,
    element::{Solid, b_block},
    material::Elasticity,
};

impl Solid<'_> {
    pub(crate) fn bbar_average(&self) -> Result<Vec<[f64; 3]>> {
        let mut average = vec![[0.0; 3]; 8];
        let mut volume = 0.0;
        for p in self.integration_points() {
            let g = self.gradients(p.xi)?;
            let w = g.det * p.weight;
            volume += w;
            for (a, v) in average.iter_mut().zip(g.dx.iter()) {
                for j in 0..3 {
                    a[j] += w * v[j];
                }
            }
        }
        if !volume.is_finite() || volume <= 0.0 {
            return Err(crate::Error::Degenerate("B-bar volume must be finite and positive".into()));
        }
        for a in &mut average {
            for v in a {
                *v /= volume;
            }
        }
        Ok(average)
    }

    pub(crate) fn bbar_stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        let average = self.bbar_average()?;
        let mut k = Matrix::zeros(24);
        for p in self.integration_points() {
            let g = self.gradients(p.xi)?;
            let blocks: Vec<_> =
                g.dx.iter()
                    .zip(&average)
                    .map(|(gradient, avg)| {
                        let mut b = b_block(gradient);
                        for row in &mut b[..3] {
                            for j in 0..3 {
                                row[j] += (avg[j] - gradient[j]) / 3.0;
                            }
                        }
                        b
                    })
                    .collect();
            for (i, bi) in blocks.iter().enumerate() {
                for (j, bj) in blocks.iter().enumerate() {
                    for a in 0..3 {
                        for b in 0..3 {
                            let mut value = 0.0;
                            for r in 0..6 {
                                for c in 0..6 {
                                    value += bi[r][a] * d[r][c] * bj[c][b];
                                }
                            }
                            k.add(3 * i + a, 3 * j + b, value * g.det * p.weight);
                        }
                    }
                }
            }
        }
        Ok(k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{element::Formulation, material::isotropic};
    use fieldcraft_model::ElementKind;

    #[test]
    fn bbar_independent_volumetric_projection_rigid_and_affine_patch() {
        for warped in [false, true] {
            let mut xyz = vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
            if warped {
                xyz[6] = [1.2, 0.9, 1.1];
                xyz[4] = [-0.1, 0., 1.];
            }
            let bbar = Solid::new(ElementKind::Hex8, &xyz, Formulation::BBar).unwrap();
            let standard = Solid::new(ElementKind::Hex8, &xyz, Formulation::Full).unwrap();
            let d = isotropic(1000., 0.3).unwrap();
            let k = bbar.stiffness(&d).unwrap();
            let full = standard.stiffness(&d).unwrap();
            assert_eq!(bbar.mass(1.0).unwrap(), standard.mass(1.0).unwrap());
            assert_eq!(bbar.lumped_mass(1.0).unwrap(), standard.lumped_mass(1.0).unwrap());
            let free: Vec<_> = (0..24).filter(|i| ![0, 1, 2, 4, 5, 11].contains(i)).collect();
            let mut restrained = Matrix::zeros(free.len());
            for (i, &a) in free.iter().enumerate() {
                for (j, &b) in free.iter().enumerate() {
                    restrained.add(i, j, k.get(a, b));
                }
            }
            crate::dense::Cholesky::new(&restrained).unwrap();
            // Independent isotropic volumetric decomposition: remove local bulk
            // divergence energy from full K, replace by its rank-one mean energy.
            let bulk = 1000. / (3. * (1. - 2. * 0.3));
            let mut integral = [0.; 24];
            let mut local = Matrix::zeros(24);
            let mut volume = 0.;
            for p in standard.integration_points() {
                let g = standard.gradients(p.xi).unwrap();
                let w = p.weight * g.det;
                volume += w;
                let div: Vec<_> = g.dx.iter().take(8).flatten().copied().collect();
                for i in 0..24 {
                    integral[i] += w * div[i];
                    for j in 0..24 {
                        local.add(i, j, w * div[i] * div[j]);
                    }
                }
            }
            for i in 0..24 {
                for j in 0..24 {
                    let expected = full.get(i, j) - bulk * (local.get(i, j) - integral[i] * integral[j] / volume);
                    assert!((k.get(i, j) - expected).abs() < 1e-10);
                    assert!((k.get(i, j) - k.get(j, i)).abs() < 1e-10);
                }
            }
            // Three translations and three infinitesimal rotations are null modes.
            for mode in 0..6 {
                let u: Vec<_> = xyz
                    .iter()
                    .flat_map(|x| match mode {
                        0 => [1., 0., 0.],
                        1 => [0., 1., 0.],
                        2 => [0., 0., 1.],
                        3 => [0., -x[2], x[1]],
                        4 => [x[2], 0., -x[0]],
                        _ => [-x[1], x[0], 0.],
                    })
                    .collect();
                for i in 0..24 {
                    assert!((0..24).map(|j| k.get(i, j) * u[j]).sum::<f64>().abs() < 1e-10);
                }
            }
            let u: Vec<_> = xyz.iter().map(|x| [0.01 * x[0] + 0.02 * x[1], -0.03 * x[1], 0.04 * x[2]]).collect();
            let strain = [0.01, -0.03, 0.04, 0.02, 0., 0.];
            let target = crate::material::stress(&d, &strain);
            for p in bbar.stresses(&d, &u).unwrap() {
                for j in 0..6 {
                    assert!((p.strain[j] - strain[j]).abs() < 1e-13);
                    assert!((p.stress[j] - target[j]).abs() < 1e-10);
                }
            }
            let flat: Vec<_> = u.iter().flatten().copied().collect();
            for i in 0..24 {
                assert!((0..24).map(|j| (k.get(i, j) - full.get(i, j)) * flat[j]).sum::<f64>().abs() < 1e-10);
            }
            println!("B-bar independent projection/rigid/affine checks: warped={warped}, volume={volume}");
        }
    }
}
