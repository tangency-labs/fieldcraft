//! Linear eight-node solid-shell: Harnau, Schweizerhof & Hauptmann, ECCOMAS 2000,
//! “On solid-shell elements with linear and quadratic shape functions ...”,
//! https://congress2.cimne.com/eccomas/proceedings/eccomas2000/pdf/826.pdf
//! Equations (1),(3): paired lower/upper faces, 24 translations; (5): shear ANS;
//! (6),(15): six membrane EAS modes; (10): four thickness EAS modes; (16):
//! thickness ANS. This initial specialization admits orthogonal affine sweeps
//! only, so the covariant/Cartesian transformations and Jacobian are constant.
//! Ten local strain parameters are statically condensed, with no added nodal DOF.
//! No plane-stress elimination, shear correction factor, or hourglass penalty.
use crate::dense::{Cholesky, Matrix};
use crate::element::{PointStress, Solid};
use crate::material::{Elasticity, stress};
use crate::{Error, Result};

type B = [[f64; 24]; 6];
type G = [[f64; 10]; 6];
const ASPECT_SLACK: f64 = 1e-12; // unit-factors: dimensionless aspect-ratio tolerance
const AFFINE_GEOMETRY_SLACK: f64 = 1e-9; // unit-factors: relative geometry tolerance

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn global(e: [f64; 6], axes: [[f64; 3]; 3]) -> [f64; 6] {
    let t = [[e[0], 0.5 * e[3], 0.5 * e[4]], [0.5 * e[3], e[1], 0.5 * e[5]], [0.5 * e[4], 0.5 * e[5], e[2]]];
    let mut a = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for r in 0..3 {
                for s in 0..3 {
                    a[i][j] += axes[r][i] * t[r][s] * axes[s][j];
                }
            }
        }
    }
    [a[0][0], a[1][1], a[2][2], 2.0 * a[0][1], 2.0 * a[0][2], 2.0 * a[1][2]]
}

impl Solid<'_> {
    /// Node pairs 0–4,1–5,2–6,3–7 define the positive thickness direction.
    pub(crate) fn solid_shell_axes(&self) -> Result<[[f64; 3]; 3]> {
        let x = self.xyz;
        if x.len() != 8 {
            return Err(Error::NotSupported("solid-shell requires hex8".into()));
        }
        let edges = [sub(x[1], x[0]), sub(x[3], x[0]), sub(x[4], x[0])];
        let lengths = edges.map(|e| dot(e, e).sqrt());
        if lengths.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(Error::Degenerate("solid-shell zero/nonfinite edge".into()));
        }
        let ratio = lengths.iter().copied().fold(0.0, f64::max) / lengths.iter().copied().fold(f64::INFINITY, f64::min);
        if ratio > 100.0 * (1.0 + ASPECT_SLACK) {
            return Err(Error::NotSupported(
                "solid-shell edge aspect ratio above 100: unresolved translational stiffness conditioning; refine in-plane mesh".into(),
            ));
        }
        let axes = std::array::from_fn(|i| edges[i].map(|v| v / lengths[i]));
        if dot(axes[0], axes[1]).abs() > 1e-10 || dot(axes[0], axes[2]).abs() > 1e-10 || dot(axes[1], axes[2]).abs() > 1e-10 {
            return Err(Error::NotSupported("solid-shell oblique sweep: orthogonal affine hex8 required".into()));
        }
        let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        for (i, c) in corners.iter().enumerate() {
            for a in 0..3 {
                let delta = dot(sub(x[i], x[0]), axes[a]) - c[a] as f64 * lengths[a];
                if delta.abs() > AFFINE_GEOMETRY_SLACK * lengths[a] {
                    return Err(Error::NotSupported("solid-shell curved, warped or variable-thickness hex8; affine sweep required".into()));
                }
            }
        }
        self.gradients([0.0; 3])?;
        Ok(axes)
    }
    fn ss_compatible(&self, p: [f64; 3], axes: [[f64; 3]; 3]) -> Result<B> {
        let dx = self.gradients(p)?.dx;
        let mut b = [[0.0; 24]; 6];
        for n in 0..8 {
            let d = axes.map(|axis| dot(dx[n], axis));
            for c in 0..3 {
                let u = axes.map(|axis| axis[c]);
                let j = 3 * n + c;
                b[0][j] = d[0] * u[0];
                b[1][j] = d[1] * u[1];
                b[2][j] = d[2] * u[2];
                b[3][j] = d[0] * u[1] + d[1] * u[0];
                b[4][j] = d[0] * u[2] + d[2] * u[0];
                b[5][j] = d[1] * u[2] + d[2] * u[1];
            }
        }
        Ok(b)
    }
    fn ss_operators(&self, p: [f64; 3], axes: [[f64; 3]; 3]) -> Result<(B, G)> {
        let [r, s, t] = p;
        let mut b = self.ss_compatible(p, axes)?;
        let rx0 = self.ss_compatible([0.0, -1.0, t], axes)?;
        let rx1 = self.ss_compatible([0.0, 1.0, t], axes)?;
        let sy0 = self.ss_compatible([-1.0, 0.0, t], axes)?;
        let sy1 = self.ss_compatible([1.0, 0.0, t], axes)?;
        b[2] = [0.0; 24];
        for [a, c] in [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]] {
            let sample = self.ss_compatible([a, c, 0.0], axes)?;
            let w = 0.25 * (1.0 + a * r) * (1.0 + c * s);
            for j in 0..24 {
                b[2][j] += w * sample[2][j];
            }
        }
        for j in 0..24 {
            b[4][j] = 0.5 * ((1.0 - s) * rx0[4][j] + (1.0 + s) * rx1[4][j]);
            b[5][j] = 0.5 * ((1.0 - r) * sy0[5][j] + (1.0 + r) * sy1[5][j]);
        }
        let mut g = [[0.0; 10]; 6];
        g[0][0] = r;
        g[1][1] = s;
        g[3][2] = r;
        g[3][3] = s;
        g[0][4] = r * s;
        g[1][5] = r * s;
        g[2][6] = t;
        g[2][7] = r * t;
        g[2][8] = s * t;
        g[2][9] = r * s * t;
        // Constant det(J) and basis scaling can be absorbed into local EAS parameters.
        for j in 0..24 {
            let v = global(std::array::from_fn(|i| b[i][j]), axes);
            for i in 0..6 {
                b[i][j] = v[i];
            }
        }
        for j in 0..10 {
            let v = global(std::array::from_fn(|i| g[i][j]), axes);
            for i in 0..6 {
                g[i][j] = v[i];
            }
        }
        Ok((b, g))
    }
    fn ss_blocks(&self, d: &Elasticity) -> Result<(Matrix, Vec<[f64; 24]>, Cholesky)> {
        let axes = self.solid_shell_axes()?;
        let mut k = Matrix::zeros(24);
        let mut h = Matrix::zeros(10);
        let mut l = vec![[0.0; 24]; 10];
        for q in self.rule().points() {
            let (b, g) = self.ss_operators(q.xi, axes)?;
            let w = self.gradients(q.xi)?.det * q.weight;
            for i in 0..6 {
                for j in 0..6 {
                    let v = w * d[i][j];
                    for a in 0..24 {
                        for c in 0..24 {
                            k.add(a, c, b[i][a] * v * b[j][c]);
                        }
                    }
                    for a in 0..10 {
                        for c in 0..24 {
                            l[a][c] += g[i][a] * v * b[j][c];
                        }
                        for c in 0..10 {
                            h.add(a, c, g[i][a] * v * g[j][c]);
                        }
                    }
                }
            }
        }
        Ok((k, l, Cholesky::new(&h)?))
    }
    pub(crate) fn solid_shell_stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        let (_, l, h) = self.ss_blocks(d)?;
        let axes = self.solid_shell_axes()?;
        let mut relax = [[0.0; 24]; 10];
        for j in 0..24 {
            let mut v: Vec<_> = l.iter().map(|r| r[j]).collect();
            h.solve(&mut v);
            for i in 0..10 {
                relax[i][j] = v[i];
            }
        }
        // Integrate the relaxed strain energy instead of subtracting two large
        // stiffness blocks. This is exactly Kuu-L^T H^-1 L in exact arithmetic.
        let mut k = Matrix::zeros(24);
        for q in self.rule().points() {
            let (mut b, g) = self.ss_operators(q.xi, axes)?;
            for i in 0..6 {
                for j in 0..24 {
                    for a in 0..10 {
                        b[i][j] -= g[i][a] * relax[a][j];
                    }
                }
            }
            let w = self.gradients(q.xi)?.det * q.weight;
            for a in 0..24 {
                for c in a..24 {
                    let mut v = 0.0;
                    for i in 0..6 {
                        for j in 0..6 {
                            v += b[i][a] * d[i][j] * b[j][c];
                        }
                    }
                    k.add(a, c, w * v);
                    if a != c {
                        k.add(c, a, w * v);
                    }
                }
            }
        }
        if k.data.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid("nonfinite solid-shell stiffness".into()));
        }
        Ok(k)
    }
    pub(crate) fn solid_shell_stresses(&self, d: &Elasticity, u: &[[f64; 3]]) -> Result<Vec<PointStress>> {
        if u.len() != 8 || u.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("solid-shell needs eight finite translations".into()));
        }
        let (_, l, h) = self.ss_blocks(d)?;
        let ue = u.as_flattened();
        let mut alpha: Vec<_> = l.iter().map(|r| -r.iter().zip(ue).map(|(a, b)| a * b).sum::<f64>()).collect();
        h.solve(&mut alpha);
        let axes = self.solid_shell_axes()?;
        let mut out = Vec::new();
        for q in self.rule().points() {
            let (b, g) = self.ss_operators(q.xi, axes)?;
            let strain = std::array::from_fn(|i| {
                b[i].iter().zip(ue).map(|(a, b)| a * b).sum::<f64>() + g[i].iter().zip(&alpha).map(|(a, b)| a * b).sum::<f64>()
            });
            out.push(PointStress { xyz: self.position(q.xi), strain, stress: stress(d, &strain) });
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Formulation;
    use crate::material::isotropic;
    use fieldcraft_model::ElementKind;
    fn brick(h: f64) -> [[f64; 3]; 8] {
        [
            [0., 0., -h / 2.],
            [2., 0., -h / 2.],
            [2., 1., -h / 2.],
            [0., 1., -h / 2.],
            [0., 0., h / 2.],
            [2., 0., h / 2.],
            [2., 1., h / 2.],
            [0., 1., h / 2.],
        ]
    }
    #[test]
    fn rigid_patch_bending_and_shear_across_thicknesses() {
        for h in [0.5, 0.05, 0.02] {
            let xyz = brick(h);
            let s = Solid::new(ElementKind::Hex8, &xyz, Formulation::SolidShell).unwrap();
            let d = isotropic(1e6, 0.3).unwrap();
            let k = s.stiffness(&d).unwrap();
            let rigid = xyz.map(|[x, y, z]| [0.1 + 0.2 * z - 0.3 * y, -0.2 + 0.3 * x - 0.1 * z, 0.4 + 0.1 * y - 0.2 * x]);
            for p in s.stresses(&d, &rigid).unwrap() {
                assert!(p.strain.iter().all(|v| v.abs() < 1e-11));
            }
            let u = xyz.map(|[x, y, z]| [0.01 * x + 0.003 * y, 0.02 * y + 0.004 * z, 0.03 * z + 0.005 * x]);
            for p in s.stresses(&d, &u).unwrap() {
                for (a, b) in p.strain.into_iter().zip([0.01, 0.02, 0.03, 0.003, 0.005, 0.004]) {
                    assert!((a - b).abs() < 1e-11);
                }
            }
            // Cylindrical bending: epsilon_xx=-kappa*z; epsilon_yy=0,
            // sigma_zz=0 after thickness EAS relaxation. D=E h^3/[12(1-nu^2)].
            let u = xyz.map(|[x, _, z]| [-x * z, 0.0, 0.5 * x * x]);
            let v = u.as_flattened();
            let energy = v.iter().zip(k.mul_vec(v)).map(|(a, b)| a * b).sum::<f64>();
            let expected = 1e6 * h.powi(3) / 12.0 / (1.0 - 0.3 * 0.3) * 2.0;
            println!("bending h={h} energy={energy} expected={expected}");
            assert!((energy / expected - 1.0).abs() < 2e-6, "h={h} bending {energy} vs {expected}");
            let shear = xyz.map(|[_, _, z]| [0.01 * z, 0.0, 0.0]);
            let v = shear.as_flattened();
            let energy = v.iter().zip(k.mul_vec(v)).map(|(a, b)| a * b).sum::<f64>();
            let expected = 1e6 / 2.6 * 0.01 * 0.01 * 2.0 * h;
            assert!((energy / expected - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn no_extra_zero_modes_and_hostile_geometry() {
        let x = brick(0.2);
        let s = Solid::new(ElementKind::Hex8, &x, Formulation::SolidShell).unwrap();
        let k = s.stiffness(&isotropic(1e6, 0.3).unwrap()).unwrap();
        let free: Vec<_> = (0..24).filter(|i| ![0, 1, 2, 4, 5, 11].contains(i)).collect();
        let mut a = Matrix::zeros(18);
        for (i, &r) in free.iter().enumerate() {
            for (j, &c) in free.iter().enumerate() {
                a.add(i, j, k.get(r, c));
            }
        }
        Cholesky::new(&a).unwrap();
        assert!(s.stresses(&isotropic(1.0, 0.0).unwrap(), &[[0.0; 3]; 7]).is_err());
        assert!(Solid::new(ElementKind::Hex8, &brick(0.0005), Formulation::SolidShell).is_err());
        let mut bad = x;
        bad[6][2] += 0.01;
        assert!(Solid::new(ElementKind::Hex8, &bad, Formulation::SolidShell).is_err());
        let mut bad = x;
        bad[6][2] = f64::NAN;
        assert!(Solid::new(ElementKind::Hex8, &bad, Formulation::SolidShell).is_err());
    }
    #[test]
    fn global_rotation_and_in_plane_node_order_are_physical() {
        let x = brick(0.2);
        let d = isotropic(1e6, 0.3).unwrap();
        let u = x.map(|[x, y, z]| [0.003 * x * y, 0.002 * y * z, -0.004 * x * z]);
        let energy = |x: &[[f64; 3]; 8], u: &[[f64; 3]; 8]| {
            let k = Solid::new(ElementKind::Hex8, x, Formulation::SolidShell).unwrap().stiffness(&d).unwrap();
            u.as_flattened().iter().zip(k.mul_vec(u.as_flattened())).map(|(a, b)| a * b).sum::<f64>()
        };
        let original = energy(&x, &u);
        let rotate = |[a, b, c]: [f64; 3]| [0.6 * a - 0.8 * c, b, 0.8 * a + 0.6 * c];
        assert!((energy(&x.map(rotate), &u.map(rotate)) / original - 1.0).abs() < 1e-11);
        for order in [[1, 2, 3, 0, 5, 6, 7, 4], [4, 7, 6, 5, 0, 3, 2, 1]] {
            assert!((energy(&order.map(|i| x[i]), &order.map(|i| u[i])) / original - 1.0).abs() < 1e-11);
        }
    }
}
