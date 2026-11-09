//! Two-node, small-strain axisymmetric conical shell.
//!
//! Simões Moita et al., Applied Sciences 13 (2023) 5325,
//! doi:10.3390/app13095325, equations (5),(6),(8)-(11):
//! e_s=t.u', e_h=u_r/r, k_s=-beta', k_h=-t_r beta/r,
//! gamma=n.u'-beta, t=(dr/ds,dz/ds), n=(-t_z,t_r).
//! This implementation uses two-point membrane/bending and one-point shear
//! integration (selective reduced penalty), not the paper's fully reduced rule
//! and not a claim of equivalence to Abaqus SAX1. Thin-shell measure 2*pi*r ds,
//! laminate ABD and density moments; only meridional/hoop-aligned plies.
//! Pole regularity follows NASA CR-2568 (1975), section 4.1.11: ur=beta=0,
//! interior Gauss integration, and L'Hopital limits for endpoint recovery.
//! This is not the report's Hermite element. A pole is a meridian endpoint,
//! never a continuation into negative radius.
use crate::dense::Matrix;
use crate::frame::tensor_to_global;
use crate::laminate::{Laminate, PlyPoint, PlyProperties};
use crate::shell::ShellStress;
use crate::{Error, Result};

#[derive(Debug, Clone)]
pub struct Axishell {
    x: [[f64; 3]; 2],
    length: f64,
    tangent: [f64; 2],
    section: Laminate,
}
impl Axishell {
    pub fn new(x: &[[f64; 3]], plies: &[PlyProperties], offset: f64) -> Result<Self> {
        let x: [[f64; 3]; 2] = x.try_into().map_err(|_| Error::Invalid("axisymmetric shell requires two ring nodes".into()))?;
        if x.iter().flatten().any(|v| !v.is_finite()) || x.iter().any(|p| p[0] < 0.0 || p[2] != 0.0) {
            return Err(Error::NotSupported(
                "axisymmetric shell requires nonnegative radius and XY meridian; negative-radius axis crossing is not supported".into(),
            ));
        }
        let length = ((x[1][0] - x[0][0]).powi(2) + (x[1][1] - x[0][1]).powi(2)).sqrt();
        if !length.is_finite() || length <= 0.0 {
            return Err(Error::Degenerate("axisymmetric shell zero length".into()));
        }
        for p in plies {
            if p.angle.sin().abs() > 1e-10 && p.angle.cos().abs() > 1e-10 {
                return Err(Error::NotSupported("axisymmetric shell off-axis ply: torsional displacement is absent; use 0/90 degree plies".into()));
            }
        }
        let section = Laminate::new(plies, offset)?;
        let tangent = [(x[1][0] - x[0][0]) / length, (x[1][1] - x[0][1]) / length];
        let radius = 0.5 * (x[0][0] + x[1][0]);
        if radius <= 0.0 {
            return Err(Error::Degenerate("axisymmetric shell meridian lies on the axis".into()));
        }
        // Hoop curvature is |dz/ds|/r, not 1/r. At a smooth spherical pole
        // both numerator and radius vanish. The facet midpoint ratio has a
        // finite limit, whereas h/min(nodal radius) wrongly excludes caps.
        if section.thickness() * tangent[1].abs() > 0.05 * radius {
            return Err(Error::NotSupported("axisymmetric thin-shell thickness times hoop curvature above 0.05".into()));
        }
        Ok(Self { x, length, tangent, section })
    }
    fn shape(&self, q: f64) -> ([f64; 2], [f64; 2], f64) {
        let n = [(1.0 - q) / 2.0, (1.0 + q) / 2.0];
        (n, [-1.0 / self.length, 1.0 / self.length], n[0] * self.x[0][0] + n[1] * self.x[1][0])
    }
    /// Rows [e_s,e_h,0,k_s,k_h,0,gamma,0]; columns ring-by-ring [ur,uz,Rz].
    pub fn operator(&self, q: f64) -> Result<[[f64; 6]; 8]> {
        if !q.is_finite() || q.abs() > 1.0 {
            return Err(Error::Invalid("axisymmetric natural coordinate outside [-1,1]".into()));
        }
        let (n, d, r) = self.shape(q);
        let [c, s] = self.tangent;
        let mut b = [[0.0; 6]; 8];
        for i in 0..2 {
            let j = 3 * i;
            b[0][j] = c * d[i];
            b[0][j + 1] = s * d[i];
            // The constrained pole shape has no radial/rotation contribution.
            // For a regular field ur=O(r), beta=O(r), the remaining shape/r
            // tends to N'/r' at r=0; no artificial epsilon radius is used.
            let hoop = if self.x[i][0] == 0.0 {
                0.0
            } else if r == 0.0 {
                d[i] / c
            } else {
                n[i] / r
            };
            b[1][j] = hoop;
            b[3][j + 2] = -d[i];
            b[4][j + 2] = -c * hoop;
            b[6][j] = -s * d[i];
            b[6][j + 1] = c * d[i];
            b[6][j + 2] = -0.5;
        }
        Ok(b)
    }
    pub fn stiffness(&self) -> Result<Matrix> {
        let mut k = Matrix::zeros(6);
        let abd = self.section.abd();
        for q in [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()] {
            let b = self.operator(q)?;
            let w = std::f64::consts::PI * self.shape(q).2 * self.length;
            for i in 0..6 {
                for j in 0..6 {
                    for a in 0..6 {
                        for z in 0..6 {
                            k.add(i, j, w * b[a][i] * abd[a][z] * b[z][j]);
                        }
                    }
                }
            }
        }
        let b = self.operator(0.0)?;
        let w = 2.0 * std::f64::consts::PI * self.shape(0.0).2 * self.length * self.section.shear()[0][0];
        for i in 0..6 {
            for j in 0..6 {
                k.add(i, j, w * b[6][i] * b[6][j]);
            }
        }
        self.finite(k)
    }
    fn finite(&self, k: Matrix) -> Result<Matrix> {
        if k.data.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("axisymmetric shell matrix overflow".into()));
        }
        Ok(k)
    }
    pub fn mass(&self, lumped: bool) -> Result<Matrix> {
        let [a, b, d] = self.section.mass_moments()?;
        let [c, s] = self.tangent;
        let block = [[a, 0.0, -b * c], [0.0, a, -b * s], [-b * c, -b * s, d]];
        let mut m = Matrix::zeros(6);
        for q in [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()] {
            let (n, _, r) = self.shape(q);
            let w = std::f64::consts::PI * r * self.length;
            for i in 0..2 {
                for j in 0..2 {
                    for u in 0..3 {
                        for v in 0..3 {
                            m.add(3 * i + u, 3 * j + v, w * n[i] * n[j] * block[u][v]);
                        }
                    }
                }
            }
        }
        if lumped {
            // Positive diagonal scaling preserves both ring mass and rotary inertia.
            let total = 2.0 * std::f64::consts::PI * self.shape(0.0).2 * self.length;
            let old = m.clone();
            m = Matrix::zeros(6);
            for u in 0..3 {
                let sum = old.data[u * 6 + u] + old.data[(3 + u) * 6 + 3 + u];
                for i in 0..2 {
                    let j = 3 * i + u;
                    m.add(j, j, old.data[j * 6 + j] * total * block[u][u] / sum);
                }
            }
        }
        self.finite(m)
    }
    pub fn body_force(&self, b: [f64; 3], gravity: bool) -> Result<Vec<f64>> {
        if b.iter().any(|v| !v.is_finite()) || b[2] != 0.0 {
            return Err(Error::NotSupported("axisymmetric shell circumferential/nonfinite body force".into()));
        }
        let moments = if gravity { self.section.mass_moments()? } else { self.section.geometric_moments() };
        let mut f = vec![0.0; 6];
        let [c, s] = self.tangent;
        for q in [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()] {
            let (n, _, r) = self.shape(q);
            let w = std::f64::consts::PI * r * self.length;
            for i in 0..2 {
                f[3 * i] += w * n[i] * moments[0] * b[0];
                f[3 * i + 1] += w * n[i] * moments[0] * b[1];
                f[3 * i + 2] -= w * n[i] * moments[1] * (c * b[0] + s * b[1]);
            }
        }
        if f.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("axisymmetric shell body load overflow".into()));
        }
        Ok(f)
    }
    /// Only a flat annular plate admits the von Karman initial-stress operator.
    pub fn plate_normal(&self) -> Result<[f64; 3]> {
        if self.tangent[1].abs() > 1e-10 {
            return Err(Error::NotSupported("curved/cylindrical/conical axisymmetric-shell buckling; only flat annular plates".into()));
        }
        Ok([-self.tangent[1], self.tangent[0], 0.0])
    }
    /// Delta e_s = (dw/ds)^2/2, Delta e_hoop=0 for axisymmetric transverse w.
    /// Thin flat-plate initial stress: integral 2*pi*r*N_s*w_i'*w_j' ds.
    /// Same restricted flat-preload theory as shell::plate_geometric_stiffness.
    pub fn plate_geometric_stiffness(&self, u: &[f64]) -> Result<Matrix> {
        let normal = self.plate_normal()?;
        let u: &[f64; 6] = u.try_into().map_err(|_| Error::Invalid("axisymmetric preload coordinates".into()))?;
        if u.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("axisymmetric preload nonfinite".into()));
        }
        let scale = u.iter().fold(0.0_f64, |a, b| a.max(b.abs())).max(self.length * f64::EPSILON);
        if (u[1] - u[4]).abs() > 1e-8 * scale || u[2].abs() * self.length > 1e-8 * scale || u[5].abs() * self.length > 1e-8 * scale {
            return Err(Error::NotSupported("axisymmetric annular buckling requires flat membrane preload; coupled pre-bending refused".into()));
        }
        let mut k = Matrix::zeros(6);
        for q in [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()] {
            let b = self.operator(q)?;
            let g = std::array::from_fn(|i| b[i].iter().zip(u).map(|(a, v)| a * v).sum());
            let ns = self.section.resultants(g)?[0];
            let (_, d, r) = self.shape(q);
            let w = std::f64::consts::PI * r * self.length * ns;
            for i in 0..2 {
                for j in 0..2 {
                    k.add(3 * i + 1, 3 * j + 1, w * d[i] * d[j] * normal[1] * normal[1]);
                }
            }
        }
        self.finite(k)
    }
    pub fn stresses(&self, u: &[f64]) -> Result<ShellStress> {
        let u: &[f64; 6] = u.try_into().map_err(|_| Error::Invalid("axisymmetric shell requires six element coordinates".into()))?;
        if u.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("axisymmetric displacement not finite".into()));
        }
        for (i, x) in self.x.iter().enumerate() {
            if x[0] == 0.0 && (u[3 * i] != 0.0 || u[3 * i + 2] != 0.0) {
                return Err(Error::Invalid("axisymmetric shell pole recovery requires ur=beta=0".into()));
            }
        }
        let [c, s] = self.tangent;
        let axes = [[c, s, 0.0], [0.0, 0.0, 1.0], [s, -c, 0.0]];
        // Paper's positive thickness normal is (-s,c), so choose the hoop -Z
        // to retain a right-handed material frame.
        let axes = [axes[0], [0.0, 0.0, -1.0], [-s, c, 0.0]];
        let mut out = ShellStress { axes, forces: vec![], bottom: vec![], middle: vec![], top: vec![], plies: vec![] };
        let rotate = |p: &PlyPoint, angle: f64| {
            let (s, c) = angle.sin_cos();
            let a = [std::array::from_fn(|j| c * axes[0][j] + s * axes[1][j]), std::array::from_fn(|j| -s * axes[0][j] + c * axes[1][j]), axes[2]];
            let [xx, yy, zz, xy, xz, yz] = p.stress;
            tensor_to_global(&a, [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        };
        for q in [-1.0, 1.0] {
            let b = self.operator(q)?;
            let g = std::array::from_fn(|i| b[i].iter().zip(u).map(|(a, v)| a * v).sum());
            let plies = self.section.recover(g)?;
            let first = plies.first().ok_or_else(|| Error::Invalid("axisymmetric shell has no ply".into()))?;
            let last = plies.last().ok_or_else(|| Error::Invalid("axisymmetric shell has no ply".into()))?;
            out.bottom.push(rotate(&first.bottom, first.angle));
            out.top.push(rotate(&last.top, last.angle));
            let (angle, p) = self.section.at(g, self.section.eccentricity())?;
            out.middle.push(rotate(&p, angle));
            out.forces.push(self.section.resultants(g)?);
            out.plies.push(plies);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::Elastic;
    fn shell(x: &[[f64; 3]], h: f64) -> Axishell {
        Axishell::new(
            x,
            &[PlyProperties {
                elastic: Elastic::Isotropic { young: 210e9, poisson: 0.3 },
                thickness: h,
                angle: 0.0,
                points: 3,
                density: Some(7800.0),
                strength: None,
            }],
            0.0,
        )
        .unwrap()
    }
    fn energy(k: &Matrix, u: &[f64; 6]) -> f64 {
        (0..6).map(|i| (0..6).map(|j| u[i] * k.data[6 * i + j] * u[j]).sum::<f64>()).sum()
    }
    #[test]
    fn axial_rigid_and_ring_membrane_have_independent_energies() {
        let s = shell(&[[1.0, 0.0, 0.0], [1.0, 2.0, 0.0]], 0.01);
        let rigid = [0.0, 0.123, 0.0, 0.0, 0.123, 0.0];
        for q in [-1.0, 0.0, 1.0] {
            for b in s.operator(q).unwrap() {
                assert!(b.iter().zip(rigid).map(|(a, u)| a * u).sum::<f64>().abs() < 1e-14);
            }
        }
        let q = 210e9 / (1.0 - 0.3 * 0.3);
        let eps = 0.001;
        let u = [eps, 0.0, 0.0, eps, 2.0 * eps, 0.0];
        let expected = 2.0 * std::f64::consts::PI * 2.0 * 0.01 * q * 2.0 * 1.3 * eps * eps;
        assert!((energy(&s.stiffness().unwrap(), &u) / expected - 1.0).abs() < 1e-12);
        for lump in [false, true] {
            let m = s.mass(lump).unwrap();
            let exact = 2.0 * std::f64::consts::PI * 2.0 * 0.01 * 7800.0 * 0.123_f64.powi(2);
            assert!((energy(&m, &rigid) / exact - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn annular_constant_curvature_patch_and_thin_limit() {
        for h in [0.01, 0.001, 0.0001] {
            let s = shell(&[[1.0, 0.0, 0.0], [1.2, 0.0, 0.0]], h);
            let a = 0.002;
            let u = [0.0, a / 2.0, a, 0.0, a * 1.2 * 1.2 / 2.0, a * 1.2];
            let exact = std::f64::consts::PI * (1.2 * 1.2 - 1.0) * 210e9 * h.powi(3) / (12.0 * (1.0 - 0.3 * 0.3)) * 2.0 * 1.3 * a * a;
            assert!((energy(&s.stiffness().unwrap(), &u) / exact - 1.0).abs() < 2e-7, "h={h}");
        }
    }
    #[test]
    fn conical_membrane_and_reversed_connectivity() {
        let x = [[1.0, 0.0, 0.0], [1.2, 0.3, 0.0]];
        let a = shell(&x, 0.01);
        let b = shell(&[x[1], x[0]], 0.01);
        let k = a.stiffness().unwrap();
        let reverse = b.stiffness().unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let ri = (i + 3) % 6;
                let rj = (j + 3) % 6;
                assert!((k.data[6 * i + j] - reverse.data[6 * ri + rj]).abs() < 1e-10 * k.data.iter().fold(0.0_f64, |s, x| s.max(x.abs())));
            }
        }
        let e = 0.001;
        let u = [e, 0.0, 0.0, 1.2 * e, 0.3 * e, 0.0];
        let q = 210e9 / (1.0 - 0.3 * 0.3);
        let reference = 2.0 * std::f64::consts::PI * 1.1 * a.length * 0.01 * q * 2.0 * 1.3 * e * e;
        assert!((energy(&k, &u) / reference - 1.0).abs() < 1e-12);
        // Global axial translation is shared at both ends even when adjacent facets turn.
        let rigid = [0.0, 0.1, 0.0, 0.0, 0.1, 0.0];
        for q in [-1.0, 0.0, 1.0] {
            for row in a.operator(q).unwrap() {
                assert!(row.iter().zip(rigid).map(|(b, u)| b * u).sum::<f64>().abs() < 1e-14);
            }
        }
    }
    #[test]
    fn regular_pole_limits_patch_mass_and_connectivity() {
        for reverse in [false, true] {
            let mut x = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
            let mut u = [0.0, 0.0, 0.0, 0.001, 0.001, 0.002];
            if reverse {
                x.reverse();
                u.rotate_left(3);
            }
            let s = shell(&x, 0.01);
            for q in [-1.0, -0.999999, 0.0, 0.999999, 1.0] {
                let b = s.operator(q).unwrap();
                let e: [f64; 8] = std::array::from_fn(|i| b[i].iter().zip(u).map(|(b, u)| b * u).sum());
                assert!((e[0] - 0.001).abs() < 1e-14);
                assert!((e[1] - 0.001).abs() < 1e-14);
                let sign = if reverse { 1.0 } else { -1.0 };
                assert!((e[3] - sign * 0.002).abs() < 1e-14);
                assert!((e[4] - sign * 0.002).abs() < 1e-14);
                assert!(e[6].abs() < 1e-14);
            }
            let rigid = [0.0, 1.0, 0.0, 0.0, 1.0, 0.0];
            for row in s.operator(if reverse { 1.0 } else { -1.0 }).unwrap() {
                assert!(row.iter().zip(rigid).map(|(b, u)| b * u).sum::<f64>().abs() < 1e-14);
            }
            for lumped in [false, true] {
                assert!((energy(&s.mass(lumped).unwrap(), &rigid) / (std::f64::consts::PI * 0.01 * 7800.0) - 1.0).abs() < 1e-12);
            }
            assert!(s.stiffness().unwrap().data.iter().all(|x| x.is_finite()));
        }
        // A smooth spherical cap remains admissible as its first facet shrinks.
        for theta in [0.1_f64, 0.01, 0.001] {
            let s = shell(&[[0.0, 1.0, 0.0], [theta.sin(), theta.cos(), 0.0]], 0.025);
            assert!(s.operator(-1.0).unwrap().iter().flatten().all(|v| v.is_finite()));
        }
    }
    #[test]
    fn pole_off_axis_layup_and_thick_shell_are_refused() {
        let p = PlyProperties {
            elastic: Elastic::Isotropic { young: 1e6, poisson: 0.3 },
            thickness: 0.01,
            angle: 0.3,
            points: 3,
            density: Some(1.0),
            strength: None,
        };
        assert!(Axishell::new(&[[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], &[p], 0.0).is_err());
        assert!(Axishell::new(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &[p], 0.0).is_err());
    }
}
