//! Linear elastic constitutive matrices.
//!
//! Stress and strain are 6-vectors in the order xx, yy, zz, xy, xz, yz, with engineering shear
//! strains (γ = 2ε), so σ = D ε. The order puts the orthotropic shear moduli G₁₂, G₁₃, G₂₃ on
//! the diagonal in the order the model lists them. Material axes are the global axes until the
//! model has material orientations.
//!
//! Orthotropic stiffness is the inverse of the compliance written with engineering constants
//! (R. M. Jones, *Mechanics of Composite Materials*, 2nd ed., Taylor & Francis 1999, §2.4):
//! S₁₁ = 1/E₁, S₁₂ = −ν₁₂/E₁, S₁₃ = −ν₁₃/E₁, S₂₃ = −ν₂₃/E₂ (ν_ji/E_j = ν_ij/E_i), and
//! 1/G on the shear diagonal.

use fieldcraft_model::Elastic;

use crate::dense::inverse3;
use crate::{Error, Result};

/// A 6 × 6 elasticity matrix, row-major.
pub type Elasticity = [[f64; 6]; 6];

/// Check isotropic constants: E > 0 and −1 < ν < 0.5.
fn check_isotropic(young: f64, poisson: f64) -> Result<()> {
    if !(young > 0.0 && young.is_finite() && poisson > -1.0 && poisson < 0.5) {
        return Err(Error::Invalid(format!("isotropic constants E = {young}, ν = {poisson} are not a stable material")));
    }
    Ok(())
}

/// The checked (E, ν) of an isotropic material, for the structural sections (shells, beams,
/// rods) that need one; `what` names the section kind in the error for other materials.
pub fn isotropic_constants(elastic: &Elastic, what: &str) -> Result<(f64, f64)> {
    let Elastic::Isotropic { young, poisson } = *elastic else {
        return Err(Error::NotSupported(format!("orthotropic {what} sections (they need material orientations)")));
    };
    check_isotropic(young, poisson)?;
    Ok((young, poisson))
}

/// The isotropic elasticity matrix.
pub fn isotropic(young: f64, poisson: f64) -> Result<Elasticity> {
    check_isotropic(young, poisson)?;
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    let mut d = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            d[i][j] = lambda;
        }
        d[i][i] = lambda + 2.0 * mu;
        d[i + 3][i + 3] = mu;
    }
    Ok(d)
}

/// Engineering constants whose compliance's determinant is below this fraction of the product of
/// its diagonal are not positive definite (or as good as singular).
const SINGULAR_COMPLIANCE: f64 = 1e-12; // unit-factors: relative to the diagonal's product

/// The orthotropic elasticity matrix from engineering constants.
pub fn orthotropic(young: [f64; 3], poisson: [f64; 3], shear: [f64; 3]) -> Result<Elasticity> {
    let [e1, e2, e3] = young;
    let [n12, n13, n23] = poisson;
    let s = [[1.0 / e1, -n12 / e1, -n13 / e1], [-n12 / e1, 1.0 / e2, -n23 / e2], [-n13 / e1, -n23 / e2, 1.0 / e3]];
    let scale = (1.0 / e1) * (1.0 / e2) * (1.0 / e3);
    let inverse = inverse3(&s).filter(|(_, det)| *det > SINGULAR_COMPLIANCE * scale);
    let Some((normal, _)) = inverse.filter(|_| shear.iter().all(|g| *g > 0.0 && g.is_finite())) else {
        return Err(Error::Invalid("the orthotropic constants do not describe a stable material".into()));
    };
    let mut d = [[0.0; 6]; 6];
    for i in 0..3 {
        d[i][..3].copy_from_slice(&normal[i]);
    }
    for (k, g) in shear.iter().enumerate() {
        d[3 + k][3 + k] = *g;
    }
    Ok(d)
}

/// The elasticity matrix of a model material's elastic behaviour.
pub fn elasticity(e: &Elastic) -> Result<Elasticity> {
    match *e {
        Elastic::Isotropic { young, poisson } => isotropic(young, poisson),
        Elastic::EngineeringConstants { young, poisson, shear } => orthotropic(young, poisson, shear),
    }
}

/// The von Mises equivalent stress of a stress 6-vector.
pub fn von_mises(s: &[f64; 6]) -> f64 {
    let [xx, yy, zz, xy, xz, yz] = *s;
    let normal = (xx - yy) * (xx - yy) + (yy - zz) * (yy - zz) + (zz - xx) * (zz - xx);
    (0.5 * normal + 3.0 * (xy * xy + xz * xz + yz * yz)).sqrt()
}

/// The principal stresses of a stress 6-vector, largest first: the eigenvalues of the symmetric
/// stress tensor in closed form (O. K. Smith, "Eigenvalues of a symmetric 3 × 3 matrix",
/// Commun. ACM 4(4):168, 1961).
pub fn principal(s: &[f64; 6]) -> [f64; 3] {
    let [xx, yy, zz, xy, xz, yz] = *s;
    let off = xy * xy + xz * xz + yz * yz;
    let q = (xx + yy + zz) / 3.0;
    let (a, b, c) = (xx - q, yy - q, zz - q);
    let p = ((a * a + b * b + c * c + 2.0 * off) / 6.0).sqrt();
    if !(p > 0.0) || off == 0.0 {
        let mut d = [xx, yy, zz];
        d.sort_by(|u, v| v.total_cmp(u));
        return d;
    }
    // det(A − qI) / (2 p³), the cosine of three times the angle, within [−1, 1].
    let det = a * (b * c - yz * yz) - xy * (xy * c - yz * xz) + xz * (xy * yz - b * xz);
    let r = (det / (2.0 * p * p * p)).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let first = q + 2.0 * p * phi.cos();
    let third = q + 2.0 * p * (phi + 2.0 * std::f64::consts::FRAC_PI_3).cos();
    [first, 3.0 * q - first - third, third]
}

/// σ = D ε.
pub fn stress(d: &Elasticity, strain: &[f64; 6]) -> [f64; 6] {
    let mut s = [0.0; 6];
    for (si, row) in s.iter_mut().zip(d) {
        *si = row.iter().zip(strain).map(|(a, b)| a * b).sum();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orthotropic_with_isotropic_constants_is_isotropic() {
        let (e, nu) = (210e3, 0.3);
        let g = e / (2.0 * (1.0 + nu));
        let a = isotropic(e, nu).unwrap();
        let b = orthotropic([e; 3], [nu; 3], [g; 3]).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                assert!((a[i][j] - b[i][j]).abs() < 1e-9 * e, "{i} {j}: {} vs {}", a[i][j], b[i][j]);
            }
        }
        assert!(isotropic(1.0, 0.5).is_err());
        assert!(orthotropic([1.0; 3], [0.9; 3], [1.0; 3]).is_err());
    }

    #[test]
    fn principal_stresses_are_the_tensor_eigenvalues() {
        let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9 * (1.0 + y.abs()));
        // Diagonal: sorted, largest first.
        assert_eq!(principal(&[1.0, 3.0, -2.0, 0.0, 0.0, 0.0]), [3.0, 1.0, -2.0]);
        // Pure shear τ in the xy plane: ±τ and 0.
        assert!(close(principal(&[0.0, 0.0, 0.0, 5.0, 0.0, 0.0]), [5.0, 0.0, -5.0]));
        // Hydrostatic: three equal.
        assert!(close(principal(&[-7.0, -7.0, -7.0, 0.0, 0.0, 0.0]), [-7.0; 3]));
        // A general tensor: the sum is the trace, and each is a root of det(σ − λI) = 0.
        let s = [10.0, -4.0, 3.0, 2.0, -1.5, 6.0];
        let l = principal(&s);
        assert!(l[0] >= l[1] && l[1] >= l[2]);
        assert!((l.iter().sum::<f64>() - (s[0] + s[1] + s[2])).abs() < 1e-9);
        for x in l {
            let (a, b, c) = (s[0] - x, s[1] - x, s[2] - x);
            let det = a * (b * c - s[5] * s[5]) - s[3] * (s[3] * c - s[5] * s[4]) + s[4] * (s[3] * s[5] - b * s[4]);
            assert!(det.abs() < 1e-7, "{x}: {det}");
        }
        // Von Mises from the principal values agrees.
        let vm = (0.5 * ((l[0] - l[1]).powi(2) + (l[1] - l[2]).powi(2) + (l[2] - l[0]).powi(2))).sqrt();
        assert!((vm - von_mises(&s)).abs() < 1e-9);
        assert!(principal(&[f64::NAN; 6]).iter().all(|v| v.is_nan()));
    }

    #[test]
    fn orthotropic_uniaxial_response() {
        // Stress along 1 alone gives strains σ/E1, −ν12 σ/E1, −ν13 σ/E1.
        let (young, poisson) = ([140e3, 10e3, 12e3], [0.3, 0.25, 0.4]);
        let d = orthotropic(young, poisson, [5e3, 4e3, 3.8e3]).unwrap();
        let strain = [1.0 / young[0], -poisson[0] / young[0], -poisson[1] / young[0], 0.0, 0.0, 0.0];
        let s = stress(&d, &strain);
        assert!((s[0] - 1.0).abs() < 1e-12 && s[1].abs() < 1e-12 && s[2].abs() < 1e-12, "{s:?}");
        assert!((von_mises(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]) - 1.0).abs() < 1e-15);
        assert!((von_mises(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]) - 3f64.sqrt()).abs() < 1e-15);
    }
}
