//! Local frames of structural elements (shells, beams, rods) and the rotation of their matrices
//! between the element's local axes and the global ones.
//!
//! A frame is a rotation R whose rows are the local axes e₁, e₂, e₃ in global components, so
//! local components are a_l = R a_g. Element degrees of freedom come in groups of three (the
//! translations, then the rotations, of each node), each rotated by R, so the element matrix is
//! K_g = Tᵀ K_l T with T = diag(R, R, …) (Cook, Malkus, Plesha & Witt 2002, §2.10 and §7.4).

use crate::dense::Matrix;

pub type Vec3 = [f64; 3];

/// The rows e₁, e₂, e₃ of a rotation from global to local components.
pub type Axes = [Vec3; 3];

pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// s a.
pub fn scale(a: Vec3, s: f64) -> Vec3 {
    a.map(|c| c * s)
}

pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// The unit vector along `a`, or `None` when it is zero or not finite.
pub fn unit(a: Vec3) -> Option<Vec3> {
    let n = norm(a);
    (n > 0.0 && n.is_finite()).then(|| a.map(|c| c / n))
}

/// R a: global components to local.
pub fn to_local(r: &Axes, a: Vec3) -> Vec3 {
    [dot(r[0], a), dot(r[1], a), dot(r[2], a)]
}

/// Rᵀ a: local components to global.
pub fn to_global(r: &Axes, a: Vec3) -> Vec3 {
    let mut g = [0.0; 3];
    for (k, axis) in r.iter().enumerate() {
        for c in 0..3 {
            g[c] += axis[c] * a[k];
        }
    }
    g
}

/// Tᵀ K T for a symmetric matrix whose degrees of freedom come in groups of three local
/// components; the result is exactly symmetric.
pub fn rotate_matrix(k: &Matrix, r: &Axes) -> Matrix {
    let n = k.n;
    let groups = n / 3;
    let mut out = Matrix::zeros(n);
    for a in 0..groups {
        for b in a..groups {
            // Rᵀ K_ab R, mirrored to (b, a): K is symmetric.
            let mut kr = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    kr[i][j] = (0..3).map(|m| k.get(3 * a + i, 3 * b + m) * r[m][j]).sum();
                }
            }
            for i in 0..3 {
                for j in 0..3 {
                    let v: f64 = (0..3).map(|m| r[m][i] * kr[m][j]).sum();
                    out.add(3 * a + i, 3 * b + j, v);
                    if a != b {
                        out.add(3 * b + j, 3 * a + i, v);
                    }
                }
            }
        }
    }
    crate::element::symmetrize(&mut out);
    out
}

/// T v: a vector of global components in groups of three to local components.
pub fn rotate_to_local(v: &[f64], r: &Axes) -> Vec<f64> {
    v.as_chunks::<3>().0.iter().flat_map(|g| to_local(r, *g)).collect()
}

/// Tᵀ v: local components in groups of three to global ones.
pub fn rotate_to_global(v: &[f64], r: &Axes) -> Vec<f64> {
    v.as_chunks::<3>().0.iter().flat_map(|g| to_global(r, *g)).collect()
}

/// A uniform translation `b` of `nodes` nodes with six degrees of freedom each (rotations
/// zero): the nodal field whose product with a consistent mass gives body-force loads.
pub fn translation(b: Vec3, nodes: usize) -> Vec<f64> {
    (0..nodes).flat_map(|_| [b[0], b[1], b[2], 0.0, 0.0, 0.0]).collect()
}

/// A symmetric tensor in local axes, [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], in global
/// components as the 6-vector xx, yy, zz, xy, xz, yz: Rᵀ S R.
pub fn tensor_to_global(r: &Axes, s: [[f64; 3]; 3]) -> [f64; 6] {
    let mut g = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            g[i][j] = (0..3).map(|a| (0..3).map(|b| r[a][i] * s[a][b] * r[b][j]).sum::<f64>()).sum();
        }
    }
    [g[0][0], g[1][1], g[2][2], g[0][1], g[0][2], g[1][2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotations_round_trip() {
        let e1 = unit([1.0, 1.0, 0.0]).unwrap();
        let e3 = [0.0, 0.0, 1.0];
        let r: Axes = [e1, cross(e3, e1), e3];
        let v = [0.3, -1.2, 2.0];
        let back = to_global(&r, to_local(&r, v));
        assert!(back.iter().zip(v).all(|(a, b)| (a - b).abs() < 1e-15));
        // A uniaxial stress along e1 is xx = yy = xy = σ/2 globally.
        let s = tensor_to_global(&r, [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        assert!((s[0] - 1.0).abs() < 1e-15 && (s[1] - 1.0).abs() < 1e-15 && (s[3] - 1.0).abs() < 1e-15 && s[2] == 0.0, "{s:?}");
        // Rotating a matrix keeps its quadratic form: (T v)ᵀ K (T v) = vᵀ (Tᵀ K T) v.
        let mut k = Matrix::zeros(6);
        for i in 0..6 {
            for j in 0..6 {
                k.data[i * 6 + j] = 1.0 / (1.0 + i as f64 + j as f64);
            }
        }
        let g = rotate_matrix(&k, &r);
        let x = [0.1, 0.2, -0.3, 0.7, 0.0, 1.0];
        let xl = rotate_to_local(&x, &r);
        let ql: f64 = (0..6).map(|i| (0..6).map(|j| xl[i] * k.get(i, j) * xl[j]).sum::<f64>()).sum();
        let qg: f64 = (0..6).map(|i| (0..6).map(|j| x[i] * g.get(i, j) * x[j]).sum::<f64>()).sum();
        assert!((ql - qg).abs() < 1e-14, "{ql} {qg}");
        assert!(unit([0.0; 3]).is_none() && unit([f64::NAN, 0.0, 0.0]).is_none());
    }
}
