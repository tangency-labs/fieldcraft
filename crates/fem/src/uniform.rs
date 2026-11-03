//! Uniform-strain (one-point) kinematics of the eight-node hexahedron and the four-node
//! tetrahedron, for explicit dynamics.
//!
//! **Hexahedron** (D. P. Flanagan, T. Belytschko, *A uniform strain hexahedron and quadrilateral
//! with orthogonal hourglass control*, IJNME 17(5):679–706, 1981): the strain is the element's
//! volume average, computed with the integrated gradients b_iI = ∫ ∂N_I/∂x_i dV, so
//! ε̄ = (1/V) Σ_I sym(u_I ⊗ b_I) and the nodal forces of a stress are f_iI = σ_ij b_jI. The
//! mean strain is exact for every linear displacement field on any hexahedron (the element
//! passes the patch test). The products det J · ∂N/∂x are polynomials of degree two at most in
//! each natural coordinate, so the 2 × 2 × 2 Gauss rule integrates b and V exactly.
//!
//! The one-point element has twelve hourglass modes that the mean strain does not see: in each
//! direction, the four base vectors h₁ = ηζ, h₂ = ξζ, h₃ = ξη, h₄ = ξηζ at the corners. The
//! hourglass shape vectors γ_α = h_α − (1/V) Σ_J (h_αJ x_iJ) b_i (Flanagan & Belytschko 1981)
//! are orthogonal to every linear velocity field, so forces resisting γ_α leave the constant
//! strain modes alone and the patch test passes on distorted elements too.
//!
//! **Tetrahedron**: the linear tetrahedron's strain is uniform, b_I = V ∇N_I.
//!
//! **Critical frequency.** The element eigenvalue inequality (T. Belytschko, W. K. Liu, B. Moran,
//! K. Elkhodary, *Nonlinear Finite Elements for Continua and Structures*, 2nd ed., Wiley 2014,
//! §6.3) bounds the mesh's highest frequency by the largest element frequency. For the uniform
//! strain stiffness K = (1/V) B̃ᵀ D B̃ we bound it in closed form, in the spirit of D. P.
//! Flanagan, T. Belytschko, *Eigenvalues and stable time steps for the uniform strain hexahedron
//! and quadrilateral*, J. Appl. Mech. 51(1):35–40, 1984: with G = Σ_I u_I ⊗ b_I,
//! uᵀ K u ≤ (k/V) ‖G‖² ≤ (k/V) λ_max(Σ_I b_I b_Iᵀ) ‖u‖², where k bounds the material's
//! energy per squared displacement gradient ([`modulus_bound`]); with the lumped masses m_I the
//! element frequency satisfies ω² ≤ k λ_max(Σ b bᵀ) / (V min m_I). The bound is exact for a
//! cube (its dilatation mode) and conservative for distorted elements.

use fieldcraft_model::Elastic;

use crate::dense::positive_inverse;
use crate::material::Elasticity;
use crate::quadrature::Rule;
use crate::shape::{self, HEX_CORNERS};
use crate::{Error, Result};

/// The integrated gradients, volume and nodal volumes of a uniform-strain element with `N` nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformStrain<const N: usize> {
    /// b_I = ∫ ∇N_I dV.
    pub b: [[f64; 3]; N],
    pub volume: f64,
    /// ∫ N_I dV: the row sums of the consistent mass per unit density (row-sum lumping).
    pub nodal_volume: [f64; N],
}

/// The one-point hexahedron.
pub type Hex8 = UniformStrain<8>;
/// The linear tetrahedron.
pub type Tet4 = UniformStrain<4>;

fn degenerate(what: &str) -> Error {
    Error::Degenerate(format!("a {what} element is inverted or degenerate (its Jacobian is not positive)"))
}

impl UniformStrain<8> {
    /// The uniform-strain hexahedron of corners `xyz` (the model's hex8 node order).
    pub fn hex8(xyz: &[[f64; 3]; 8]) -> Result<Hex8> {
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        let mut b = [[0.0; 3]; 8];
        let mut nodal_volume = [0.0; 8];
        let mut volume = 0.0;
        for p in Rule::Hex(2).points() {
            let s = shape::hex8(p.xi);
            let mut j = [[0.0; 3]; 3];
            for (d, x) in s.derivatives().iter().zip(xyz) {
                for a in 0..3 {
                    for c in 0..3 {
                        j[a][c] += d[a] * x[c];
                    }
                }
            }
            let (inv, det) = positive_inverse(&j).ok_or_else(|| degenerate("hex8"))?;
            let w = det * p.weight;
            volume += w;
            for (((bi, vi), d), n) in b.iter_mut().zip(nodal_volume.iter_mut()).zip(s.derivatives()).zip(s.values()) {
                for c in 0..3 {
                    bi[c] += w * (inv[c][0] * d[0] + inv[c][1] * d[1] + inv[c][2] * d[2]);
                }
                *vi += w * n;
            }
        }
        Ok(UniformStrain { b, volume, nodal_volume })
    }

    /// The hourglass shape vectors γ_α (α = 1..4) of the element with corners `xyz`.
    pub fn hourglass_vectors(&self, xyz: &[[f64; 3]; 8]) -> [[f64; 8]; 4] {
        let mut gamma = [[0.0; 8]; 4];
        for (alpha, g) in gamma.iter_mut().enumerate() {
            let mut h = [0.0f64; 8];
            for (hk, c) in h.iter_mut().zip(&HEX_CORNERS) {
                *hk = match alpha {
                    0 => c[1] * c[2],
                    1 => c[0] * c[2],
                    2 => c[0] * c[1],
                    _ => c[0] * c[1] * c[2],
                };
            }
            // (h · x_i) for each coordinate direction i, from −0 as `Iterator::sum` starts (loops:
            // the array maps were not inlined in the explicit solver's element loop).
            let mut hx = [-0.0f64; 3];
            for (hj, x) in h.iter().zip(xyz) {
                for i in 0..3 {
                    hx[i] += hj * x[i];
                }
            }
            for (k, gk) in g.iter_mut().enumerate() {
                *gk = h[k] - (hx[0] * self.b[k][0] + hx[1] * self.b[k][1] + hx[2] * self.b[k][2]) / self.volume;
            }
        }
        gamma
    }
}

/// The hexahedron's volume as a form in its corner positions (in the model's hex8 node order):
/// with C_IJK = ∫ ε_abc ∂_aN_I ∂_bN_J ∂_cN_K dξ (the determinant of the Jacobian expanded over
/// the corners), which is antisymmetric in I, J, K, V = Σ_{I<J<K} C_IJK x_I · (x_J × x_K). The 24
/// nonzero coefficients with I < J < K are ±1/12; this table holds them times 12 (the tests
/// integrate them again).
const VOLUME_FORM: [(usize, usize, usize, f64); 24] = [
    (0, 1, 2, -1.0),
    (0, 1, 3, -1.0),
    (0, 1, 4, 1.0),
    (0, 1, 5, 1.0),
    (0, 2, 3, -1.0),
    (0, 3, 4, -1.0),
    (0, 3, 7, -1.0),
    (0, 4, 5, -1.0),
    (0, 4, 7, 1.0),
    (1, 2, 3, -1.0),
    (1, 2, 5, 1.0),
    (1, 2, 6, 1.0),
    (1, 4, 5, -1.0),
    (1, 5, 6, -1.0),
    (2, 3, 6, 1.0),
    (2, 3, 7, 1.0),
    (2, 5, 6, -1.0),
    (2, 6, 7, -1.0),
    (3, 4, 7, 1.0),
    (3, 6, 7, -1.0),
    (4, 5, 6, 1.0),
    (4, 5, 7, 1.0),
    (4, 6, 7, 1.0),
    (5, 6, 7, 1.0),
];

/// The index of pair p < q among the 28 pairs of a hexahedron's corners taken in order (0 1, 0 2,
/// …, 6 7), as [`Hex8::hex8_deformed`] lists their cross products.
const fn pair_index(p: usize, q: usize) -> usize {
    // pairs before p: 7 + 6 + … + (8 − p); then q − p − 1 within p's
    p * (15 - p) / 2 + q - p - 1
}

/// 12 b_I of the volume form ([`VOLUME_FORM`]) term by term: for each corner I, its nine terms
/// (a pair's index ([`pair_index`]) and the sign it is added with) in the order the form adds them,
/// term (i, j, k, c) adding c x_j × x_k to b_i, c x_k × x_i to b_j and c x_i × x_j to b_k. x_k × x_i
/// is −(x_i × x_k): the same two products subtracted the other way, equal but for a zero's sign,
/// and a zero of either sign added to a sum that starts at +0 leaves it as it was, so this is the
/// form's sum bit for bit.
const B_TERMS: [[(usize, f64); 9]; 8] = b_terms();

const fn b_terms() -> [[(usize, f64); 9]; 8] {
    let mut out = [[(0, 0.0); 9]; 8];
    let mut count = [0usize; 8];
    let mut m = 0;
    while m < VOLUME_FORM.len() {
        let (i, j, k, c) = VOLUME_FORM[m];
        out[i][count[i]] = (pair_index(j, k), c);
        count[i] += 1;
        out[j][count[j]] = (pair_index(i, k), -c);
        count[j] += 1;
        out[k][count[k]] = (pair_index(i, j), c);
        count[k] += 1;
        m += 1;
    }
    out
}

impl UniformStrain<8> {
    /// The integrated gradients and volume of a hexahedron in closed form: b_I = ∂V/∂x_I of the
    /// volume form ([`VOLUME_FORM`]; D. P. Flanagan, T. Belytschko 1981 give it written out),
    /// so b_I gathers C_IJK x_J × x_K, and V = ⅓ Σ_I b_I · x_I. A few hundred operations instead of a Gauss loop, for the
    /// current configuration of large deformation. The element is refused when its volume is not
    /// positive. `nodal_volume` is left at the mean V/8: an explicit solver keeps the lumped
    /// masses of the undeformed element.
    pub fn hex8_deformed(xyz: &[[f64; 3]; 8]) -> Result<Hex8> {
        let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        // x_p × x_q of each pair p < q, worked out once: the 24 terms of the form take 72 cross
        // products of these 28.
        // Written out in [`pair_index`] order (a loop over a table zeroed the array first).
        let x = |p: usize, q: usize| cross(xyz[p], xyz[q]);
        let pair: [[f64; 3]; 28] = [
            x(0, 1),
            x(0, 2),
            x(0, 3),
            x(0, 4),
            x(0, 5),
            x(0, 6),
            x(0, 7),
            x(1, 2),
            x(1, 3),
            x(1, 4),
            x(1, 5),
            x(1, 6),
            x(1, 7),
            x(2, 3),
            x(2, 4),
            x(2, 5),
            x(2, 6),
            x(2, 7),
            x(3, 4),
            x(3, 5),
            x(3, 6),
            x(3, 7),
            x(4, 5),
            x(4, 6),
            x(4, 7),
            x(5, 6),
            x(5, 7),
            x(6, 7),
        ];
        // Each b_I summed on its own, its terms in the order the form adds them ([`B_TERMS`]),
        // so the sums stay in registers and every b_I is the same bit for bit.
        let mut b = [[0.0f64; 3]; 8];
        for (bi, terms) in b.iter_mut().zip(&B_TERMS) {
            let mut sum = [0.0f64; 3];
            for &(k, s) in terms {
                let v = pair[k];
                for d in 0..3 {
                    sum[d] += s * v[d];
                }
            }
            for d in 0..3 {
                bi[d] = sum[d] / 12.0;
            }
        }
        let volume: f64 = b.iter().zip(xyz).map(|(bi, x)| bi[0] * x[0] + bi[1] * x[1] + bi[2] * x[2]).sum::<f64>() / 3.0;
        if !(volume > 0.0 && volume.is_finite() && b.iter().flatten().all(|v| v.is_finite())) {
            return Err(degenerate("hex8"));
        }
        Ok(UniformStrain { b, volume, nodal_volume: [volume / 8.0; 8] })
    }
}

impl UniformStrain<4> {
    /// The linear tetrahedron of corners `xyz` (the model's tet4 node order).
    pub fn tet4(xyz: &[[f64; 3]; 4]) -> Result<Tet4> {
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        // J[a][c] = ∂x_c/∂ξ_a: the edges from corner 1.
        let j = fieldcraft_model::tet4_jacobian(xyz);
        let (inv, det) = positive_inverse(&j).ok_or_else(|| degenerate("tet4"))?;
        let volume = det / 6.0;
        let s = shape::tet4([0.25; 3]);
        let mut b = [[0.0; 3]; 4];
        for (bi, d) in b.iter_mut().zip(s.derivatives()) {
            for c in 0..3 {
                bi[c] = volume * (inv[c][0] * d[0] + inv[c][1] * d[1] + inv[c][2] * d[2]);
            }
        }
        Ok(UniformStrain { b, volume, nodal_volume: [volume / 4.0; 4] })
    }
}

impl<const N: usize> UniformStrain<N> {
    /// The displacement gradient (1/V) Σ_I u_I ⊗ b_I: `[i][j]` = ∂u_i/∂x_j.
    pub fn gradient(&self, u: &[[f64; 3]; N]) -> [[f64; 3]; 3] {
        let mut g = [[0.0; 3]; 3];
        for (ui, bi) in u.iter().zip(&self.b) {
            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += ui[i] * bi[j];
                }
            }
        }
        g.map(|r| r.map(|v| v / self.volume))
    }

    /// The mean strain of nodal displacements (or its rate, of nodal velocities): xx, yy, zz,
    /// xy, xz, yz with engineering shears.
    pub fn strain(&self, u: &[[f64; 3]; N]) -> [f64; 6] {
        let g = self.gradient(u);
        [g[0][0], g[1][1], g[2][2], g[0][1] + g[1][0], g[0][2] + g[2][0], g[1][2] + g[2][1]]
    }

    /// The trace of the mean gradient (1/V) Σ_I u_I · b_I: the volumetric strain of nodal
    /// displacements, or its rate of nodal velocities.
    pub fn divergence(&self, u: &[[f64; 3]; N]) -> f64 {
        let mut d = [0.0; 3];
        for (ui, bi) in u.iter().zip(&self.b) {
            for i in 0..3 {
                d[i] += ui[i] * bi[i];
            }
        }
        (d[0] + d[1] + d[2]) / self.volume
    }

    /// The nodal forces of a uniform pressure p (positive in compression): f_I = −p b_I.
    pub fn pressure_forces(&self, p: f64) -> [[f64; 3]; N] {
        self.b.map(|bi| bi.map(|c| -p * c))
    }

    /// The nodal forces that balance a uniform stress (xx, yy, zz, xy, xz, yz): f_iI = σ_ij b_jI,
    /// the internal force ∫ Bᵀ σ dV.
    pub fn forces(&self, s: &[f64; 6]) -> [[f64; 3]; N] {
        let sigma = [[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]];
        self.b.map(|bi| [0, 1, 2].map(|i| sigma[i][0] * bi[0] + sigma[i][1] * bi[1] + sigma[i][2] * bi[2]))
    }

    /// The largest eigenvalue of Σ_I b_I b_Iᵀ.
    pub fn gradient_norm(&self) -> f64 {
        let mut m = [[0.0; 3]; 3];
        for bi in &self.b {
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] += bi[i] * bi[j];
                }
            }
        }
        largest_eigenvalue(&m)
    }

    /// An upper bound of the element's highest natural frequency squared (ω²) with the
    /// row-sum lumped mass, for a material whose [`modulus_bound`] is `modulus`.
    pub fn frequency_bound(&self, modulus: f64, density: f64) -> f64 {
        let lightest = self.nodal_volume.iter().copied().fold(f64::INFINITY, f64::min);
        modulus * self.gradient_norm() / (self.volume * density * lightest)
    }
}

/// The largest eigenvalue of a symmetric 3 × 3 matrix, in closed form (O. K. Smith, *Eigenvalues
/// of a symmetric 3 × 3 matrix*, Commun. ACM 4(4):168, 1961).
pub fn largest_eigenvalue(a: &[[f64; 3]; 3]) -> f64 {
    let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
    let q = (a[0][0] + a[1][1] + a[2][2]) / 3.0;
    let p2 = (a[0][0] - q).powi(2) + (a[1][1] - q).powi(2) + (a[2][2] - q).powi(2) + 2.0 * off;
    if !(p2 > 1e-30 * q * q) {
        return a[0][0].max(a[1][1]).max(a[2][2]);
    }
    let p = (p2 / 6.0).sqrt();
    let mut b = *a;
    for (i, row) in b.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = (*v - if i == j { q } else { 0.0 }) / p;
        }
    }
    let r = (crate::dense::det3(&b) / 2.0).clamp(-1.0, 1.0);
    q + 2.0 * p * (r.acos() / 3.0).cos()
}

/// A bound k of a material's strain energy density per squared displacement gradient:
/// εᵀ D ε ≤ k ‖∇u‖² for every displacement gradient. Isotropic: λ(tr ∇u)² + 2μ‖sym ∇u‖² ≤
/// max(3λ + 2μ, 2μ)‖∇u‖². Otherwise the largest eigenvalue of D (bounded by its largest
/// absolute row sum) times 2, because the engineering shear strains satisfy γ² ≤ 2(g_ij² + g_ji²).
pub fn modulus_bound(elastic: &Elastic, d: &Elasticity) -> f64 {
    match *elastic {
        Elastic::Isotropic { young, poisson } => {
            let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
            let mu = young / (2.0 * (1.0 + poisson));
            (3.0 * lambda + 2.0 * mu).max(2.0 * mu)
        }
        Elastic::EngineeringConstants { .. } => 2.0 * d.iter().map(|row| row.iter().map(|v| v.abs()).sum::<f64>()).fold(0.0, f64::max),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::isotropic;

    /// A distorted hexahedron (every face warped) and a tetrahedron.
    fn distorted() -> [[f64; 3]; 8] {
        [[0.0, 0.0, 0.0], [2.1, 0.2, -0.1], [2.3, 1.9, 0.2], [-0.2, 1.7, 0.1], [0.1, -0.1, 1.8], [2.0, 0.1, 2.2], [2.4, 2.2, 1.9], [0.2, 1.8, 2.1]]
    }

    const TET: [[f64; 3]; 4] = [[0.1, 0.0, 0.2], [1.3, 0.1, 0.0], [0.2, 1.1, 0.1], [0.3, 0.2, 0.9]];

    /// The volume form's b_I summed term by term as the form lists them, three cross products a
    /// term: the reference [`Hex8::hex8_deformed`] reorganises.
    fn form_gradients(xyz: &[[f64; 3]; 8]) -> [[f64; 3]; 8] {
        let mut b = [[0.0f64; 3]; 8];
        let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        for (i, j, k, c) in VOLUME_FORM {
            let (xi, xj, xk) = (xyz[i], xyz[j], xyz[k]);
            for (target, v) in [(i, cross(xj, xk)), (j, cross(xk, xi)), (k, cross(xi, xj))] {
                for d in 0..3 {
                    b[target][d] += c * v[d];
                }
            }
        }
        b.map(|v| v.map(|c| c / 12.0))
    }

    /// The shared cross products and per-corner sums give the form's gradients bit for bit, on
    /// distorted and random hexahedra and on axis-aligned ones, whose many equal products make
    /// zeros of both signs.
    #[test]
    fn deformed_gradients_are_the_forms_bit_for_bit() {
        let mut state = 7u64;
        let mut random = move || {
            state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut shapes = vec![distorted(), HEX_CORNERS, HEX_CORNERS.map(|c| c.map(|v| 2.0 * v - 1.0))];
        for _ in 0..200 {
            shapes.push(HEX_CORNERS.map(|c| c.map(|v| 3.0 * v + 0.2 * (random() - 0.5))));
        }
        for xyz in shapes {
            let bits = |b: [[f64; 3]; 8]| b.map(|v| v.map(f64::to_bits));
            assert_eq!(bits(Hex8::hex8_deformed(&xyz).unwrap().b), bits(form_gradients(&xyz)), "{xyz:?}");
        }
        let pairs: Vec<(usize, usize)> = (0..8).flat_map(|p| (p + 1..8).map(move |q| (p, q))).collect();
        assert!(pairs.iter().enumerate().all(|(n, (p, q))| pair_index(*p, *q) == n));
    }

    #[test]
    fn mean_strain_is_exact_for_linear_fields() {
        let a = [[0.003, -0.001, 0.002], [0.004, -0.002, 0.0005], [-0.0015, 0.001, 0.0025]];
        let field = |x: &[f64; 3]| [0, 1, 2].map(|i| 0.7 + a[i][0] * x[0] + a[i][1] * x[1] + a[i][2] * x[2]);
        let expected = [a[0][0], a[1][1], a[2][2], a[0][1] + a[1][0], a[0][2] + a[2][0], a[1][2] + a[2][1]];
        let hex = distorted();
        let h = Hex8::hex8(&hex).unwrap();
        let t = Tet4::tet4(&TET).unwrap();
        for (got, b, x) in
            [(h.strain(&hex.map(|x| field(&x))), h.b.as_slice(), hex.as_slice()), (t.strain(&TET.map(|x| field(&x))), t.b.as_slice(), TET.as_slice())]
        {
            assert!(got.iter().zip(expected).all(|(g, e)| (g - e).abs() < 1e-15), "{got:?}");
            // Σ_I b_iI = 0 and Σ_I b_iI x_jI = V δ_ij.
            let v: f64 = (0..3).map(|i| b.iter().zip(x).map(|(bi, xi)| bi[i] * xi[i]).sum::<f64>()).sum::<f64>() / 3.0;
            for i in 0..3 {
                assert!(b.iter().map(|bi| bi[i]).sum::<f64>().abs() < 1e-14);
                for j in 0..3 {
                    let s: f64 = b.iter().zip(x).map(|(bi, xi)| bi[i] * xi[j]).sum();
                    assert!((s - if i == j { v } else { 0.0 }).abs() < 1e-14, "{i} {j}: {s}");
                }
            }
        }
        assert!((h.nodal_volume.iter().sum::<f64>() - h.volume).abs() < 1e-14);
        // The volume of the distorted hexahedron agrees with the fem element's.
        let solid = crate::Solid::new(fieldcraft_model::ElementKind::Hex8, &hex, crate::Formulation::Full).unwrap();
        assert!((solid.volume().unwrap() - h.volume).abs() < 1e-14);
        let mut inverted = hex;
        inverted.swap(0, 1);
        inverted.swap(4, 5);
        assert!(Hex8::hex8(&inverted).is_err());
        assert!(Tet4::tet4(&[TET[1], TET[0], TET[2], TET[3]]).is_err());
        assert!(Hex8::hex8(&[[f64::NAN; 3]; 8]).is_err());
    }

    #[test]
    fn the_volume_form_is_the_integrated_one() {
        let mut c = [[[0.0f64; 8]; 8]; 8];
        for p in Rule::Hex(2).points() {
            let s = shape::hex8(p.xi);
            let d = s.derivatives();
            for i in 0..8 {
                for j in 0..8 {
                    for k in 0..8 {
                        let m = [[d[i][0], d[j][0], d[k][0]], [d[i][1], d[j][1], d[k][1]], [d[i][2], d[j][2], d[k][2]]];
                        c[i][j][k] += p.weight * crate::dense::det3(&m);
                    }
                }
            }
        }
        let mut table = Vec::new();
        for i in 0..8 {
            for j in i + 1..8 {
                for k in j + 1..8 {
                    if (12.0 * c[i][j][k]).abs() > 1e-9 {
                        assert!((12.0 * c[i][j][k] - (12.0 * c[i][j][k]).round()).abs() < 1e-12);
                        table.push((i, j, k, (12.0 * c[i][j][k]).round()));
                    }
                }
            }
        }
        assert_eq!(table, VOLUME_FORM.to_vec());
    }

    #[test]
    fn closed_form_gradients_equal_the_integrated_ones() {
        let hex = distorted();
        let (a, b) = (Hex8::hex8(&hex).unwrap(), Hex8::hex8_deformed(&hex).unwrap());
        assert!((a.volume - b.volume).abs() < 1e-14 * a.volume, "{} vs {}", a.volume, b.volume);
        for (x, y) in a.b.iter().flatten().zip(b.b.iter().flatten()) {
            assert!((x - y).abs() < 1e-14, "{x} vs {y}");
        }
        let mut inverted = hex;
        inverted.swap(0, 1);
        inverted.swap(2, 3);
        inverted.swap(4, 5);
        inverted.swap(6, 7);
        assert!(Hex8::hex8_deformed(&inverted).is_err());
        assert!(Hex8::hex8_deformed(&[[f64::NAN; 3]; 8]).is_err());
    }

    #[test]
    fn hourglass_vectors_are_orthogonal_to_linear_fields() {
        let hex = distorted();
        let h = Hex8::hex8(&hex).unwrap();
        for g in h.hourglass_vectors(&hex) {
            assert!(g.iter().sum::<f64>().abs() < 1e-14);
            for i in 0..3 {
                let s: f64 = g.iter().zip(&hex).map(|(gk, x)| gk * x[i]).sum();
                assert!(s.abs() < 1e-13, "{s}");
            }
        }
        // On a parallelepiped the shape vectors are the base vectors themselves.
        let cube: [[f64; 3]; 8] = HEX_CORNERS.map(|c| [c[0] + 0.3 * c[1], c[1], c[2] + 0.2 * c[0]]);
        let hc = Hex8::hex8(&cube).unwrap();
        let g = hc.hourglass_vectors(&cube);
        assert!(g[3].iter().zip(HEX_CORNERS).all(|(v, c)| (v - c[0] * c[1] * c[2]).abs() < 1e-14), "{g:?}");
    }

    /// The largest eigenvalue of M⁻¹K by power iteration on M^{-1/2} K M^{-1/2}.
    fn largest_frequency<const N: usize>(e: &UniformStrain<N>, d: &Elasticity, density: f64) -> f64 {
        let n = 3 * N;
        let mut k = vec![0.0; n * n];
        for a in 0..n {
            let mut u = [[0.0; 3]; N];
            u[a / 3][a % 3] = 1.0;
            let f = e.forces(&crate::material::stress(d, &e.strain(&u)));
            for b in 0..n {
                k[b * n + a] = f[b / 3][b % 3];
            }
        }
        let m: Vec<f64> = (0..n).map(|i| density * e.nodal_volume[i / 3]).collect();
        let mut x: Vec<f64> = (0..n).map(|i| 1.0 + (i as f64 * 0.37).sin()).collect();
        let mut lambda = 0.0;
        for _ in 0..20_000 {
            let y: Vec<f64> = (0..n).map(|i| (0..n).map(|j| k[i * n + j] * x[j] / m[j].sqrt()).sum::<f64>() / m[i].sqrt()).collect();
            let norm = y.iter().map(|v| v * v).sum::<f64>().sqrt();
            lambda = x.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>() / x.iter().map(|v| v * v).sum::<f64>();
            x = y.iter().map(|v| v / norm).collect();
        }
        lambda
    }

    #[test]
    fn frequency_bound_holds_and_is_exact_for_a_cube() {
        let cube: [[f64; 3]; 8] = HEX_CORNERS.map(|c| c.map(|v| 0.5 * (v + 1.0)));
        for nu in [0.0, 0.3, 0.45] {
            let elastic = Elastic::Isotropic { young: 1000.0, poisson: nu };
            let d = isotropic(1000.0, nu).unwrap();
            let k = modulus_bound(&elastic, &d);
            let h = Hex8::hex8(&cube).unwrap();
            let (bound, actual) = (h.frequency_bound(k, 2.0), largest_frequency(&h, &d, 2.0));
            assert!((actual - bound).abs() < 1e-6 * bound, "ν = {nu}: {actual} vs {bound}");
            let hex = distorted();
            let h = Hex8::hex8(&hex).unwrap();
            let t = Tet4::tet4(&TET).unwrap();
            for (bound, actual) in
                [(h.frequency_bound(k, 2.0), largest_frequency(&h, &d, 2.0)), (t.frequency_bound(k, 2.0), largest_frequency(&t, &d, 2.0))]
            {
                assert!(actual <= bound * (1.0 + 1e-9) && actual > 0.3 * bound, "ν = {nu}: {actual} vs {bound}");
            }
        }
        let ortho = Elastic::EngineeringConstants { young: [140e3, 10e3, 12e3], poisson: [0.3, 0.25, 0.4], shear: [5e3, 4e3, 3.8e3] };
        let d = crate::material::elasticity(&ortho).unwrap();
        let h = Hex8::hex8(&distorted()).unwrap();
        assert!(largest_frequency(&h, &d, 1.0) <= h.frequency_bound(modulus_bound(&ortho, &d), 1.0));
    }

    #[test]
    fn closed_form_eigenvalue() {
        assert_eq!(largest_eigenvalue(&[[2.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 1.0]]), 5.0);
        // [[2, 1, 0], [1, 2, 0], [0, 0, 1]] has eigenvalues 3, 1, 1.
        assert!((largest_eigenvalue(&[[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1.0]]) - 3.0).abs() < 1e-14);
        assert!((largest_eigenvalue(&[[4.0, 1.0, 2.0], [1.0, 3.0, 0.5], [2.0, 0.5, 5.0]]) - 6.831_254_430_698_204).abs() < 1e-9);
    }
}
