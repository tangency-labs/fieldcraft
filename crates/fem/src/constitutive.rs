//! The material point of finite-strain analysis: from the deformation gradient F (and the
//! material's history) to the stress and its consistent tangent.
//!
//! **Stress measure.** A material returns the **Kirchhoff stress** τ = J σ (J = det F, σ the
//! Cauchy stress) and the **spatial tangent** c of the Lie derivative of τ, L_v τ = c : d, with
//! d the rate of deformation. Both are 6-vectors and 6 × 6 matrices in the order xx, yy, zz,
//! xy, xz, yz, c taking engineering shears like the small-strain elasticity matrices (so the
//! linear elasticity D is c at F = I). This is the total Lagrangian formulation written with
//! spatial quantities: τ = F S Fᵀ and c_ijkl = F_iI F_jJ F_kK F_lL C_IJKL push the second
//! Piola–Kirchhoff stress S and its material tangent C = ∂S/∂E forward, so the internal forces
//! ∫ Bₓᵀ τ dV₀ and the stiffness ∫ (Bₓᵀ c Bₓ + ∇ₓN·τ∇ₓN) dV₀, with the shape function gradients
//! of the current configuration and integrals over the reference volume, are exactly the total
//! Lagrangian ∫ B₀ᵀ S dV₀ and its linearization (J. Bonet, R. D. Wood, *Nonlinear Continuum
//! Mechanics for Finite Element Analysis*, 2nd ed., Cambridge 2008, §7.4 and §9.2; G. A.
//! Holzapfel, *Nonlinear Solid Mechanics*, Wiley 2000, §6.6). Finite-strain plasticity (the
//! exponential return map of J. C. Simo, CMAME 99:61–112, 1992) gives τ and c directly, which is
//! why the interface is spatial.
//!
//! **History.** A material keeps [`MaterialModel::history_len`] numbers per integration point.
//! An evaluation reads the history committed at the end of the last converged increment and
//! writes the updated one to a separate slice; the solver commits it when the increment
//! converges, so a rejected iteration or a cut-back increment leaves no trace. The elastic
//! materials here keep none.
//!
//! **Materials:**
//! - **St Venant–Kirchhoff**: S = D E with E = ½(FᵀF − I) the Green–Lagrange strain and D the
//!   small-strain elasticity (isotropic or orthotropic in the global axes): an elastic material
//!   in a nonlinear step. Exact for large rotations with small strains; it softens without bound
//!   in large compression, so large strains want a hyperelastic model.
//! - **Neo-Hookean** (compressible, with the isochoric–volumetric split of P. J. Flory,
//!   *Thermodynamic relations for high elastic materials*, Trans. Faraday Soc. 57:829–838, 1961):
//!   W = μ/2 (Ī₁ − 3) + U(J), Ī₁ = J^(−2/3) tr b, U = κ/2 (J − 1)². Then
//!   τ = μ dev b̄ + J U′ I with b̄ = J^(−2/3) b, and
//!   c = ⅔ μ tr b̄ ℙ − ⅔ (τ_iso ⊗ I + I ⊗ τ_iso) + J (p̃ I ⊗ I − 2p 𝕀), with p = U′(J),
//!   p̃ = p + J U″(J), 𝕀 the symmetric fourth-order identity and ℙ = 𝕀 − ⅓ I ⊗ I (Holzapfel 2000,
//!   §6.6, eqs. 6.165–6.169 and 6.193, in Kirchhoff form). The split is what the mixed (hybrid)
//!   elements need: they take the isochoric part from [`MaterialModel::isochoric`] and the
//!   volumetric one from [`MaterialModel::volumetric`] with the element's own dilatation.
//! - **Mooney–Rivlin and Yeoh**: decoupled polynomials
//!   W̄ = C₁₀ (Ī₁ − 3) + C₀₁ (Ī₂ − 3) + C₂₀ (Ī₁ − 3)² + C₃₀ (Ī₁ − 3)³ with the neo-Hookean volume; their fictitious stress and
//!   elasticity from Holzapfel 2000, eqs. 6.193 and 6.168 ([`MaterialModel::Polynomial`]).
//! - **J2 plasticity** ([`crate::elastoplastic`]): Simo's multiplicative plasticity with the
//!   exponential return map in the principal logarithmic strains (isotropic hardening),
//!   Lagrangian logarithmic-strain return (finite kinematic hardening), or the
//!   additive radial return (small strain, [`MaterialModel::small_strain`]); decoupled too, the
//!   return deviatoric and the Hencky volume K/2 (ln J)².

use std::sync::Arc;

use fieldcraft_model::{Elastic, Hyperelastic};

use crate::dense::det3;
use crate::elastoplastic;
use crate::material::{Elasticity, elasticity, stress};
use crate::plasticity::J2;
use crate::rotation::{IDENTITY, Mat3, mul, transpose, voigt};
use crate::{Error, Result};

/// A 6 × 6 tangent in Voigt order, engineering-shear columns.
pub type Tangent = [[f64; 6]; 6];

/// The (row, column) of each Voigt component.
pub const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// The Kirchhoff stress and the spatial tangent of a material point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Response {
    pub stress: [f64; 6],
    pub tangent: Tangent,
}

/// The constitutive model of a solid in a nonlinear step.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialModel {
    StVenantKirchhoff(Box<Elasticity>),
    NeoHookean {
        shear: f64,
        bulk: f64,
    },
    /// A decoupled polynomial in the isochoric invariants,
    /// W̄ = C₁₀ (Ī₁ − 3) + C₀₁ (Ī₂ − 3) + C₂₀ (Ī₁ − 3)² + C₃₀ (Ī₁ − 3)³ with `c` = (C₁₀, C₀₁, C₂₀, C₃₀), and U = κ/2 (J − 1)²:
    /// Mooney–Rivlin and Yeoh.
    Polynomial {
        name: &'static str,
        c: [f64; 4],
        bulk: f64,
    },
    /// Von Mises plasticity (shared: one per material, whatever the element count).
    Plastic(Arc<J2>),
}

/// The 6 × 6 matrix of the symmetric fourth-order identity 𝕀 (engineering-shear columns).
pub(crate) fn identity4() -> Tangent {
    let mut i4 = [[0.0; 6]; 6];
    for k in 0..6 {
        i4[k][k] = if k < 3 { 1.0 } else { 0.5 };
    }
    i4
}

/// The matrix T with (F A Fᵀ) = T a for a symmetric A of Voigt vector a (tensor components):
/// T[(ij)][(IJ)] = F_iI F_jJ + F_iJ F_jI for I ≠ J and F_iI F_jI for I = J. A tangent pushes
/// forward as T C Tᵀ.
pub fn push_forward(f: &Mat3) -> Tangent {
    let mut t = [[0.0; 6]; 6];
    for (a, &(i, j)) in PAIRS.iter().enumerate() {
        for (b, &(ii, jj)) in PAIRS.iter().enumerate() {
            t[a][b] = if ii == jj { f[i][ii] * f[j][ii] } else { f[i][ii] * f[j][jj] + f[i][jj] * f[j][ii] };
        }
    }
    t
}

/// T C Tᵀ.
pub(crate) fn sandwich(t: &Tangent, c: &Tangent) -> Tangent {
    let mut tc = [[0.0; 6]; 6];
    for a in 0..6 {
        for b in 0..6 {
            tc[a][b] = (0..6).map(|k| t[a][k] * c[k][b]).sum();
        }
    }
    let mut out = [[0.0; 6]; 6];
    for a in 0..6 {
        for b in 0..6 {
            out[a][b] = (0..6).map(|k| tc[a][k] * t[b][k]).sum();
        }
    }
    out
}

/// det F, refusing an inverted or collapsed point.
fn jacobian(f: &Mat3) -> Result<f64> {
    let j = det3(f);
    if j > 0.0 && j.is_finite() { Ok(j) } else { Err(Error::Degenerate(format!("the material is inverted (J = det F = {j:.3e} is not positive)"))) }
}

impl MaterialModel {
    /// The model of a material's elastic behaviour in a nonlinear step: the hyperelastic option
    /// when the material has one (its isotropic constants give μ and κ), St Venant–Kirchhoff
    /// otherwise.
    pub fn from_model(elastic: &Elastic, hyperelastic: Option<&Hyperelastic>) -> Result<MaterialModel> {
        match hyperelastic {
            None => Ok(MaterialModel::StVenantKirchhoff(Box::new(elasticity(elastic)?))),
            Some(h) => {
                let (young, poisson) = crate::material::isotropic_constants(elastic, h.name())?;
                let (shear, bulk) = Hyperelastic::moduli(young, poisson);
                Ok(Self::hyperelastic(h, shear, bulk))
            }
        }
    }

    /// The hyperelastic model of an option with small-strain moduli μ and κ.
    fn hyperelastic(h: &Hyperelastic, shear: f64, bulk: f64) -> MaterialModel {
        match h {
            Hyperelastic::NeoHookean {} => MaterialModel::NeoHookean { shear, bulk },
            _ => MaterialModel::Polynomial { name: h.name(), c: h.constants(shear), bulk },
        }
    }

    /// The model of an element whose elasticity is `d` (the small-strain one) and whose material
    /// has the large-strain option `hyperelastic`: the hyperelastic models take μ = D₄₄ and
    /// κ = D₁₂ + ⅔ μ of an isotropic `d`, St Venant–Kirchhoff `d` itself.
    pub fn from_elasticity(d: &Elasticity, hyperelastic: Option<&Hyperelastic>) -> MaterialModel {
        match hyperelastic {
            None => MaterialModel::StVenantKirchhoff(Box::new(*d)),
            Some(h) => Self::hyperelastic(h, d[3][3], d[0][1] + 2.0 / 3.0 * d[3][3]),
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            MaterialModel::StVenantKirchhoff(_) => "St Venant–Kirchhoff",
            MaterialModel::NeoHookean { .. } => "neo-Hookean",
            MaterialModel::Polynomial { name, .. } => name,
            MaterialModel::Plastic(_) => "J2 plastic",
        }
    }

    /// History variables per integration point.
    pub fn history_len(&self) -> usize {
        match self {
            MaterialModel::Plastic(_) => elastoplastic::HISTORY,
            _ => 0,
        }
    }

    /// Whether the response depends on the path (the material keeps history), so the solver
    /// treats it with the inelastic kinematics and predictors.
    pub fn is_path_dependent(&self) -> bool {
        self.history_len() > 0
    }

    /// The equivalent plastic strain in one point's history; `None` for an elastic material.
    pub fn plastic_strain(&self, history: &[f64]) -> Option<f64> {
        match self {
            MaterialModel::Plastic(_) => Some(elastoplastic::value(history, elastoplastic::PLASTIC_STRAIN)),
            _ => None,
        }
    }

    /// Whether the strain energy splits into isochoric and volumetric parts (what mixed
    /// elements need).
    pub fn is_decoupled(&self) -> bool {
        !matches!(self, MaterialModel::StVenantKirchhoff(_))
    }

    /// The Kirchhoff stress and spatial tangent at deformation gradient `f` and (for
    /// temperature-dependent plasticity) temperature `temperature`, from the history
    /// `committed` (written, updated, to `updated`).
    pub fn response(&self, f: &Mat3, temperature: Option<f64>, committed: &[f64], updated: &mut [f64]) -> Result<Response> {
        updated.iter_mut().zip(committed).for_each(|(u, c)| *u = *c);
        let j = jacobian(f)?;
        match self {
            MaterialModel::StVenantKirchhoff(d) => {
                let c = mul(&transpose(f), f);
                let e = [0.5 * (c[0][0] - 1.0), 0.5 * (c[1][1] - 1.0), 0.5 * (c[2][2] - 1.0), c[0][1], c[0][2], c[1][2]];
                let s = stress(d, &e);
                let t = push_forward(f);
                Ok(Response { stress: std::array::from_fn(|a| (0..6).map(|b| t[a][b] * s[b]).sum()), tangent: sandwich(&t, d) })
            }
            MaterialModel::NeoHookean { .. } | MaterialModel::Polynomial { .. } | MaterialModel::Plastic(_) => {
                let iso = self.isochoric(f, temperature, committed, updated)?;
                let (p, dp) = self.volumetric(j, false).ok_or_else(|| Error::Invalid("volumetric response".into()))?;
                Ok(add_pressure(&iso, j, p, j * dp))
            }
        }
    }

    /// The isochoric part of a decoupled material (τ_iso = dev τ̄ and its tangent), from the
    /// history `committed`; the history that changes is written to `updated`, which the caller
    /// starts as a copy of `committed`.
    pub fn isochoric(&self, f: &Mat3, temperature: Option<f64>, committed: &[f64], updated: &mut [f64]) -> Result<Response> {
        let mu = match self {
            MaterialModel::NeoHookean { shear, .. } => *shear,
            MaterialModel::Plastic(j2) => return elastoplastic::finite_isochoric(j2, f, temperature, committed, updated),
            MaterialModel::Polynomial { c, .. } => return polynomial_isochoric(c, f),
            MaterialModel::StVenantKirchhoff(_) => {
                return Err(Error::NotSupported(format!("the isochoric–volumetric split of a {} material", self.name())));
            }
        };
        let j = jacobian(f)?;
        let scale = libm::pow(j, -2.0 / 3.0);
        let b = mul(f, &transpose(f));
        let bbar = voigt(&b).map(|v| scale * v);
        let tr = bbar[0] + bbar[1] + bbar[2];
        let mut tau = bbar.map(|v| mu * v);
        for v in tau.iter_mut().take(3) {
            *v -= mu * tr / 3.0;
        }
        let i4 = identity4();
        let mut c = [[0.0; 6]; 6];
        for a in 0..6 {
            for b in 0..6 {
                let (ia, ib) = (if a < 3 { 1.0 } else { 0.0 }, if b < 3 { 1.0 } else { 0.0 });
                let projector = i4[a][b] - ia * ib / 3.0;
                c[a][b] = 2.0 / 3.0 * mu * tr * projector - 2.0 / 3.0 * (tau[a] * ib + ia * tau[b]);
            }
        }
        Ok(Response { stress: tau, tangent: c })
    }

    /// U′(J) and U″(J) of a decoupled material's volumetric energy at dilatation J (with
    /// `small`, J = 1 + tr ε of small strain, where every material's is K/2 (J − 1)²).
    pub fn volumetric(&self, j: f64, small: bool) -> Option<(f64, f64)> {
        match self {
            MaterialModel::NeoHookean { bulk, .. } | MaterialModel::Polynomial { bulk, .. } => Some((bulk * (j - 1.0), *bulk)),
            MaterialModel::Plastic(j2) => elastoplastic::volumetric(j2, j, small),
            MaterialModel::StVenantKirchhoff(_) => None,
        }
    }

    /// Small strain: the stress and its tangent at the (mechanical) strain `strain`
    /// (engineering shears), from the history `committed` (written, updated, to `updated`);
    /// `bulk` false keeps the deviatoric part only (what a mixed element takes). St
    /// Venant–Kirchhoff is linear elasticity here; hyperelastic materials have no small-strain
    /// form.
    pub fn small_strain(&self, strain: &[f64; 6], temperature: Option<f64>, committed: &[f64], updated: &mut [f64], bulk: bool) -> Result<Response> {
        updated.iter_mut().zip(committed).for_each(|(u, c)| *u = *c);
        match self {
            MaterialModel::Plastic(j2) => elastoplastic::small(j2, strain, temperature, committed, updated, bulk),
            MaterialModel::StVenantKirchhoff(d) if bulk => Ok(Response { stress: stress(d, strain), tangent: **d }),
            _ => Err(Error::NotSupported(format!("the {} material at small strain{}", self.name(), if bulk { "" } else { " in mixed elements" }))),
        }
    }
}

/// The isochoric response of the polynomial energy W̄(Ī₁, Ī₂) with `c` = (C₁₀, C₀₁, C₂₀, C₃₀)
/// (Holzapfel 2000, §6.5–6.6, in Kirchhoff form): the fictitious stress
/// τ̄ = 2 (Ψ₁ + Ī₁ Ψ₂) b̄ − 2 Ψ₂ b̄², τ_iso = dev τ̄, and with the fictitious elasticity
/// c̄ = 4 (Ψ₁₁ + Ψ₂) b̄ ⊗ b̄ − 4 Ψ₂ b̄ ⊙ b̄ (eq. 6.193 pushed forward, Ψ₁₂ = Ψ₂₂ = 0 here),
/// c_iso = ℙ : c̄ : ℙ + ⅔ tr τ̄ ℙ − ⅔ (τ_iso ⊗ I + I ⊗ τ_iso) (eq. 6.168). With C₁₀ alone it is
/// the neo-Hookean response.
fn polynomial_isochoric(c: &[f64; 4], f: &Mat3) -> Result<Response> {
    let [c10, c01, c20, c30] = *c;
    let j = jacobian(f)?;
    let scale = libm::pow(j, -2.0 / 3.0);
    let b = mul(f, &transpose(f)).map(|row| row.map(|v| scale * v));
    let b2 = mul(&b, &b);
    let i1 = b[0][0] + b[1][1] + b[2][2];
    let x = i1 - 3.0;
    let (psi1, psi2, psi11) = (c10 + 2.0 * c20 * x + 3.0 * c30 * x * x, c01, 2.0 * c20 + 6.0 * c30 * x);
    let (bv, b2v) = (voigt(&b), voigt(&b2));
    let fict: [f64; 6] = std::array::from_fn(|a| 2.0 * (psi1 + i1 * psi2) * bv[a] - 2.0 * psi2 * b2v[a]);
    let tr = fict[0] + fict[1] + fict[2];
    let mut tau = fict;
    for v in tau.iter_mut().take(3) {
        *v -= tr / 3.0;
    }
    let mut cbar = [[0.0; 6]; 6];
    for (a, &(i, jj)) in PAIRS.iter().enumerate() {
        for (bb, &(k, l)) in PAIRS.iter().enumerate() {
            cbar[a][bb] = 4.0 * (psi11 + psi2) * bv[a] * bv[bb] - 2.0 * psi2 * (b[i][k] * b[jj][l] + b[i][l] * b[jj][k]);
        }
    }
    // ℙ : c̄ : ℙ = c̄ − ⅓ I ⊗ (I : c̄) − ⅓ (c̄ : I) ⊗ I + ⅑ (I : c̄ : I) I ⊗ I.
    let row: [f64; 6] = std::array::from_fn(|a| (0..3).map(|bb| cbar[a][bb]).sum());
    let col: [f64; 6] = std::array::from_fn(|bb| (0..3).map(|a| cbar[a][bb]).sum());
    let all: f64 = (0..3).map(|a| row[a]).sum();
    let i4 = identity4();
    let mut out = [[0.0; 6]; 6];
    for a in 0..6 {
        for bb in 0..6 {
            let (ia, ib) = (if a < 3 { 1.0 } else { 0.0 }, if bb < 3 { 1.0 } else { 0.0 });
            let projected = cbar[a][bb] - ia * col[bb] / 3.0 - row[a] * ib / 3.0 + ia * ib * all / 9.0;
            let projector = i4[a][bb] - ia * ib / 3.0;
            out[a][bb] = projected + 2.0 / 3.0 * tr * projector - 2.0 / 3.0 * (tau[a] * ib + ia * tau[bb]);
        }
    }
    Ok(Response { stress: tau, tangent: out })
}

/// An isochoric response with a pressure p acting on dilatation J: τ + J p I, and the tangent
/// c + J (p̃ I ⊗ I − 2p 𝕀) with p̃ = p + `j_dp` (J dp/dJ, zero when p is held).
pub fn add_pressure(iso: &Response, j: f64, p: f64, j_dp: f64) -> Response {
    let mut out = *iso;
    let i4 = identity4();
    for a in 0..3 {
        out.stress[a] += j * p;
    }
    for a in 0..6 {
        for b in 0..6 {
            let ii = if a < 3 && b < 3 { 1.0 } else { 0.0 };
            out.tangent[a][b] += j * ((p + j_dp) * ii - 2.0 * p * i4[a][b]);
        }
    }
    out
}

/// The logarithmic (Hencky) strain ln V = ½ ln b of a deformation gradient, as a strain 6-vector
/// with engineering shears; `None` when F is not finite or inverted.
pub fn logarithmic_strain(f: &Mat3) -> Option<[f64; 6]> {
    let (l, v) = crate::dense::symmetric_eigen(&mul(f, &transpose(f)))?;
    if !l.iter().all(|x| *x > 0.0) {
        return None;
    }
    let ln = l.map(|x| 0.5 * libm::log(x));
    let e = |i: usize, j: usize| (0..3).map(|k| ln[k] * v[i][k] * v[j][k]).sum::<f64>();
    Some([e(0, 0), e(1, 1), e(2, 2), 2.0 * e(0, 1), 2.0 * e(0, 2), 2.0 * e(1, 2)])
}

/// The deformation gradient F = I + Σ u_a ⊗ ∇₀N_a of nodal displacements and reference
/// gradients.
pub fn deformation_gradient(u: &[[f64; 3]], gradients: &[[f64; 3]]) -> Mat3 {
    let mut f = IDENTITY;
    for (ua, g) in u.iter().zip(gradients) {
        for i in 0..3 {
            for jj in 0..3 {
                f[i][jj] += ua[i] * g[jj];
            }
        }
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::inverse3;
    use crate::material::isotropic;

    fn deformation() -> Mat3 {
        [[1.2, 0.15, -0.05], [-0.1, 0.9, 0.2], [0.07, -0.12, 1.05]]
    }

    /// S = F⁻¹ τ F⁻ᵀ.
    fn pull_back_stress(f: &Mat3, tau: &[f64; 6]) -> [f64; 6] {
        crate::rotation::rotate(&inverse3(f).unwrap().0, tau)
    }

    /// The tangent's consistency: dS = C : dE with C the pull-back of c, checked by central
    /// differences of S along several directions dF (dE = sym(Fᵀ dF), engineering shears).
    fn check_tangent(m: &MaterialModel, f: &Mat3) {
        let r = m.response(f, None, &[], &mut []).unwrap();
        let (inv, _) = inverse3(f).unwrap();
        let back = sandwich(&push_forward(&inv), &r.tangent);
        let scale = back.iter().flatten().fold(0.0f64, |a, v| a.max(v.abs()));
        let h = 1e-6;
        for k in 0..9 {
            let mut df = [[0.0; 3]; 3];
            df[k / 3][k % 3] = 1.0;
            df[(k + 1) % 3][(k + 2) % 3] += 0.3;
            let at = |s: f64| {
                let fp: Mat3 = std::array::from_fn(|i| std::array::from_fn(|j| f[i][j] + s * df[i][j]));
                pull_back_stress(&fp, &m.response(&fp, None, &[], &mut []).unwrap().stress)
            };
            let (sp, sm) = (at(h), at(-h));
            let ftdf = mul(&transpose(f), &df);
            let de = [ftdf[0][0], ftdf[1][1], ftdf[2][2], ftdf[0][1] + ftdf[1][0], ftdf[0][2] + ftdf[2][0], ftdf[1][2] + ftdf[2][1]];
            for a in 0..6 {
                let numeric = (sp[a] - sm[a]) / (2.0 * h);
                let analytic: f64 = (0..6).map(|b| back[a][b] * de[b]).sum();
                assert!((numeric - analytic).abs() < 1e-6 * scale, "{} {k} {a}: {numeric} vs {analytic}", m.name());
            }
        }
    }

    #[test]
    fn tangents_are_consistent_with_the_stress() {
        let d = isotropic(1000.0, 0.3).unwrap();
        for m in [MaterialModel::StVenantKirchhoff(Box::new(d)), MaterialModel::NeoHookean { shear: 300.0, bulk: 5000.0 }] {
            check_tangent(&m, &deformation());
            check_tangent(&m, &IDENTITY);
            // Symmetric (hyperelastic): major symmetry of c.
            let c = m.response(&deformation(), None, &[], &mut []).unwrap().tangent;
            for a in 0..6 {
                for b in 0..6 {
                    assert!((c[a][b] - c[b][a]).abs() < 1e-9 * 5000.0, "{a} {b}");
                }
            }
        }
        // Mooney–Rivlin and Yeoh: consistent and symmetric; with C₁₀ alone the neo-Hookean model.
        for m in [
            MaterialModel::Polynomial { name: "Mooney–Rivlin", c: [120.0, 30.0, 0.0, 0.0], bulk: 5000.0 },
            MaterialModel::Polynomial { name: "Yeoh", c: [150.0, 0.0, -12.0, 1.5], bulk: 5000.0 },
        ] {
            check_tangent(&m, &deformation());
            check_tangent(&m, &IDENTITY);
            let c = m.response(&deformation(), None, &[], &mut []).unwrap().tangent;
            assert!((0..36).all(|k| (c[k / 6][k % 6] - c[k % 6][k / 6]).abs() < 1e-9 * 5000.0), "{} symmetric", m.name());
        }
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 5000.0 };
        let poly = MaterialModel::Polynomial { name: "neo", c: [150.0, 0.0, 0.0, 0.0], bulk: 5000.0 };
        let (a, b) = (neo.response(&deformation(), None, &[], &mut []).unwrap(), poly.response(&deformation(), None, &[], &mut []).unwrap());
        assert!(a.stress.iter().zip(&b.stress).all(|(x, y)| (x - y).abs() < 1e-11 * 300.0));
        assert!(a.tangent.iter().flatten().zip(b.tangent.iter().flatten()).all(|(x, y)| (x - y).abs() < 1e-10 * 5000.0));
        let ortho = crate::material::orthotropic([1000.0, 300.0, 500.0], [0.3, 0.2, 0.25], [100.0, 120.0, 90.0]).unwrap();
        check_tangent(&MaterialModel::StVenantKirchhoff(Box::new(ortho)), &deformation());
    }

    #[test]
    fn small_strain_limit_is_linear_elasticity_and_rotations_do_not_stress() {
        let (e, nu) = (1000.0, 0.3);
        let d = isotropic(e, nu).unwrap();
        let neo = MaterialModel::from_model(&Elastic::Isotropic { young: e, poisson: nu }, Some(&Hyperelastic::NeoHookean {})).unwrap();
        for m in [MaterialModel::StVenantKirchhoff(Box::new(d)), neo] {
            let c = m.response(&IDENTITY, None, &[], &mut []).unwrap().tangent;
            for a in 0..6 {
                for b in 0..6 {
                    assert!((c[a][b] - d[a][b]).abs() < 1e-9 * e, "{}: {a} {b}", m.name());
                }
            }
            // A rigid rotation (90° about an inclined axis) leaves the stress at zero.
            let (s, cth) = (1.0f64, 0.0f64);
            let axis = [1.0 / 3f64.sqrt(); 3];
            let k = crate::rotation::skew(&axis);
            let k2 = mul(&k, &k);
            let rot: Mat3 = std::array::from_fn(|i| std::array::from_fn(|j| IDENTITY[i][j] + s * k[i][j] + (1.0 - cth) * k2[i][j]));
            let tau = m.response(&rot, None, &[], &mut []).unwrap().stress;
            assert!(tau.iter().all(|v| v.abs() < 1e-12 * e), "{}: {tau:?}", m.name());
        }
    }

    #[test]
    fn uniaxial_stretch_of_neo_hookean_and_log_strain() {
        // Stretch λ along x with the lateral stretch μ chosen so that σ_yy = 0: check σ_xx
        // against W directly, σ = (μ/J^(5/3)) (λ² − (λ² + 2 l²)/3) + κ (J − 1).
        let (mu, kappa) = (300.0, 5000.0);
        let m = MaterialModel::NeoHookean { shear: mu, bulk: kappa };
        let (lam, l) = (1.5f64, 0.85f64);
        let f = [[lam, 0.0, 0.0], [0.0, l, 0.0], [0.0, 0.0, l]];
        let j = lam * l * l;
        let tau = m.response(&f, None, &[], &mut []).unwrap().stress;
        let i1 = lam * lam + 2.0 * l * l;
        let sxx = mu / j.powf(5.0 / 3.0) * (lam * lam - i1 / 3.0) + kappa * (j - 1.0);
        let syy = mu / j.powf(5.0 / 3.0) * (l * l - i1 / 3.0) + kappa * (j - 1.0);
        assert!((tau[0] / j - sxx).abs() < 1e-10 * mu && (tau[1] / j - syy).abs() < 1e-10 * mu, "{tau:?}");
        let e = logarithmic_strain(&f).unwrap();
        assert!((e[0] - lam.ln()).abs() < 1e-14 && (e[1] - l.ln()).abs() < 1e-14 && e[3].abs() < 1e-14, "{e:?}");
        assert!(m.response(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], None, &[], &mut []).is_err(), "inverted");
        assert!(MaterialModel::StVenantKirchhoff(Box::new(isotropic(1.0, 0.3).unwrap())).isochoric(&IDENTITY, None, &[], &mut []).is_err());
    }
}
