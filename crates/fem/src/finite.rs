//! Solid and 2D elements of nonlinear steps: total Lagrangian (finite strain) or with small
//! strain and an inelastic material; internal forces and the consistent tangent stiffness for
//! nodal displacements u, and the Cauchy stress and the strain at the integration points.
//!
//! **Formulation.** Total Lagrangian: the integrals run over the reference (undeformed) element
//! with its own quadrature, F = I + Σ u_a ⊗ ∇₀N_a at each point, and the material returns the
//! Kirchhoff stress τ and its spatial tangent c ([`crate::constitutive`]). With the shape
//! function gradients of the current configuration ∇ₓN_a = F⁻ᵀ ∇₀N_a,
//! - internal forces f_a = ∫ Bₓ,aᵀ τ dV₀,
//! - tangent K_ab = ∫ (Bₓ,aᵀ c Bₓ,b + (∇ₓN_a · τ ∇ₓN_b) I) dV₀ (material plus geometric, or
//!   initial-stress, stiffness),
//!
//! with Bₓ the small-strain operator of the current gradients. This is the push-forward of the
//! textbook ∫ B₀ᵀ S dV₀ and ∫ (B₀ᵀ C B₀ + ∇₀N·S∇₀N I) dV₀, identical term by term (Bonet & Wood
//! 2008, §9.2–9.4; K.-J. Bathe, *Finite Element Procedures*, 2nd ed. 2014, §6.3). A rigid motion
//! leaves τ at zero exactly, whatever its size.
//!
//! **Hybrid (mixed u–p) solids** with a decoupled material W = W̄(C̄) + U(J): the three-field
//! formulation of J. C. Simo, R. L. Taylor, K. S. Pister, *Variational and projection methods
//! for the volume constraint in finite deformation elasto-plasticity*, CMAME 51:177–208, 1985,
//! with the element's discontinuous pressure functions N_p (the linear hybrid elements' fields,
//! [`crate::hybrid`]): the dilatation θ = N_p·θ̂ is the projection of J, M θ̂ = ∫ N_p J dV₀ with
//! M = ∫ N_p N_pᵀ dV₀; the pressure p = N_p·p̂ the projection of U′(θ), M p̂ = ∫ N_p U′(θ) dV₀.
//! Then τ = τ_iso(F) + J p I, and with G = ∫ N_p J (∇ₓN)ᵀ dV₀ and H = ∫ N_p U″(θ) N_pᵀ dV₀,
//! K = ∫ (Bₓᵀ c* Bₓ + ∇ₓN·τ∇ₓN I) dV₀ + Gᵀ M⁻¹ H M⁻¹ G with c* = c_iso + J p (I ⊗ I − 2𝕀),
//! symmetric; at u = 0 it is the linear hybrid stiffness K_dev + κ G M⁻¹ Gᵀ. The constant
//! pressure case is the mean dilatation method of J. C. Nagtegaal, D. M. Parks, J. R. Rice
//! (CMAME 4:153–177, 1974) at finite strain. St Venant–Kirchhoff has no such split, so hybrid
//! elements of an elastic material are refused in nonlinear steps; J2 plasticity has it exactly
//! (the return is deviatoric, the volume Hencky's K/2 (ln J)², [`crate::elastoplastic`]).
//!
//! **Thermal strain** ([`Temperature`]): at finite strain the multiplicative split F = F_m F_θ
//! with the thermal stretches F_θ = diag(1 + αᵢ (T − T_ref)) at each point (T interpolated from
//! the nodes; V. A. Lubarda, *Constitutive theories based on the multiplicative decomposition of
//! deformation gradient: thermoelasticity, elastoplasticity, and biomechanics*, Appl. Mech. Rev.
//! 57(2):95–108, 2004): the material sees
//! F_m, the stress per reference volume is J_θ τ_m(F_m), and since F_θ is fixed over an
//! iteration the spatial tangent is J_θ c_m, so the loops run unchanged on the thermally
//! expanded volume J_θ dV₀ (the mixed elements project det F_m on it). The linear stretch
//! 1 + α ΔT is the engineering definition of α, so a free body grows by exactly α ΔT. At small
//! strain the split is additive, ε_m = ε − α ΔT. Reported strains are the total ones.
//!
//! **Small strain** ([`Kinematics::Small`], a step without nonlinear geometry with an inelastic
//! material): the same loops on the linear kinematics, ε = sym ∇u with the reference gradients,
//! f = ∫ Bᵀ σ dV₀ and K = ∫ Bᵀ 𝔻 B dV₀ with the material's consistent tangent 𝔻, and no
//! geometric stiffness. Mixed elements project the dilatation 1 + tr ε and the pressure K tr ε
//! the same way, which is the linear hybrid element of [`crate::hybrid`] while the material is
//! elastic.
//!
//! **2D elements**: plane strain (F₃₃ = 1) and axisymmetric (F₃₃ = r/R, the hoop stretch, with
//! the hoop column N/r of the current radius in Bₓ and a geometric term τ_θθ N_a N_b / r² on the
//! radial motion; Bonet & Wood 2008, §9.5). Plane stress (whose thickness change needs a local
//! iteration for σ_zz = 0) is refused.
//!
//! **Incompatible modes** (hex8, quad4): the enhanced assumed strain form of J. C. Simo,
//! F. Armero, *Geometrically non-linear enhanced strain mixed methods and the method of
//! incompatible modes*, IJNME 33(7):1413–1449, 1992, in which the incompatible-mode element is
//! the case of an enhanced deformation gradient F = I + Σ u_a ⊗ ∇₀N_a + Σ α_k ⊗ G̃_k with the
//! linear element's mapped mode gradients G̃_k = det J₀/det J · J₀⁻ᵀ ∂P_k/∂ξ (Taylor, Beresford
//! & Wilson 1976; in 2D with the mean-shifted hoop column of the axisymmetric element). The
//! modes join the loops above as extra "nodes": h = ∫ B̃ₓᵀ τ dV₀ are their forces and the
//! tangent over (u, α) has the same material and geometric terms. The modes α are element
//! unknowns, stored after the material history (the last converged α), and solved at every
//! evaluation by a local Newton iteration from them, α ← α − K_αα⁻¹ h, at most
//! `MAX_MODE_ITERATIONS` times to `MODE_TOLERANCE`; the material restarts from its
//! committed history at every local iteration, so plasticity is consistent too. The element
//! returns f* = f_u − K_uα K_αα⁻¹ h and K* = K_uu − K_uα K_αα⁻¹ K_αu, the consistent tangent at
//! the solved α (implicit function theorem). At u = 0 K* is the linear element's stiffness. A
//! K_αα that is not positive definite, or modes that do not converge, are an element error
//! (the increment cuts back).
//!
//! In large compression enhanced strain elements develop spurious hourglass modes (P. Wriggers,
//! S. Reese, *A note on enhanced strain methods for large deformations*, CMAME 135:201–209,
//! 1996). No stabilisation is adopted here (Glaser & Armero 1997 and Reese, Wriggers & Reddy
//! 2000 are the published remedies); instead, with nonlinear geometry and an elastic material
//! every tangent evaluation checks that K* is positive definite on the element's hourglass
//! modes (the patterns ξη, ηζ, ζξ, ξηζ of hex8 and ξη of quad4 along each axis, less their
//! linear part in the current configuration, so no constant-strain stiffness enters) and
//! refuses the state with an element error when it is not, so a step never converges to a
//! state that has silently lost its hourglass stiffness. The error cuts the increment back
//! like any element error; the loss belongs to the state, so smaller increments approach it
//! again and the step stops at the minimum increment, naming the element. The check is an
//! element-level criterion: it fires where the coarse element itself loses stiffness, not
//! where the material does (a neo-Hookean material stays strongly elliptic in compression), and
//! it would fire on displacement-based elements at large compression too (a fully integrated
//! hex8 cube loses that stiffness at the same stretch). (With plasticity a fully plastic
//! element may legitimately soften in such a mode, so the check is left to the global
//! tangent's negative pivots.)
//! Plane stress stays refused, with or without the modes.

use crate::constitutive::{MaterialModel, Response, add_pressure, logarithmic_strain};
use crate::dense::{Cholesky, Matrix, inverse3};
use crate::element::{Formulation, PointStress, Solid, accumulate_symmetric, b_block};
use crate::hybrid::pressure_basis;
use crate::plane::{Column, Plane, Theory, sub_block};
use crate::rotation::{IDENTITY, Mat3};
use crate::{Error, Result};

/// Internal forces (element degrees of freedom, node by node) and, when asked for, the tangent.
#[derive(Debug, Clone, PartialEq)]
pub struct Finite {
    pub force: Vec<f64>,
    pub tangent: Option<Matrix>,
}

/// The kinematics of an evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kinematics {
    /// Total Lagrangian, large displacements, rotations and strains.
    Finite,
    /// Small displacements and strains (ε = sym ∇u, the reference gradients, no geometric
    /// stiffness): only the material is nonlinear.
    Small,
}

/// The temperatures of an element in a nonlinear step: the linear expansion coefficients along
/// X, Y, Z (zero for a material that does not expand), the reference temperature at which the
/// thermal strain is zero, and the temperature at each node.
#[derive(Debug, Clone, PartialEq)]
pub struct Temperature {
    pub alpha: [f64; 3],
    pub reference: f64,
    pub nodal: Vec<f64>,
}

impl Temperature {
    fn check(&self, nodes: usize) -> Result<()> {
        if self.nodal.len() != nodes {
            return Err(Error::Invalid(format!("{} temperatures for a {nodes}-node element", self.nodal.len())));
        }
        if self.alpha.iter().chain(&self.nodal).chain([&self.reference]).any(|v| !v.is_finite()) {
            return Err(Error::Invalid("temperatures and expansion coefficients must be finite".into()));
        }
        Ok(())
    }
}

/// The temperature at a point and its thermal stretches ϑᵢ = 1 + αᵢ (T − T_ref).
#[derive(Debug, Clone, Copy)]
struct Heat {
    temperature: Option<f64>,
    stretch: [f64; 3],
}

impl Heat {
    /// At a point of shape function values `n`.
    fn at(t: Option<&Temperature>, n: &[f64]) -> Result<Heat> {
        let Some(t) = t else { return Ok(Heat { temperature: None, stretch: [1.0; 3] }) };
        let temperature: f64 = t.nodal.iter().zip(n).map(|(a, b)| a * b).sum();
        let stretch = t.alpha.map(|a| 1.0 + a * (temperature - t.reference));
        if !stretch.iter().all(|v| *v > 0.0 && v.is_finite()) {
            return Err(Error::Degenerate(format!("the thermal stretches {stretch:?} at T = {temperature} are not positive")));
        }
        Ok(Heat { temperature: Some(temperature), stretch })
    }

    /// J_θ = ϑ₁ ϑ₂ ϑ₃.
    fn volume(&self) -> f64 {
        self.stretch.iter().product()
    }

    /// The mechanical part F_m = F F_θ⁻¹ of a deformation gradient.
    fn mechanical(&self, f: &Mat3) -> Mat3 {
        std::array::from_fn(|i| std::array::from_fn(|j| f[i][j] / self.stretch[j]))
    }

    /// The mechanical part ε − α (T − T_ref) of a small strain.
    fn mechanical_strain(&self, e: &[f64; 6]) -> [f64; 6] {
        std::array::from_fn(|k| if k < 3 { e[k] - (self.stretch[k] - 1.0) } else { e[k] })
    }
}

/// What one integration point gives the element.
struct PointKinematics {
    /// The point's share of the reference volume (of the thermally expanded volume J_θ dV₀ at
    /// finite strain).
    weight: f64,
    /// The deformation gradient, and its mechanical part F_m = F F_θ⁻¹ the material sees.
    total: Mat3,
    f: Mat3,
    /// det F_m, or the dilatation 1 + tr ε_m of small strain (what a mixed element projects).
    j: f64,
    /// The mechanical small strain ε − ε_th (engineering shears); zero at finite strain.
    strain: [f64; 6],
    temperature: Option<f64>,
    /// ∇ₓN of every node (∇₀N at small strain).
    gradients: Vec<[f64; 3]>,
    /// Reference position.
    x0: [f64; 3],
}

/// The checked history slice of one integration point.
fn point_history(history: &[f64], len: usize, point: usize) -> &[f64] {
    history.get(point * len..(point + 1) * len).unwrap_or_default()
}

/// The small strain of a displacement gradient F − I (engineering shears).
fn small_strain(f: &Mat3) -> [f64; 6] {
    [f[0][0] - 1.0, f[1][1] - 1.0, f[2][2] - 1.0, f[0][1] + f[1][0], f[0][2] + f[2][0], f[1][2] + f[2][1]]
}

/// The response of one point's material with its history slices.
fn point_response(m: &MaterialModel, kin: Kinematics, p: &PointKinematics, committed: &[f64], updated: &mut [f64]) -> Result<Response> {
    match kin {
        Kinematics::Finite => m.response(&p.f, p.temperature, committed, updated),
        Kinematics::Small => m.small_strain(&p.strain, p.temperature, committed, updated, true),
    }
}

/// Add f_a += w Bₓ,aᵀ τ and, with `k`, the material tangent and (with `geometric`) the
/// geometric one.
fn add_point(force: &mut [f64], k: Option<&mut Matrix>, gradients: &[[f64; 3]], r: &Response, w: f64, geometric: bool) {
    for (a, g) in gradients.iter().enumerate() {
        let b = b_block(g);
        for c in 0..3 {
            if let Some(x) = force.get_mut(3 * a + c) {
                *x += w * (0..6).map(|s| b[s][c] * r.stress[s]).sum::<f64>();
            }
        }
    }
    let Some(k) = k else { return };
    accumulate_symmetric(&mut k.data, gradients, &r.tangent, w);
    if !geometric {
        return;
    }
    add_geometric(k, gradients, &r.stress, w);
}

/// K_ab += w (∇N_a · τ ∇N_b) I, the geometric (initial-stress) term.
fn add_geometric(k: &mut Matrix, gradients: &[[f64; 3]], stress: &[f64; 6], w: f64) {
    let tau = crate::rotation::tensor(stress);
    let n = k.n;
    for (a, ga) in gradients.iter().enumerate() {
        let tg: [f64; 3] = std::array::from_fn(|i| (0..3).map(|jj| tau[i][jj] * ga[jj]).sum());
        for (b, gb) in gradients.iter().enumerate() {
            let s = w * (tg[0] * gb[0] + tg[1] * gb[1] + tg[2] * gb[2]);
            for c in 0..3 {
                if let Some(x) = k.data.get_mut((3 * a + c) * n + 3 * b + c) {
                    *x += s;
                }
            }
        }
    }
}

fn inverted(kind: &str, j: f64) -> Error {
    Error::Degenerate(format!("a {kind} element is inverted at an integration point (J = det F = {j:.3e})"))
}

/// One displacement vector per node, all finite.
fn check_displacements<const N: usize>(u: &[[f64; N]], nodes: usize) -> Result<()> {
    if u.len() != nodes {
        return Err(Error::Invalid(format!("{} displacements for a {nodes}-node element", u.len())));
    }
    if u.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("displacements must be finite".into()));
    }
    Ok(())
}

/// Check the element's history slices against the material, its points and its enhanced
/// modes (stored after the points' blocks).
fn check_history(m: &MaterialModel, points: usize, modes: usize, committed: &[f64], updated: &[f64]) -> Result<usize> {
    let len = m.history_len();
    let expected = points * len + modes;
    if committed.len() != expected || updated.len() != expected {
        return Err(Error::Invalid(format!("{} history values for {points} points of {len} and {modes} enhanced modes", committed.len())));
    }
    Ok(len)
}

/// At most this many local Newton iterations solve an element's enhanced modes.
const MAX_MODE_ITERATIONS: usize = 25;
/// The enhanced modes are solved when their Newton correction is at most this fraction of the
/// largest nodal displacement or mode amplitude, plus `MODE_ABSOLUTE`.
const MODE_TOLERANCE: f64 = 1e-11;
/// The absolute part of the modes' tolerance, as a fraction of the element's size: above the
/// round-off of the Newton correction (about 1e-16 of the size, whatever the displacements), so
/// an element that barely moves still converges.
const MODE_ABSOLUTE: f64 = 1e-13;

/// The enhanced modes α of an incompatible-mode element at nodal displacements u (J. C. Simo,
/// F. Armero, IJNME 33:1413–1449, 1992, §4–5): starting from the last converged `alpha`, the
/// local Newton iteration α ← α − K_αα⁻¹ h solves h(u, α) = 0, and the condensed forces and
/// tangent are f* = f_u − K_uα K_αα⁻¹ h and K* = K_uu − K_uα K_αα⁻¹ K_αu, the consistent
/// tangent at the solved α. `evaluate` returns the forces and the tangent over the `nodal`
/// degrees of freedom and the modes after them; `scale` is the largest nodal displacement and
/// `floor` the absolute tolerance.
/// Returns the condensed forces and tangent (when `tangent`), what `evaluate` returned with
/// them, and the solved α. With `hourglass` (the element's hourglass modes and its name)
/// and `tangent`, the solved state is checked with [`check_hourglass`].
fn enhanced<P>(
    nodal: usize,
    alpha: &[f64],
    scale: f64,
    floor: f64,
    tangent: bool,
    hourglass: Option<(&[Vec<f64>], &str)>,
    mut evaluate: impl FnMut(&[f64]) -> Result<(Finite, P)>,
) -> Result<(Finite, P, Vec<f64>)> {
    let modes = alpha.len();
    let mut alpha = alpha.to_vec();
    for _ in 0..MAX_MODE_ITERATIONS {
        let (finite, points) = evaluate(&alpha)?;
        let k = finite.tangent.ok_or_else(|| Error::Invalid("the enhanced modes need the element tangent".into()))?;
        let chol = Cholesky::new(&sub_block(&k, nodal, nodal, modes))
            .map_err(|_| Error::Degenerate("the enhanced modes' stiffness K_αα of an incompatible-mode element is not positive definite".into()))?;
        let mut d = finite.force.get(nodal..).unwrap_or_default().to_vec();
        chol.solve(&mut d);
        if d.len() != modes || d.iter().any(|v| !v.is_finite()) {
            return Err(Error::Degenerate("the enhanced modes of an incompatible-mode element are not finite".into()));
        }
        let size = alpha.iter().fold(scale, |a, v| a.max(v.abs()));
        let step = d.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        alpha.iter_mut().zip(&d).for_each(|(a, v)| *a -= v);
        if step > MODE_TOLERANCE * size + floor {
            continue;
        }
        // f* = f_u − K_uα d with d = K_αα⁻¹ h, what is left of h.
        let force: Vec<f64> = (0..nodal)
            .map(|i| finite.force.get(i).copied().unwrap_or(0.0) - d.iter().enumerate().map(|(a, v)| k.get(i, nodal + a) * v).sum::<f64>())
            .collect();
        let condensed = tangent.then(|| {
            let mut out = sub_block(&k, 0, 0, nodal);
            for j in 0..nodal {
                let mut col: Vec<f64> = (0..modes).map(|a| k.get(nodal + a, j)).collect();
                chol.solve(&mut col);
                for i in 0..nodal {
                    out.add(i, j, -col.iter().enumerate().map(|(a, v)| k.get(i, nodal + a) * v).sum::<f64>());
                }
            }
            out
        });
        if let (Some((patterns, kind)), Some(condensed)) = (hourglass, condensed.as_ref()) {
            check_hourglass(condensed, patterns, kind)?;
        }
        return Ok((Finite { force, tangent: condensed }, points, alpha));
    }
    Err(Error::Degenerate(format!("the enhanced modes of an incompatible-mode element did not converge in {MAX_MODE_ITERATIONS} local iterations")))
}

/// The nonlinear evaluation of an incompatible-mode element whose `modes` follow the material
/// history in `committed` and `updated`: the modes solved and condensed by [`enhanced`] from
/// the stored ones, `evaluate` giving the uncondensed forces and tangent over the nodes and
/// the modes (as extra nodes after them in its `u`) with the material history only. The
/// hourglass check runs on elastic materials at finite strain when the tangent is asked for;
/// `hourglass` gives the element's hourglass modes at `u` and its name then.
#[allow(clippy::too_many_arguments)]
fn condensed<const N: usize, P>(
    m: &MaterialModel,
    u: &[[f64; N]],
    xyz: &[[f64; 3]],
    kin: Kinematics,
    modes: usize,
    committed: &[f64],
    updated: &mut [f64],
    tangent: bool,
    hourglass: impl FnOnce() -> Result<(Vec<Vec<f64>>, String)>,
    mut evaluate: impl FnMut(&[[f64; N]], &[f64], &mut [f64]) -> Result<(Finite, P)>,
) -> Result<(Finite, P)> {
    let split = committed.len().saturating_sub(modes);
    let (material, stored) = committed.split_at(split);
    let (updated, alpha_out) = updated.split_at_mut(split);
    let check = (tangent && kin == Kinematics::Finite && !m.is_path_dependent()).then(hourglass).transpose()?;
    let check = check.as_ref().map(|(patterns, kind)| (patterns.as_slice(), kind.as_str()));
    let mut extended = u.to_vec();
    let (scale, floor) = mode_scales(u, xyz);
    let (finite, points, alpha) = enhanced(N * u.len(), stored, scale, floor, tangent, check, |alpha| {
        extended.truncate(u.len());
        extended.extend(alpha.as_chunks::<N>().0);
        evaluate(&extended, material, updated)
    })?;
    if alpha_out.len() != alpha.len() {
        return Err(Error::Invalid(format!("{} enhanced modes stored for {}", alpha_out.len(), alpha.len())));
    }
    alpha_out.copy_from_slice(&alpha);
    Ok((finite, points))
}

/// The scales of the modes' tolerance: the largest nodal displacement, and `MODE_ABSOLUTE` of
/// the element's size.
fn mode_scales<const N: usize>(u: &[[f64; N]], xyz: &[[f64; 3]]) -> (f64, f64) {
    let size = (0..3)
        .map(|k| {
            let (lo, hi) = xyz.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p[k]), hi.max(p[k])));
            hi - lo
        })
        .fold(0.0f64, f64::max);
    (u.iter().flatten().fold(0.0f64, |a, v| a.max(v.abs())), MODE_ABSOLUTE * size)
}

/// The hourglass modes of a linear element in its current configuration: for each parent
/// pattern h (its values at the corners, which sum to zero) and each direction, the nodal
/// displacement field h minus its linear part. With the corners' current positions x_I, the
/// integrated gradients b_I = ∫ ∇N_I dv and the volume v, the mean gradient of a field u is
/// (1/v) Σ_I u_I ⊗ b_I, exact for every linear field because Σ_I b_I = 0 and
/// Σ_I x_I ⊗ b_I = v I (D. P. Flanagan, T. Belytschko, IJNME 17(5):679–706, 1981). The mode
/// Γ_I = h_I − (x_I − x̄) · g with g = (1/v) Σ_J h_J b_J has zero mean gradient and zero mean:
/// it is the part of h that no linear field reaches, so K restricted to the modes holds no
/// constant-strain stiffness. (Flanagan and Belytschko's γ_α = h_α − (1/v) Σ_i (h_α · x_i) b_i
/// is the same projection transposed: the covector that extracts these modes, orthogonal to
/// every linear field.) On a parallelepiped Γ = h.
fn hourglass_modes<const D: usize>(x: &[[f64; D]], b: &[[f64; D]], volume: f64, bases: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = x.len() as f64;
    let centre: [f64; D] = std::array::from_fn(|i| x.iter().map(|p| p[i]).sum::<f64>() / n);
    let mut out = Vec::with_capacity(bases.len() * D);
    for h in bases {
        let g: [f64; D] = std::array::from_fn(|i| h.iter().zip(b).map(|(hj, bj)| hj * bj[i]).sum::<f64>() / volume);
        let mode: Vec<f64> = h.iter().zip(x).map(|(hj, p)| hj - (0..D).map(|i| (p[i] - centre[i]) * g[i]).sum::<f64>()).collect();
        for dir in 0..D {
            out.push(mode.iter().flat_map(|m| std::array::from_fn::<f64, D, _>(|k| if k == dir { *m } else { 0.0 })).collect());
        }
    }
    out
}

/// The current positions of the first `N` corners `xyz` displaced by `u` (in `D` dimensions).
fn current<const N: usize, const D: usize>(xyz: &[[f64; 3]], u: &[[f64; D]]) -> [[f64; D]; N] {
    std::array::from_fn(|i| {
        let (p, d) = (xyz.get(i).copied().unwrap_or_default(), u.get(i).copied().unwrap_or([0.0; D]));
        std::array::from_fn(|k| p[k] + d[k])
    })
}

/// The hourglass modes of hex8 with corners `xyz` displaced by `u` ([`hourglass_modes`] of the
/// patterns ξη, ηζ, ζξ and ξηζ along X, Y and Z).
fn hex8_hourglass(xyz: &[[f64; 3]], u: &[[f64; 3]]) -> Result<Vec<Vec<f64>>> {
    let x: [[f64; 3]; 8] = current(xyz, u);
    let uniform = crate::uniform::Hex8::hex8_deformed(&x)?;
    let patterns: [fn(&[f64; 3]) -> f64; 4] = [|c| c[0] * c[1], |c| c[1] * c[2], |c| c[2] * c[0], |c| c[0] * c[1] * c[2]];
    let bases: Vec<Vec<f64>> = patterns.iter().map(|p| crate::shape::HEX_CORNERS.iter().map(p).collect()).collect();
    Ok(hourglass_modes(&x, &uniform.b, uniform.volume, &bases))
}

/// The hourglass modes of quad4 with corners `xyz` (in its plane) displaced by `u`
/// ([`hourglass_modes`] of the pattern ξη along x and y). The bilinear quadrilateral's integrated
/// gradients are exact in closed form, b_I = ½ (y_{I+1} − y_{I−1}, x_{I−1} − x_{I+1}), and its
/// area is their sum Σ_I x_I b_Ix; the axisymmetric element uses the same (r, z) geometry.
fn quad4_hourglass(xyz: &[[f64; 3]], u: &[[f64; 2]]) -> Result<Vec<Vec<f64>>> {
    let x: [[f64; 2]; 4] = current(xyz, u);
    let b: [[f64; 2]; 4] = std::array::from_fn(|i| {
        let (next, prev) = (x[(i + 1) % 4], x[(i + 3) % 4]);
        [0.5 * (next[1] - prev[1]), 0.5 * (prev[0] - next[0])]
    });
    let area: f64 = x.iter().zip(&b).map(|(p, q)| p[0] * q[0]).sum();
    if !(area > 0.0 && area.is_finite()) {
        return Err(Error::Degenerate(format!("a 2D quad4 element has a current area of {area:.3e}")));
    }
    Ok(hourglass_modes(&x, &b, area, &[vec![1.0, -1.0, 1.0, -1.0]]))
}

/// Whether `k` (its leading block over the patterns' degrees of freedom) is positive definite
/// on the span of `patterns`, by the symmetric part of Pᵀ K P.
fn positive_on(k: &Matrix, patterns: &[Vec<f64>]) -> bool {
    let n = patterns.len();
    let mut h = Matrix::zeros(n);
    for (i, a) in patterns.iter().enumerate() {
        for (j, b) in patterns.iter().enumerate() {
            let mut v = 0.0;
            for (r, x) in a.iter().enumerate().filter(|(_, x)| **x != 0.0) {
                for (c, y) in b.iter().enumerate().filter(|(_, y)| **y != 0.0) {
                    v += x * y * 0.5 * (k.get(r, c) + k.get(c, r));
                }
            }
            h.add(i, j, v);
        }
    }
    Cholesky::new(&h).is_ok()
}

/// Refuse a condensed tangent K* that is not positive definite on the element's hourglass
/// modes ([`hourglass_modes`]): the spurious instability of enhanced strain elements in large
/// compression (P. Wriggers, S. Reese, CMAME 135:201–209, 1996). It is an element-level
/// criterion: a fully integrated hex8 loses the same stiffness at large compression too, and an
/// incompatible-mode element stops there rather than risk a converged state that has
/// hourglassed. The error is `Degenerate`, so the increment cuts back.
fn check_hourglass(condensed: &Matrix, patterns: &[Vec<f64>], kind: &str) -> Result<()> {
    if positive_on(condensed, patterns) {
        return Ok(());
    }
    Err(Error::Degenerate(format!(
        "the {kind} incompatible-mode element lost its stiffness against hourglassing (the instability of enhanced strain elements in large compression, Wriggers & Reese 1996); refine the mesh there or use full integration, hybrid or quadratic elements"
    )))
}

impl Solid<'_> {
    /// Refuse what has no form here.
    fn check_nonlinear(&self, m: &MaterialModel, u: &[[f64; 3]]) -> Result<()> {
        check_displacements(u, self.ty.nodes)?;
        match self.formulation {
            Formulation::Hybrid { .. } if !m.is_decoupled() => Err(Error::NotSupported(format!(
                "hybrid elements of a {} material in nonlinear steps (the mixed formulation needs a volumetric–isochoric split of the strain energy; give the material the neo-Hookean option)",
                m.name()
            ))),
            _ => Ok(()),
        }
    }

    /// The length of the element's history in a nonlinear step: the material's per integration
    /// point, then the enhanced modes α of incompatible-mode elements.
    pub fn nonlinear_history_len(&self, m: &MaterialModel) -> usize {
        self.rule().points().len() * m.history_len() + self.dofs().internal.incompatible_modes
    }

    /// F, J and the gradients at every stiffness point: the current ones at finite strain, the
    /// reference ones at small strain. With incompatible modes `u` carries the modes α after
    /// the nodes, and their mapped gradients follow the nodes' (the enhanced F).
    fn kinematics(&self, u: &[[f64; 3]], kin: Kinematics, thermal: Option<&Temperature>) -> Result<Vec<PointKinematics>> {
        let n = self.ty.nodes;
        if let Some(t) = thermal {
            t.check(n)?;
        }
        let j0 = if self.formulation == Formulation::IncompatibleModes { Some(self.centre_jacobian()?) } else { None };
        let mut out = Vec::new();
        let mut reference: Vec<[f64; 3]> = Vec::with_capacity(n + 3);
        for p in self.rule().points() {
            let g = self.gradients(p.xi)?;
            reference.clear();
            reference.extend_from_slice(g.dx.get(..n).unwrap_or_default());
            if let Some(j0) = &j0 {
                reference.extend_from_slice(&self.incompatible_gradients(p.xi, g.det, j0));
            }
            let reference = reference.as_slice();
            let total = crate::constitutive::deformation_gradient(u, reference);
            let (weight, x0) = (g.det * p.weight, self.position_of(g.shape.values()));
            let heat = Heat::at(thermal, g.shape.values())?;
            if kin == Kinematics::Small {
                let strain = heat.mechanical_strain(&small_strain(&total));
                let j = 1.0 + strain[0] + strain[1] + strain[2];
                let gradients = reference.to_vec();
                out.push(PointKinematics { weight, total, f: total, j, strain, temperature: heat.temperature, gradients, x0 });
                continue;
            }
            let (inv, det) = inverse3(&total).ok_or_else(|| inverted(self.ty.kind.name(), 0.0))?;
            let j = det / heat.volume();
            if !(det > 0.0 && j.is_finite()) {
                return Err(inverted(self.ty.kind.name(), det));
            }
            let gradients = reference.iter().map(|d| std::array::from_fn(|i| (0..3).map(|jj| d[jj] * inv[jj][i]).sum())).collect();
            let f = heat.mechanical(&total);
            out.push(PointKinematics { weight: weight * heat.volume(), total, f, j, strain: [0.0; 6], temperature: heat.temperature, gradients, x0 });
        }
        Ok(out)
    }

    /// The effective response (Kirchhoff stress and tangent, or the stress of small strain) at
    /// every point, the forces, and the tangent when `tangent`.
    #[allow(clippy::too_many_arguments)]
    fn nonlinear_points(
        &self,
        m: &MaterialModel,
        u: &[[f64; 3]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> Result<(Finite, Vec<(PointKinematics, Response)>)> {
        if self.formulation == Formulation::SolidShell {
            return Err(Error::NotSupported("solid-shell nonlinear analysis".into()));
        }
        self.check_nonlinear(m, u)?;
        let modes = self.dofs().internal.incompatible_modes;
        check_history(m, self.rule().points().len(), modes, committed, updated)?;
        if modes == 0 {
            return self.evaluate_points(m, u, kin, thermal, committed, updated, tangent);
        }
        let hourglass = || Ok((hex8_hourglass(self.xyz, u)?, self.ty.kind.name().to_string()));
        condensed(m, u, self.xyz, kin, modes, committed, updated, tangent, hourglass, |u, material, updated| {
            self.evaluate_points(m, u, kin, thermal, material, updated, true)
        })
    }

    /// [`Self::nonlinear_points`] over the nodes and, with incompatible modes, the modes α
    /// after them in `u` (the forces and the tangent over both, uncondensed), with the material
    /// history only.
    #[allow(clippy::too_many_arguments)]
    fn evaluate_points(
        &self,
        m: &MaterialModel,
        u: &[[f64; 3]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> Result<(Finite, Vec<(PointKinematics, Response)>)> {
        let points = self.kinematics(u, kin, thermal)?;
        let len = m.history_len();
        let ndof = 3 * u.len();
        let mut force = vec![0.0; ndof];
        let mut k = tangent.then(|| Matrix::zeros(ndof));
        let responses: Vec<Response> = if let Formulation::Hybrid { .. } = self.formulation {
            self.mixed_responses(m, &points, kin, committed, updated, k.as_mut())?
        } else {
            let mut out = Vec::with_capacity(points.len());
            for (i, p) in points.iter().enumerate() {
                let updated_point = updated.get_mut(i * len..(i + 1) * len).unwrap_or_default();
                out.push(point_response(m, kin, p, point_history(committed, len, i), updated_point)?);
            }
            out
        };
        for (p, r) in points.iter().zip(&responses) {
            add_point(&mut force, k.as_mut(), &p.gradients, r, p.weight, kin == Kinematics::Finite);
        }
        Ok((Finite { force, tangent: k }, points.into_iter().zip(responses).collect()))
    }

    /// The three-field mixed responses τ_iso + J p I with their tangents (small strain:
    /// s + p I), and the dilatation stiffness Gᵀ M⁻¹ H M⁻¹ G added to `k`.
    fn mixed_responses(
        &self,
        m: &MaterialModel,
        points: &[PointKinematics],
        kin: Kinematics,
        committed: &[f64],
        updated: &mut [f64],
        k: Option<&mut Matrix>,
    ) -> Result<Vec<Response>> {
        let small = kin == Kinematics::Small;
        let len = m.history_len();
        let (np, frame) = (self.np(), self.frame());
        let basis: Vec<[f64; 4]> = points.iter().map(|p| pressure_basis(np, p.x0, &frame)).collect();
        let mut mm = Matrix::zeros(np);
        let mut theta = vec![0.0; np];
        for (p, nb) in points.iter().zip(&basis) {
            for a in 0..np {
                theta[a] += p.weight * nb[a] * p.j;
                for b in 0..np {
                    mm.add(a, b, p.weight * nb[a] * nb[b]);
                }
            }
        }
        let chol = Cholesky::new(&mm)?;
        chol.solve(&mut theta);
        let mut pressure = vec![0.0; np];
        let mut h = Matrix::zeros(np);
        for (p, nb) in points.iter().zip(&basis) {
            let th: f64 = (0..np).map(|a| nb[a] * theta[a]).sum();
            if !(th > 0.0 && th.is_finite()) {
                return Err(inverted(self.ty.kind.name(), th));
            }
            let (dp, ddp) = m.volumetric(th, small).ok_or_else(|| Error::Invalid("volumetric response".into()))?;
            for a in 0..np {
                pressure[a] += p.weight * nb[a] * dp;
                for b in 0..np {
                    h.add(a, b, p.weight * ddp * nb[a] * nb[b]);
                }
            }
        }
        chol.solve(&mut pressure);
        let mut out = Vec::with_capacity(points.len());
        for (i, (p, nb)) in points.iter().zip(&basis).enumerate() {
            let pr: f64 = (0..np).map(|a| nb[a] * pressure[a]).sum();
            let (c, u) = (point_history(committed, len, i), updated.get_mut(i * len..(i + 1) * len).unwrap_or_default());
            out.push(if small {
                let mut r = m.small_strain(&p.strain, p.temperature, c, u, false)?;
                r.stress.iter_mut().take(3).for_each(|s| *s += pr);
                r
            } else {
                add_pressure(&m.isochoric(&p.f, p.temperature, c, u)?, p.j, pr, 0.0)
            });
        }
        if let Some(k) = k {
            // A = M⁻¹ G, np × ndof; K += Aᵀ H A. At small strain G = ∫ N_p ∇₀N dV₀.
            let ndof = k.n;
            let mut a_cols = vec![vec![0.0; np]; ndof];
            for (p, nb) in points.iter().zip(&basis) {
                let jw = p.weight * if small { 1.0 } else { p.j };
                for (node, g) in p.gradients.iter().enumerate() {
                    for c in 0..3 {
                        if let Some(col) = a_cols.get_mut(3 * node + c) {
                            for (x, n) in col.iter_mut().zip(nb) {
                                *x += jw * n * g[c];
                            }
                        }
                    }
                }
            }
            for col in &mut a_cols {
                chol.solve(col);
            }
            let ha: Vec<Vec<f64>> = a_cols.iter().map(|col| h.mul_vec(col)).collect();
            for i in 0..ndof {
                for jj in i..ndof {
                    let v: f64 = a_cols.get(i).zip(ha.get(jj)).map_or(0.0, |(x, y)| x.iter().zip(y).map(|(p, q)| p * q).sum());
                    k.add(i, jj, v);
                    if i != jj {
                        k.add(jj, i, v);
                    }
                }
            }
        }
        Ok(out)
    }

    /// Internal forces and (with `tangent`) the tangent stiffness at nodal displacements `u`,
    /// with the element's history `committed` at the start of the increment written, updated,
    /// to `updated`: points × the material's history length, then the enhanced modes α of an
    /// incompatible-mode element ([`Self::nonlinear_history_len`] values in all).
    pub fn finite_strain(&self, m: &MaterialModel, u: &[[f64; 3]], committed: &[f64], updated: &mut [f64], tangent: bool) -> Result<Finite> {
        Ok(self.nonlinear_points(m, u, Kinematics::Finite, None, committed, updated, tangent)?.0)
    }

    /// As [`Self::finite_strain`] with either kinematics and the element's temperatures when
    /// it has some; with small strain an inelastic material on the linear kinematics,
    /// f = ∫ Bᵀ σ dV and K = ∫ Bᵀ 𝔻 B dV with the material's consistent tangent 𝔻.
    #[allow(clippy::too_many_arguments)]
    pub fn nonlinear(
        &self,
        m: &MaterialModel,
        u: &[[f64; 3]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> Result<Finite> {
        Ok(self.nonlinear_points(m, u, kin, thermal, committed, updated, tangent)?.0)
    }

    /// The Cauchy stress and the total strain (logarithmic at finite strain, ε at small strain)
    /// at every integration point (at its reference position).
    pub fn nonlinear_stresses(
        &self,
        m: &MaterialModel,
        u: &[[f64; 3]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
    ) -> Result<Vec<PointStress>> {
        let mut scratch = committed.to_vec();
        let (_, points) = self.nonlinear_points(m, u, kin, thermal, committed, &mut scratch, false)?;
        points
            .iter()
            .map(|(p, r)| {
                if kin == Kinematics::Small {
                    return Ok(PointStress { xyz: p.x0, strain: small_strain(&p.total), stress: r.stress });
                }
                let strain = logarithmic_strain(&p.total).ok_or_else(|| inverted(self.ty.kind.name(), p.j))?;
                Ok(PointStress { xyz: p.x0, strain, stress: r.stress.map(|v| v / p.j) })
            })
            .collect()
    }
}

impl Solid<'_> {
    /// The initial-stress (geometric) stiffness of a stress state given at the stiffness
    /// integration points: K_σ,ab = ∫ (∇N_a · σ ∇N_b) I dV over the reference element, the
    /// geometric term of the tangent with the reference gradients (linear buckling: Bathe 2014,
    /// §6.8.2; the same term as [`Solid::nonlinear`]'s at u = 0 with τ = σ).
    pub fn geometric_stiffness(&self, stress: &[[f64; 6]]) -> Result<Matrix> {
        if self.formulation == Formulation::SolidShell {
            return Err(Error::NotSupported("solid-shell buckling or prestressed perturbation".into()));
        }
        let points = self.rule().points();
        if stress.len() != points.len() {
            return Err(Error::Invalid(format!("{} stresses for {} integration points", stress.len(), points.len())));
        }
        let n = self.ty.nodes;
        let mut k = Matrix::zeros(3 * n);
        for (p, s) in points.iter().zip(stress) {
            if s.iter().any(|v| !v.is_finite()) {
                return Err(Error::Invalid("stresses must be finite".into()));
            }
            let g = self.gradients(p.xi)?;
            add_geometric(&mut k, g.dx.get(..n).unwrap_or_default(), s, g.det * p.weight);
        }
        Ok(k)
    }
}

impl Plane<'_> {
    /// The initial-stress stiffness of a 2D element (u, v per node) for the stresses of
    /// [`Plane::stresses_with`] at its stiffness points: the in-plane term and, axisymmetric,
    /// the hoop stress's σ_θθ N_a N_b / r² on the radial motion.
    pub fn geometric_stiffness(&self, stress: &[[f64; 6]]) -> Result<Matrix> {
        let rule = self.rule().points();
        if stress.len() != rule.len() {
            return Err(Error::Invalid(format!("{} stresses for {} integration points", stress.len(), rule.len())));
        }
        let n = self.ty.nodes;
        let mut k = Matrix::zeros(2 * n);
        for (p, s) in rule.iter().zip(stress) {
            if s.iter().any(|v| !v.is_finite()) {
                return Err(Error::Invalid("stresses must be finite".into()));
            }
            let at = self.at(p)?;
            let cols = at.g.get(..n).unwrap_or_default();
            for (a, ga) in cols.iter().enumerate() {
                for (b, gb) in cols.iter().enumerate() {
                    let v = ga[0] * gb[0] * s[0] + ga[1] * gb[1] * s[1] + (ga[0] * gb[1] + ga[1] * gb[0]) * s[3];
                    k.add(2 * a, 2 * b, at.dv * (v + s[2] * ga[2] * gb[2]));
                    k.add(2 * a + 1, 2 * b + 1, at.dv * v);
                }
            }
        }
        Ok(k)
    }
}

/// The kinematics of one point of a 2D element in a nonlinear step: F, det F_m, the volume share, the strain columns (current at finite strain), the small strain,
/// the temperature and the reference position.
struct PlaneKinematics {
    f: Mat3,
    /// det F_m (small strain: 1), and the dilatation the mixed element averages (small strain:
    /// 1 + tr ε).
    j: f64,
    dilatation: f64,
    dv: f64,
    current: Vec<[f64; 3]>,
    strain: [f64; 6],
    heat: Heat,
    x: [f64; 3],
}

/// One integration point of a 2D element in a nonlinear step: the total deformation gradient,
/// det F_m, the reference position and the response.
struct PlanePoint {
    f: Mat3,
    j: f64,
    x: [f64; 3],
    response: Response,
}

impl Plane<'_> {
    fn check_nonlinear(&self, m: &MaterialModel, u: &[[f64; 2]], kin: Kinematics) -> Result<()> {
        check_displacements(u, self.ty.nodes)?;
        if matches!(self.formulation, Formulation::Hybrid { .. }) && !m.is_decoupled() {
            return Err(Error::NotSupported(format!(
                "hybrid elements of a {} material in nonlinear steps (the mixed formulation needs a volumetric–isochoric split of the strain energy; give the material the neo-Hookean option)",
                m.name()
            )));
        }
        if let Theory::PlaneStress { .. } = self.theory {
            let why = match kin {
                Kinematics::Finite => "the thickness change needs a local iteration",
                Kinematics::Small => "σ_zz = 0 needs a local iteration of the inelastic material",
            };
            return Err(Error::NotSupported(format!(
                "plane stress elements in nonlinear steps ({why}; use plane strain, axisymmetric or solid elements)"
            )));
        }
        Ok(())
    }

    /// The length of the element's history in a nonlinear step, as [`Solid::nonlinear_history_len`].
    pub fn nonlinear_history_len(&self, m: &MaterialModel) -> usize {
        self.rule().points().len() * m.history_len() + self.dofs().internal.incompatible_modes
    }

    /// Per point: F, J, the reference position, the small strain and the response; the forces
    /// and the tangent.
    #[allow(clippy::too_many_arguments)]
    fn nonlinear_points(
        &self,
        m: &MaterialModel,
        u: &[[f64; 2]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> Result<(Finite, Vec<PlanePoint>)> {
        self.check_nonlinear(m, u, kin)?;
        if let Some(t) = thermal {
            t.check(self.ty.nodes)?;
        }
        let modes = self.dofs().internal.incompatible_modes;
        check_history(m, self.rule().points().len(), modes, committed, updated)?;
        if modes == 0 {
            return self.evaluate_points(m, u, None, kin, thermal, committed, updated, tangent);
        }
        let columns = self.incompatible_columns(&self.points()?)?;
        let hourglass = || Ok((quad4_hourglass(self.xyz, u)?, format!("2D {}", self.ty.kind.name())));
        condensed(m, u, self.xyz, kin, modes, committed, updated, tangent, hourglass, |u, material, updated| {
            self.evaluate_points(m, u, Some(&columns), kin, thermal, material, updated, true)
        })
    }

    /// [`Self::nonlinear_points`] over the nodes and, with the incompatible modes' strain
    /// `columns` at the points, the modes α after them in `u` (uncondensed), with the
    /// material history only.
    #[allow(clippy::too_many_arguments)]
    fn evaluate_points(
        &self,
        m: &MaterialModel,
        u: &[[f64; 2]],
        columns: Option<&[[Column; 2]]>,
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> Result<(Finite, Vec<PlanePoint>)> {
        let n = self.ty.nodes;
        let rule = self.rule().points();
        let len = m.history_len();
        let ndof = 2 * u.len();
        let mut force = vec![0.0; ndof];
        let mut k = tangent.then(|| Matrix::zeros(ndof));
        let axisymmetric = self.theory == Theory::Axisymmetric;
        let small = kin == Kinematics::Small;
        // The kinematics of every point first: the mixed element's pressure needs them all.
        let mut kinematics = Vec::with_capacity(rule.len());
        let mut cols: Vec<Column> = Vec::with_capacity(n + 2);
        for (i, p) in rule.iter().enumerate() {
            let at = self.at(p)?;
            cols.clear();
            cols.extend_from_slice(at.g.get(..n).unwrap_or_default());
            if let Some(m) = columns.and_then(|c| c.get(i)) {
                cols.extend_from_slice(m);
            }
            let cols = cols.as_slice();
            let mut f = IDENTITY;
            for (ua, g) in u.iter().zip(cols) {
                for r in 0..2 {
                    for c in 0..2 {
                        f[r][c] += ua[r] * g[c];
                    }
                }
                if axisymmetric {
                    f[2][2] += ua[0] * g[2];
                }
            }
            let heat = Heat::at(thermal, at.shape.values())?;
            let point = if small {
                let strain = heat.mechanical_strain(&small_strain(&f));
                let volume = 1.0 + strain[0] + strain[1] + strain[2];
                PlaneKinematics { f, j: 1.0, dilatation: volume, dv: at.dv, current: cols.to_vec(), strain, heat, x: at.x }
            } else {
                let det2 = f[0][0] * f[1][1] - f[0][1] * f[1][0];
                let j = det2 * f[2][2];
                if !(j > 0.0 && j.is_finite() && f[2][2] > 0.0) {
                    return Err(inverted(&format!("2D {}", self.ty.kind.name()), j));
                }
                let inv = [[f[1][1] / det2, -f[0][1] / det2], [-f[1][0] / det2, f[0][0] / det2]];
                let current: Vec<[f64; 3]> =
                    cols.iter().map(|g| [g[0] * inv[0][0] + g[1] * inv[1][0], g[0] * inv[0][1] + g[1] * inv[1][1], g[2] / f[2][2]]).collect();
                let j = j / heat.volume();
                PlaneKinematics { f, j, dilatation: j, dv: at.dv * heat.volume(), current, strain: [0.0; 6], heat, x: at.x }
            };
            kinematics.push(point);
        }
        let responses: Vec<Response> = if let Formulation::Hybrid { .. } = self.formulation {
            self.mean_dilatation(m, &kinematics, small, committed, updated, k.as_mut())?
        } else {
            let mut out = Vec::with_capacity(kinematics.len());
            for (i, p) in kinematics.iter().enumerate() {
                let updated_point = updated.get_mut(i * len..(i + 1) * len).unwrap_or_default();
                let history = point_history(committed, len, i);
                out.push(if small {
                    m.small_strain(&p.strain, p.heat.temperature, history, updated_point, true)?
                } else {
                    m.response(&p.heat.mechanical(&p.f), p.heat.temperature, history, updated_point)?
                });
            }
            out
        };
        let mut out = Vec::with_capacity(rule.len());
        for (p, r) in kinematics.into_iter().zip(responses) {
            let (dv, current) = (p.dv, &p.current);
            let tau = [r.stress[0], r.stress[1], r.stress[2], r.stress[3]];
            for (a, g) in current.iter().enumerate() {
                if let Some(x) = force.get_mut(2 * a) {
                    *x += dv * (g[0] * tau[0] + g[2] * tau[2] + g[1] * tau[3]);
                }
                if let Some(x) = force.get_mut(2 * a + 1) {
                    *x += dv * (g[1] * tau[1] + g[0] * tau[3]);
                }
            }
            if let Some(k) = k.as_mut() {
                let c4: [[f64; 4]; 4] = std::array::from_fn(|a| std::array::from_fn(|b| r.tangent[a][b]));
                crate::plane::accumulate(&mut k.data, current, &c4, dv);
                if !small {
                    for (a, ga) in current.iter().enumerate() {
                        for (b, gb) in current.iter().enumerate() {
                            let s = ga[0] * gb[0] * tau[0] + ga[1] * gb[1] * tau[1] + (ga[0] * gb[1] + ga[1] * gb[0]) * tau[3];
                            k.add(2 * a, 2 * b, dv * (s + tau[2] * ga[2] * gb[2]));
                            k.add(2 * a + 1, 2 * b + 1, dv * s);
                        }
                    }
                }
            }
            out.push(PlanePoint { f: p.f, j: p.j, x: p.x, response: r });
        }
        Ok((Finite { force, tangent: k }, out))
    }

    /// The mean dilatation responses of a hybrid 2D element (one constant pressure): θ = ∫ J dV₀
    /// / V₀ (hoop stretch included), p = U′(θ), τ = τ_iso + J p I (small strain: s + p I), and
    /// the dilatation stiffness U″(θ) g gᵀ / V₀ added to `k`, g = ∫ J ∇ₓ·N dV₀ (small strain
    /// ∫ ∇₀·N dV₀), the divergence taking the hoop column N/r too (the constant-pressure case of
    /// the solids' three-field form, Simo, Taylor & Pister 1985).
    fn mean_dilatation(
        &self,
        m: &MaterialModel,
        points: &[PlaneKinematics],
        small: bool,
        committed: &[f64],
        updated: &mut [f64],
        k: Option<&mut Matrix>,
    ) -> Result<Vec<Response>> {
        let len = m.history_len();
        let volume: f64 = points.iter().map(|p| p.dv).sum();
        if !(volume > 0.0 && volume.is_finite()) {
            return Err(inverted(&format!("2D {}", self.ty.kind.name()), volume));
        }
        let theta = points.iter().map(|p| p.dv * p.dilatation).sum::<f64>() / volume;
        if !(theta > 0.0 && theta.is_finite()) {
            return Err(inverted(&format!("2D {}", self.ty.kind.name()), theta));
        }
        let (pressure, stiffness) = m.volumetric(theta, small).ok_or_else(|| Error::Invalid("volumetric response".into()))?;
        let mut out = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            let (c, u) = (point_history(committed, len, i), updated.get_mut(i * len..(i + 1) * len).unwrap_or_default());
            out.push(if small {
                let mut r = m.small_strain(&p.strain, p.heat.temperature, c, u, false)?;
                r.stress.iter_mut().take(3).for_each(|s| *s += pressure);
                r
            } else {
                add_pressure(&m.isochoric(&p.heat.mechanical(&p.f), p.heat.temperature, c, u)?, p.j, pressure, 0.0)
            });
        }
        if let Some(k) = k {
            let mut g = vec![0.0; k.n];
            for p in points {
                crate::plane::divergence(&mut g, &p.current, p.dv * p.j);
            }
            crate::plane::add_outer(k, &g, stiffness / volume);
        }
        Ok(out)
    }

    /// Internal forces (u, v per node) and (with `tangent`) the tangent stiffness at nodal
    /// displacements `u`, as for [`Solid::finite_strain`].
    pub fn finite_strain(&self, m: &MaterialModel, u: &[[f64; 2]], committed: &[f64], updated: &mut [f64], tangent: bool) -> Result<Finite> {
        Ok(self.nonlinear_points(m, u, Kinematics::Finite, None, committed, updated, tangent)?.0)
    }

    /// As [`Solid::nonlinear`], on the 2D element.
    #[allow(clippy::too_many_arguments)]
    pub fn nonlinear(
        &self,
        m: &MaterialModel,
        u: &[[f64; 2]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> Result<Finite> {
        Ok(self.nonlinear_points(m, u, kin, thermal, committed, updated, tangent)?.0)
    }

    /// The Cauchy stress and the total strain at every integration point, as the solid
    /// elements' 6-vectors (axisymmetric: rr, axial, hoop, r-axial).
    pub fn nonlinear_stresses(
        &self,
        m: &MaterialModel,
        u: &[[f64; 2]],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
    ) -> Result<Vec<PointStress>> {
        let mut scratch = committed.to_vec();
        let (_, points) = self.nonlinear_points(m, u, kin, thermal, committed, &mut scratch, false)?;
        points
            .iter()
            .map(|p| {
                if kin == Kinematics::Small {
                    return Ok(PointStress { xyz: p.x, strain: small_strain(&p.f), stress: p.response.stress });
                }
                let strain = logarithmic_strain(&p.f).ok_or_else(|| inverted(&format!("2D {}", self.ty.kind.name()), p.j))?;
                Ok(PointStress { xyz: p.x, strain, stress: p.response.stress.map(|v| v / p.j) })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::isotropic;
    use fieldcraft_model::ElementKind;

    fn hex20() -> Vec<[f64; 3]> {
        let mut x: Vec<[f64; 3]> = crate::shape::HEX_CORNERS.iter().map(|c| [c[0] + 1.0, 0.5 * c[1] + 1.0, 0.8 * c[2]]).collect();
        let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
        for (a, b) in edges {
            x.push(std::array::from_fn(|k| 0.5 * (x[a][k] + x[b][k])));
        }
        // Distort a corner and a mid-side node.
        x[6][0] += 0.2;
        x[9][2] += 0.05;
        x
    }

    fn tet10() -> Vec<[f64; 3]> {
        let c = [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.2, 1.0, 0.1], [0.1, 0.2, 0.9]];
        let mut x = c.to_vec();
        for (a, b) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
            x.push(std::array::from_fn(|k| 0.5 * (c[a][k] + c[b][k])));
        }
        x
    }

    fn displacement(n: usize, scale: f64) -> Vec<[f64; 3]> {
        (0..n).map(|i| [scale * (i as f64 * 1.3).sin(), scale * (i as f64 * 0.7).cos(), scale * (i as f64 * 2.1).sin()]).collect()
    }

    /// The hex8 of corners [0, 1]³.
    fn unit_cube() -> Vec<[f64; 3]> {
        crate::shape::HEX_CORNERS.iter().map(|c| c.map(|v| 0.5 * (v + 1.0))).collect()
    }

    /// Central differences of the internal forces against the tangent.
    fn check_tangent(solid: &Solid, m: &MaterialModel, u: &[[f64; 3]]) {
        check_tangent_from(solid, m, u, Kinematics::Finite, &[], 1e-6);
    }

    /// The same from a committed history, with either kinematics.
    fn check_tangent_from(solid: &Solid, m: &MaterialModel, u: &[[f64; 3]], kin: Kinematics, committed: &[f64], tol: f64) {
        let mut scratch = committed.to_vec();
        let k = solid.nonlinear_points(m, u, kin, None, committed, &mut scratch, true).unwrap().0.tangent.unwrap();
        let scale = k.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        let n = 3 * u.len();
        let h = if kin == Kinematics::Small { 1e-9 } else { 1e-6 };
        for j in 0..n {
            let shifted = |s: f64| {
                let mut v = u.to_vec();
                v[j / 3][j % 3] += s;
                solid.nonlinear_points(m, &v, kin, None, committed, &mut committed.to_vec(), false).unwrap().0.force
            };
            let (p, q) = (shifted(h), shifted(-h));
            for i in 0..n {
                let numeric = (p[i] - q[i]) / (2.0 * h);
                assert!((numeric - k.get(i, j)).abs() < tol * scale, "{:?} {kin:?} {i} {j}: {numeric} vs {}", solid.ty.kind, k.get(i, j));
            }
        }
        for i in 0..n {
            for j in 0..n {
                assert!((k.get(i, j) - k.get(j, i)).abs() <= 1e-12 * scale, "symmetric");
            }
        }
    }

    #[test]
    fn solid_tangents_are_the_derivative_of_the_forces() {
        let d = isotropic(1000.0, 0.3).unwrap();
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 3000.0 };
        let svk = MaterialModel::StVenantKirchhoff(Box::new(d));
        let (x20, x10) = (hex20(), tet10());
        let x8: Vec<[f64; 3]> = x20[..8].to_vec();
        let x6 = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.1], [0., 0., 1.], [1.1, 0., 1.], [0., 1., 1.2]];
        let x4 = vec![x10[0], x10[1], x10[2], x10[3]];
        let x5 = vec![[0., 0., 0.], [1., 0., 0.1], [1.1, 1., 0.], [0., 1., 0.], [0.5, 0.4, 0.9]];
        for (kind, x) in [
            (ElementKind::Tet4, &x4),
            (ElementKind::Tet10, &x10),
            (ElementKind::Wedge6, &x6),
            (ElementKind::Pyramid5, &x5),
            (ElementKind::Hex8, &x8),
            (ElementKind::Hex20, &x20),
        ] {
            for m in [&svk, &neo] {
                let solid = Solid::new(kind, x, Formulation::Full).unwrap();
                check_tangent(&solid, m, &displacement(x.len(), 0.08));
            }
        }
        let reduced = Solid::new(ElementKind::Hex20, &x20, Formulation::Reduced).unwrap();
        check_tangent(&reduced, &neo, &displacement(20, 0.08));
        for (kind, x, np) in [(ElementKind::Hex8, &x8, 1), (ElementKind::Hex20, &x20, 4), (ElementKind::Tet10, &x10, 1)] {
            let hybrid = Solid::new(kind, x, Formulation::Hybrid { pressures: np }).unwrap();
            check_tangent(&hybrid, &neo, &displacement(x.len(), 0.08));
            let mr = MaterialModel::Polynomial { name: "Mooney–Rivlin", c: [100.0, 50.0, 0.0, 0.0], bulk: 3000.0 };
            let yeoh = MaterialModel::Polynomial { name: "Yeoh", c: [150.0, 0.0, -10.0, 2.0], bulk: 3000.0 };
            check_tangent(&hybrid, &mr, &displacement(x.len(), 0.08));
            check_tangent(&hybrid, &yeoh, &displacement(x.len(), 0.08));
            assert!(hybrid.finite_strain(&svk, &displacement(x.len(), 0.0), &[], &mut [], false).is_err(), "SVK has no split");
        }
        // Incompatible modes: the condensed tangent, the modes re-solved at every evaluation,
        // from rest and from stored modes.
        let im = Solid::new(ElementKind::Hex8, &x8, Formulation::IncompatibleModes).unwrap();
        for m in [&svk, &neo] {
            check_tangent_from(&im, m, &displacement(8, 0.08), Kinematics::Finite, &[0.0; 9], 1e-6);
            let mut stored = [0.0; 9];
            im.finite_strain(m, &displacement(8, 0.05), &[0.0; 9], &mut stored, false).unwrap();
            assert!(stored.iter().any(|a| *a != 0.0), "the modes are stored");
            check_tangent_from(&im, m, &displacement(8, 0.08), Kinematics::Finite, &stored, 1e-6);
        }
    }

    /// The incompatible-mode elements' local iteration: the modes solve h = 0 whatever they start
    /// from, hostile input is an error, and the history must hold the modes.
    #[test]
    fn incompatible_modes_solve_their_equations_and_refuse_hostile_input() {
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 3000.0 };
        let x8: Vec<[f64; 3]> = hex20()[..8].to_vec();
        let im = Solid::new(ElementKind::Hex8, &x8, Formulation::IncompatibleModes).unwrap();
        let u = displacement(8, 0.08);
        let (mut a, mut b) = ([0.0; 9], [0.0; 9]);
        let fa = im.finite_strain(&neo, &u, &[0.0; 9], &mut a, false).unwrap().force;
        let fb = im.finite_strain(&neo, &u, &[0.01; 9], &mut b, false).unwrap().force;
        for (x, y) in fa.iter().zip(&fb).chain(a.iter().zip(&b)) {
            assert!((x - y).abs() < 1e-9 * (1.0 + x.abs()), "the same state from another start: {x} vs {y}");
        }
        // The modes' own equations h = 0 at the stored modes.
        let mut extended = u.clone();
        extended.extend(a.as_chunks::<3>().0);
        let full = im.evaluate_points(&neo, &extended, Kinematics::Finite, None, &[], &mut [], false).unwrap().0;
        assert!(full.force[24..].iter().all(|h| h.abs() < 1e-9), "{:?}", &full.force[24..]);
        let mut bad = u.clone();
        bad[3][1] = f64::NAN;
        assert!(im.finite_strain(&neo, &bad, &[0.0; 9], &mut [0.0; 9], true).is_err());
        assert!(im.finite_strain(&neo, &u, &[f64::INFINITY; 9], &mut [0.0; 9], true).is_err());
        assert!(matches!(im.finite_strain(&neo, &u, &[], &mut [], true), Err(Error::Invalid(_))), "no room for the modes");
        let flip: Vec<[f64; 3]> = x8.iter().map(|p| [-2.0 * p[0], 0.0, 0.0]).collect();
        assert!(matches!(im.finite_strain(&neo, &flip, &[0.0; 9], &mut [0.0; 9], true), Err(Error::Degenerate(_))));
        assert_eq!(im.nonlinear_history_len(&neo), 9);
    }

    #[test]
    fn at_zero_displacement_the_tangent_is_the_linear_stiffness() {
        let (e, nu) = (1000.0, 0.45);
        let d = isotropic(e, nu).unwrap();
        let neo = MaterialModel::NeoHookean { shear: e / (2.0 * (1.0 + nu)), bulk: e / (3.0 * (1.0 - 2.0 * nu)) };
        let x20 = hex20();
        // Incompatible modes: the linear element's condensed stiffness.
        let x8: Vec<[f64; 3]> = x20[..8].to_vec();
        let im = Solid::new(ElementKind::Hex8, &x8, Formulation::IncompatibleModes).unwrap();
        let k = im.finite_strain(&neo, &[[0.0; 3]; 8], &[0.0; 9], &mut [0.0; 9], true).unwrap().tangent.unwrap();
        for (a, b) in k.data.iter().zip(&im.stiffness(&d).unwrap().data) {
            assert!((a - b).abs() < 1e-9 * e, "incompatible modes: {a} vs {b}");
        }
        for (f, m) in [
            (Formulation::Full, MaterialModel::StVenantKirchhoff(Box::new(d))),
            (Formulation::Full, neo.clone()),
            (Formulation::Hybrid { pressures: 4 }, neo),
        ] {
            let solid = Solid::new(ElementKind::Hex20, &x20, f).unwrap();
            let r = solid.finite_strain(&m, &vec![[0.0; 3]; 20], &[], &mut [], true).unwrap();
            let linear = solid.stiffness(&d).unwrap();
            let k = r.tangent.unwrap();
            for (a, b) in k.data.iter().zip(&linear.data) {
                assert!((a - b).abs() < 1e-9 * e, "{f:?}: {a} vs {b}");
            }
            assert!(r.force.iter().all(|v| v.abs() < 1e-12 * e));
        }
    }

    /// J2 plasticity in full and hybrid elements, at finite and small strain: the tangent is
    /// the derivative of the forces from a plastic history too, and the plastic points are
    /// really plastic.
    #[test]
    fn plastic_tangents_are_the_derivative_of_the_forces() {
        for hardening in [fieldcraft_model::Hardening::Isotropic, fieldcraft_model::Hardening::Kinematic] {
            check_plastic_solid_tangents(hardening);
        }
    }

    fn check_plastic_solid_tangents(hardening: fieldcraft_model::Hardening) {
        use crate::elastoplastic::{HISTORY, PLASTIC_STRAIN};
        use fieldcraft_model::{Elastic, Plastic};
        let j2 = crate::plasticity::J2::from_model(
            &Elastic::Isotropic { young: 1000.0, poisson: 0.3 },
            Some(&Plastic::Bilinear { yield_stress: 5.0, hardening_modulus: 50.0, hardening }),
            None,
        )
        .unwrap();
        let m = MaterialModel::Plastic(std::sync::Arc::new(j2));
        let (x20, x10) = (hex20(), tet10());
        let x8: Vec<[f64; 3]> = x20[..8].to_vec();
        let x4 = x10[..4].to_vec();
        let x6 = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]];
        for (kind, x, f) in [
            (ElementKind::Tet4, &x4, Formulation::Full),
            (ElementKind::Wedge6, &x6, Formulation::Full),
            (ElementKind::Hex20, &x20, Formulation::Full),
            (ElementKind::Tet10, &x10, Formulation::Hybrid { pressures: 1 }),
            (ElementKind::Hex8, &x8, Formulation::Full),
            (ElementKind::Tet10, &x10, Formulation::Full),
            (ElementKind::Hex20, &x20, Formulation::Reduced),
            (ElementKind::Hex8, &x8, Formulation::Hybrid { pressures: 1 }),
            (ElementKind::Hex20, &x20, Formulation::Hybrid { pressures: 4 }),
        ] {
            let solid = Solid::new(kind, x, f).unwrap();
            let points = solid.rule().points().len();
            for (kin, scale) in [(Kinematics::Finite, 0.05), (Kinematics::Small, 0.05)] {
                // A first increment from the virgin state, committed, then a second from it.
                let mut committed = vec![0.0; points * HISTORY];
                let first = displacement(x.len(), 0.5 * scale);
                solid.nonlinear_points(&m, &first, kin, None, &vec![0.0; points * HISTORY], &mut committed, false).unwrap();
                assert!((0..points).all(|p| committed[p * HISTORY + PLASTIC_STRAIN] > 0.0), "{kind:?} {f:?} {kin:?} yields everywhere");
                check_tangent_from(&solid, &m, &displacement(x.len(), scale), kin, &committed, 1e-6);
            }
        }
        // Incompatible modes: the modes are re-solved from the committed ones, the material from
        // its committed history at every local iteration.
        let im = Solid::new(ElementKind::Hex8, &x8, Formulation::IncompatibleModes).unwrap();
        let len = 8 * HISTORY + 9;
        for kin in [Kinematics::Finite, Kinematics::Small] {
            let mut committed = vec![0.0; len];
            im.nonlinear_points(&m, &displacement(8, 0.025), kin, None, &vec![0.0; len], &mut committed, false).unwrap();
            assert!((0..8).all(|p| committed[p * HISTORY + PLASTIC_STRAIN] > 0.0), "incompatible modes {kin:?} yield everywhere");
            check_tangent_from(&im, &m, &displacement(8, 0.05), kin, &committed, 1e-6);
        }
    }

    #[test]
    fn plastic_plane_strain_and_axisymmetric_tangents_are_the_derivative_of_the_forces() {
        for hardening in [fieldcraft_model::Hardening::Isotropic, fieldcraft_model::Hardening::Kinematic] {
            check_plastic_plane_tangents(hardening);
        }
    }

    fn check_plastic_plane_tangents(hardening: fieldcraft_model::Hardening) {
        use crate::elastoplastic::{HISTORY, PLASTIC_STRAIN};
        use fieldcraft_model::{Elastic, Plastic};
        let j2 = crate::plasticity::J2::from_model(
            &Elastic::Isotropic { young: 1000.0, poisson: 0.3 },
            Some(&Plastic::Bilinear { yield_stress: 5.0, hardening_modulus: 50.0, hardening }),
            None,
        )
        .unwrap();
        let m = MaterialModel::Plastic(std::sync::Arc::new(j2));
        let quad8 = [
            [1.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [2.2, 1.0, 0.0],
            [0.9, 1.1, 0.0],
            [1.5, 0.0, 0.0],
            [2.1, 0.55, 0.0],
            [1.5, 1.1, 0.0],
            [0.95, 0.5, 0.0],
        ];
        for theory in [Theory::PlaneStrain { thickness: 2.0 }, Theory::Axisymmetric] {
            for (kind, x) in [(ElementKind::Quad8, &quad8[..]), (ElementKind::Tri3, &quad8[..3])] {
                let plane = Plane::new(kind, x, Formulation::Full, theory).unwrap();
                let points = plane.rule().points().len();
                for kin in [Kinematics::Finite, Kinematics::Small] {
                    let u = |scale: f64| -> Vec<[f64; 2]> { displacement(x.len(), scale).iter().map(|v| [v[0], v[1]]).collect() };
                    let mut committed = vec![0.0; points * HISTORY];
                    plane.nonlinear_points(&m, &u(0.03), kin, None, &vec![0.0; points * HISTORY], &mut committed, false).unwrap();
                    assert!((0..points).all(|p| committed[p * HISTORY + PLASTIC_STRAIN] > 0.0), "{theory:?} {kind:?} {kin:?}");
                    let at = u(0.06);
                    let k = plane.nonlinear_points(&m, &at, kin, None, &committed, &mut committed.clone(), true).unwrap().0.tangent.unwrap();
                    let scale = k.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                    let h = 1e-8;
                    for j in 0..2 * x.len() {
                        let shifted = |s: f64| {
                            let mut v = at.clone();
                            v[j / 2][j % 2] += s;
                            plane.nonlinear_points(&m, &v, kin, None, &committed, &mut committed.clone(), false).unwrap().0.force
                        };
                        let (p, q) = (shifted(h), shifted(-h));
                        for i in 0..2 * x.len() {
                            let numeric = (p[i] - q[i]) / (2.0 * h);
                            assert!(
                                (numeric - k.get(i, j)).abs() < 1e-6 * scale,
                                "{theory:?} {kind:?} {kin:?} {i} {j}: {numeric} vs {}",
                                k.get(i, j)
                            );
                        }
                    }
                }
            }
        }
        let stress = Plane::new(ElementKind::Quad4, &quad8[..4], Formulation::Full, Theory::PlaneStress { thickness: 1.0 }).unwrap();
        assert!(matches!(
            stress.nonlinear(&m, &[[0.0; 2]; 4], Kinematics::Small, None, &[0.0; 4 * HISTORY], &mut [0.0; 4 * HISTORY], true),
            Err(Error::NotSupported(_))
        ));
    }

    #[test]
    fn logarithmic_kinematic_simple_shear_is_smooth_and_step_convergent() {
        use crate::elastoplastic::{HISTORY, PLASTIC_STRAIN};
        use fieldcraft_model::{Elastic, Hardening, Plastic};
        let j2 = crate::plasticity::J2::from_model(
            &Elastic::Isotropic { young: 200_000.0, poisson: 0.3 },
            Some(&Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 2000.0, hardening: Hardening::Kinematic }),
            None,
        )
        .unwrap();
        let m = MaterialModel::Plastic(std::sync::Arc::new(j2));
        let x: Vec<_> = crate::shape::HEX_CORNERS.iter().map(|p| p.map(|v| 0.5 * (v + 1.0))).collect();
        let solid = Solid::new(ElementKind::Hex8, &x, Formulation::Full).unwrap();
        let mut terminal = Vec::new();
        for steps in [50, 100, 200] {
            let mut committed = vec![0.0; 8 * HISTORY];
            let mut curve = Vec::new();
            for step in 1..=steps {
                let gamma = 5.0 * step as f64 / steps as f64;
                let u: Vec<_> = x.iter().map(|p| [gamma * p[1], 0.0, 0.0]).collect();
                let mut updated = committed.clone();
                let (element, points) = solid.nonlinear_points(&m, &u, Kinematics::Finite, None, &committed, &mut updated, false).unwrap();
                let shear = points[0].1.stress[3];
                for (_, r) in &points {
                    assert!((r.stress[3] - shear).abs() < 1e-8);
                    assert!((r.stress[0] + r.stress[1] + r.stress[2]).abs() < 1e-8, "isochoric");
                }
                for axis in 0..3 {
                    assert!(element.force.as_chunks::<3>().0.iter().map(|f| f[axis]).sum::<f64>().abs() < 1e-8, "balanced forces");
                }
                assert!(updated[PLASTIC_STRAIN] >= committed[PLASTIC_STRAIN]);
                if steps == 50 && step == steps {
                    check_tangent_from(&solid, &m, &u, Kinematics::Finite, &committed, 1e-6);
                }
                curve.push(shear);
                committed = updated;
            }
            // The logarithmic model peaks smoothly, then decreases; monotonicity to γ = 5
            // is not a property even of its elastic Hencky limit.
            let peak = curve.iter().enumerate().max_by(|(_, a), (_, b)| a.total_cmp(b)).unwrap().0;
            assert!((2.6..2.9).contains(&(5.0 * (peak + 1) as f64 / steps as f64)));
            assert!(curve[..=peak].windows(2).all(|w| w[1] >= w[0]));
            assert!(curve[peak..].windows(2).all(|w| w[1] <= w[0]));
            terminal.push(*curve.last().unwrap());
        }
        assert!((terminal[2] - terminal[1]).abs() < (terminal[1] - terminal[0]).abs());
        assert!((terminal[2] - 870.981_087_088_643).abs() < 1e-3);
    }

    /// A temperature field: the tangent is the derivative of the forces (finite and small strain,
    /// full and hybrid, elastic and plastic), and a free expansion u = (F_θ − I) X of a uniform
    /// temperature stresses nothing.
    #[test]
    fn thermal_tangents_and_free_expansion() {
        use crate::elastoplastic::HISTORY;
        use fieldcraft_model::{Elastic, Hardening, Plastic};
        let j2 = crate::plasticity::J2::from_model(
            &Elastic::Isotropic { young: 1000.0, poisson: 0.3 },
            Some(&Plastic::Bilinear { yield_stress: 5.0, hardening_modulus: 50.0, hardening: Hardening::Isotropic }),
            None,
        )
        .unwrap();
        let plastic = MaterialModel::Plastic(std::sync::Arc::new(j2));
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 3000.0 };
        let x20 = hex20();
        let field = |n: usize| Temperature {
            alpha: [1e-3, 2e-3, 1.5e-3],
            reference: 20.0,
            nodal: (0..n).map(|i| 20.0 + 40.0 * (i as f64 * 0.37).sin()).collect(),
        };
        for f in [Formulation::Full, Formulation::Hybrid { pressures: 4 }] {
            let solid = Solid::new(ElementKind::Hex20, &x20, f).unwrap();
            let points = solid.rule().points().len();
            let t = field(20);
            for (m, kin, len) in [(&neo, Kinematics::Finite, 0), (&plastic, Kinematics::Finite, HISTORY), (&plastic, Kinematics::Small, HISTORY)] {
                let u = displacement(20, 0.03);
                let committed = vec![0.0; points * len];
                let k = solid.nonlinear(m, &u, kin, Some(&t), &committed, &mut committed.clone(), true).unwrap().tangent.unwrap();
                let scale = k.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                let h = if kin == Kinematics::Small { 1e-9 } else { 1e-7 };
                for j in 0..60 {
                    let shifted = |s: f64| {
                        let mut v = u.clone();
                        v[j / 3][j % 3] += s;
                        solid.nonlinear(m, &v, kin, Some(&t), &committed, &mut committed.clone(), false).unwrap().force
                    };
                    let (p, q) = (shifted(h), shifted(-h));
                    for i in 0..60 {
                        let numeric = (p[i] - q[i]) / (2.0 * h);
                        assert!((numeric - k.get(i, j)).abs() < 1e-6 * scale, "{f:?} {} {kin:?} {i} {j}: {numeric} vs {}", m.name(), k.get(i, j));
                    }
                }
            }
            // Free expansion by a uniform 150 degrees: no force, no stress.
            let uniform = Temperature { alpha: [1e-3, 2e-3, 1.5e-3], reference: 20.0, nodal: vec![170.0; 20] };
            let stretch = uniform.alpha.map(|a| 1.0 + 150.0 * a);
            let u: Vec<[f64; 3]> = x20.iter().map(|x| std::array::from_fn(|k| (stretch[k] - 1.0) * x[k])).collect();
            let r = solid.nonlinear(&neo, &u, Kinematics::Finite, Some(&uniform), &[], &mut [], false).unwrap();
            assert!(r.force.iter().all(|v| v.abs() < 1e-10 * 3000.0), "{f:?}: {:?}", r.force);
            let s = solid.nonlinear_stresses(&neo, &u, Kinematics::Finite, Some(&uniform), &[]).unwrap();
            assert!(s.iter().all(|p| p.stress.iter().all(|v| v.abs() < 1e-10 * 3000.0)), "{f:?}");
            assert!((s[0].strain[1] - stretch[1].ln()).abs() < 1e-12, "total logarithmic strain: {:?}", s[0].strain);
        }
    }

    #[test]
    fn rigid_rotation_strains_and_stresses_nothing() {
        let x20 = hex20();
        let rot = [[0.0, -0.6, 0.8], [0.6, 0.64, 0.48], [-0.8, 0.48, 0.36]];
        let u: Vec<[f64; 3]> = x20.iter().map(|p| std::array::from_fn(|i| (0..3).map(|jj| rot[i][jj] * p[jj]).sum::<f64>() + 5.0 - p[i])).collect();
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 3000.0 };
        for f in [Formulation::Full, Formulation::Hybrid { pressures: 4 }] {
            let solid = Solid::new(ElementKind::Hex20, &x20, f).unwrap();
            let s = solid.nonlinear_stresses(&neo, &u, Kinematics::Finite, None, &[]).unwrap();
            assert!(s.iter().all(|p| p.stress.iter().chain(&p.strain).all(|v| v.abs() < 1e-12 * 3000.0)), "{f:?}: {:?}", s[0]);
            let r = solid.finite_strain(&neo, &u, &[], &mut [], false).unwrap();
            assert!(r.force.iter().all(|v| v.abs() < 1e-11 * 3000.0));
        }
        // A rotation that inverts the element is refused.
        let flip: Vec<[f64; 3]> = x20.iter().map(|p| [-2.0 * p[0], 0.0, 0.0]).collect();
        let solid = Solid::new(ElementKind::Hex20, &x20, Formulation::Full).unwrap();
        assert!(matches!(solid.finite_strain(&neo, &flip, &[], &mut [], false), Err(Error::Degenerate(_))));
    }

    #[test]
    fn plane_strain_and_axisymmetric_tangents_are_the_derivative_of_the_forces() {
        let quad8 = [
            [1.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [2.2, 1.0, 0.0],
            [0.9, 1.1, 0.0],
            [1.5, 0.0, 0.0],
            [2.1, 0.55, 0.0],
            [1.5, 1.1, 0.0],
            [0.95, 0.5, 0.0],
        ];
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 3000.0 };
        let svk = MaterialModel::StVenantKirchhoff(Box::new(isotropic(1000.0, 0.3).unwrap()));
        for theory in [Theory::PlaneStrain { thickness: 2.0 }, Theory::Axisymmetric] {
            for (kind, x) in [(ElementKind::Quad8, &quad8[..]), (ElementKind::Quad4, &quad8[..4]), (ElementKind::Tri3, &quad8[..3])] {
                let plane = Plane::new(kind, x, Formulation::Full, theory).unwrap();
                for m in [&neo, &svk] {
                    let u: Vec<[f64; 2]> = displacement(x.len(), 0.06).iter().map(|v| [v[0], v[1]]).collect();
                    let k = plane.finite_strain(m, &u, &[], &mut [], true).unwrap().tangent.unwrap();
                    let scale = k.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                    let h = 1e-6;
                    for j in 0..2 * x.len() {
                        let shifted = |s: f64| {
                            let mut v = u.clone();
                            v[j / 2][j % 2] += s;
                            plane.finite_strain(m, &v, &[], &mut [], false).unwrap().force
                        };
                        let (p, q) = (shifted(h), shifted(-h));
                        for i in 0..2 * x.len() {
                            let numeric = (p[i] - q[i]) / (2.0 * h);
                            assert!((numeric - k.get(i, j)).abs() < 1e-6 * scale, "{theory:?} {kind:?} {i} {j}: {numeric} vs {}", k.get(i, j));
                        }
                    }
                }
                // At rest the tangent is the linear stiffness.
                let d = isotropic(1000.0, 0.3).unwrap();
                let k = plane.finite_strain(&svk, &vec![[0.0; 2]; x.len()], &[], &mut [], true).unwrap().tangent.unwrap();
                for (a, b) in k.data.iter().zip(&plane.stiffness(&d).unwrap().data) {
                    assert!((a - b).abs() < 1e-9 * 1000.0);
                }
            }
        }
        let stress = Plane::new(ElementKind::Quad4, &quad8[..4], Formulation::Full, Theory::PlaneStress { thickness: 1.0 }).unwrap();
        assert!(matches!(stress.finite_strain(&svk, &[[0.0; 2]; 4], &[], &mut [], true), Err(Error::NotSupported(_))));
        let stress = Plane::new(ElementKind::Quad4, &quad8[..4], Formulation::IncompatibleModes, Theory::PlaneStress { thickness: 1.0 }).unwrap();
        assert!(matches!(stress.finite_strain(&svk, &[[0.0; 2]; 4], &[0.0; 4], &mut [0.0; 4], true), Err(Error::NotSupported(_))));
    }

    /// Incompatible-mode quad4, plane strain and axisymmetric: the condensed tangent is the
    /// derivative of the forces (the modes re-solved, from rest and from stored modes), elastic
    /// and plastic, and at rest it is the linear element's stiffness.
    #[test]
    fn incompatible_quad4_tangents_are_the_derivative_of_the_forces() {
        use crate::elastoplastic::HISTORY;
        use fieldcraft_model::{Elastic, Hardening, Plastic};
        let quad4 = [[1.0, 0.0, 0.0], [2.0, 0.1, 0.0], [2.2, 1.0, 0.0], [0.9, 1.1, 0.0]];
        let j2 = crate::plasticity::J2::from_model(
            &Elastic::Isotropic { young: 1000.0, poisson: 0.3 },
            Some(&Plastic::Bilinear { yield_stress: 5.0, hardening_modulus: 50.0, hardening: Hardening::Isotropic }),
            None,
        )
        .unwrap();
        let plastic = MaterialModel::Plastic(std::sync::Arc::new(j2));
        let neo = MaterialModel::NeoHookean { shear: 300.0, bulk: 3000.0 };
        let svk = MaterialModel::StVenantKirchhoff(Box::new(isotropic(1000.0, 0.3).unwrap()));
        let u = |scale: f64| -> Vec<[f64; 2]> { displacement(4, scale).iter().map(|v| [v[0], v[1]]).collect() };
        for theory in [Theory::PlaneStrain { thickness: 2.0 }, Theory::Axisymmetric] {
            let plane = Plane::new(ElementKind::Quad4, &quad4, Formulation::IncompatibleModes, theory).unwrap();
            for (m, kin) in [(&neo, Kinematics::Finite), (&svk, Kinematics::Finite), (&plastic, Kinematics::Finite), (&plastic, Kinematics::Small)] {
                let len = plane.nonlinear_history_len(m);
                let mut committed = vec![0.0; len];
                plane.nonlinear_points(m, &u(0.03), kin, None, &vec![0.0; len], &mut committed, false).unwrap();
                assert!(committed[len - 4..].iter().any(|a| *a != 0.0), "the modes are stored");
                let at = u(0.06);
                let k = plane.nonlinear_points(m, &at, kin, None, &committed, &mut committed.clone(), true).unwrap().0.tangent.unwrap();
                let scale = k.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                let h = if kin == Kinematics::Small { 1e-8 } else { 1e-6 };
                for j in 0..8 {
                    let shifted = |s: f64| {
                        let mut v = at.clone();
                        v[j / 2][j % 2] += s;
                        plane.nonlinear_points(m, &v, kin, None, &committed, &mut committed.clone(), false).unwrap().0.force
                    };
                    let (p, q) = (shifted(h), shifted(-h));
                    for i in 0..8 {
                        let numeric = (p[i] - q[i]) / (2.0 * h);
                        assert!(
                            (numeric - k.get(i, j)).abs() < 1e-6 * scale,
                            "{theory:?} {} {kin:?} {i} {j}: {numeric} vs {}",
                            m.name(),
                            k.get(i, j)
                        );
                    }
                }
                assert_eq!(len, if m.history_len() > 0 { 4 * HISTORY + 4 } else { 4 });
            }
            let d = isotropic(1000.0, 0.3).unwrap();
            let k = plane.finite_strain(&svk, &[[0.0; 2]; 4], &[0.0; 4], &mut [0.0; 4], true).unwrap().tangent.unwrap();
            for (a, b) in k.data.iter().zip(&plane.stiffness(&d).unwrap().data) {
                assert!((a - b).abs() < 1e-9 * 1000.0, "{theory:?}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn hybrid_2d_tangents_are_the_derivative_of_the_forces() {
        let quad8 = [
            [1.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [2.2, 1.0, 0.0],
            [0.9, 1.1, 0.0],
            [1.5, 0.0, 0.0],
            [2.1, 0.55, 0.0],
            [1.5, 1.1, 0.0],
            [0.95, 0.5, 0.0],
        ];
        let (shear, bulk) = (300.0, 30000.0);
        let neo = MaterialModel::NeoHookean { shear, bulk };
        let mr = MaterialModel::Polynomial { name: "Mooney–Rivlin", c: [100.0, 50.0, 0.0, 0.0], bulk };
        let hybrid = Formulation::Hybrid { pressures: 1 };
        for theory in [Theory::PlaneStrain { thickness: 2.0 }, Theory::Axisymmetric] {
            for (kind, x) in [(ElementKind::Quad8, &quad8[..]), (ElementKind::Quad4, &quad8[..4])] {
                let plane = Plane::new(kind, x, hybrid, theory).unwrap();
                for m in [&neo, &mr] {
                    let u: Vec<[f64; 2]> = displacement(x.len(), 0.06).iter().map(|v| [v[0], v[1]]).collect();
                    let k = plane.finite_strain(m, &u, &[], &mut [], true).unwrap().tangent.unwrap();
                    let scale = k.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                    let h = 1e-6;
                    for j in 0..2 * x.len() {
                        let shifted = |s: f64| {
                            let mut v = u.clone();
                            v[j / 2][j % 2] += s;
                            plane.finite_strain(m, &v, &[], &mut [], false).unwrap().force
                        };
                        let (p, q) = (shifted(h), shifted(-h));
                        for i in 0..2 * x.len() {
                            let numeric = (p[i] - q[i]) / (2.0 * h);
                            assert!((numeric - k.get(i, j)).abs() < 1e-6 * scale, "{theory:?} {kind:?} {i} {j}: {numeric} vs {}", k.get(i, j));
                        }
                    }
                }
                // At rest the tangent is the linear hybrid stiffness of the same moduli.
                let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
                let poisson = (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear));
                let d = isotropic(young, poisson).unwrap();
                let k = plane.finite_strain(&neo, &vec![[0.0; 2]; x.len()], &[], &mut [], true).unwrap().tangent.unwrap();
                let linear = plane.stiffness(&d).unwrap();
                let scale = linear.data.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                for (a, b) in k.data.iter().zip(&linear.data) {
                    assert!((a - b).abs() < 1e-9 * scale, "{a} {b}");
                }
            }
        }
        // The linear triangle and plane stress have none; St Venant–Kirchhoff has no split.
        assert!(Plane::new(ElementKind::Tri3, &quad8[..3], hybrid, Theory::PlaneStrain { thickness: 1.0 }).is_err());
        assert!(Plane::new(ElementKind::Quad4, &quad8[..4], hybrid, Theory::PlaneStress { thickness: 1.0 }).is_err());
        let svk = MaterialModel::StVenantKirchhoff(Box::new(isotropic(1000.0, 0.3).unwrap()));
        let plane = Plane::new(ElementKind::Quad4, &quad8[..4], hybrid, Theory::Axisymmetric).unwrap();
        assert!(matches!(plane.finite_strain(&svk, &[[0.0; 2]; 4], &[], &mut [], true), Err(Error::NotSupported(_))));
    }

    /// A nearly incompressible neo-Hookean cube in homogeneous compression, free laterally: the
    /// incompatible-mode element keeps its hourglass stiffness to a stretch of 0.6 and stops with
    /// an error below 0.58, where the full-integration element at the same state has lost it
    /// too.
    #[test]
    fn hourglass_instability_in_compression_is_an_error() {
        let x = unit_cube();
        let im = Solid::new(ElementKind::Hex8, &x, Formulation::IncompatibleModes).unwrap();
        let full = Solid::new(ElementKind::Hex8, &x, Formulation::Full).unwrap();
        for bulk in [10.0, 1000.0] {
            let m = MaterialModel::NeoHookean { shear: 1.0, bulk };
            let at = |stretch: f64| -> Vec<[f64; 3]> {
                let lateral = 1.0 / stretch.sqrt();
                x.iter().map(|p| [(lateral - 1.0) * p[0], (lateral - 1.0) * p[1], (stretch - 1.0) * p[2]]).collect()
            };
            im.finite_strain(&m, &at(0.6), &[0.0; 9], &mut [0.0; 9], true).unwrap();
            let err = im.finite_strain(&m, &at(0.57), &[0.0; 9], &mut [0.0; 9], true).unwrap_err();
            assert!(matches!(&err, Error::Degenerate(s) if s.contains("hourglass")), "{err}");
            // The forces alone (no tangent) are not checked: stresses at a converged state.
            im.finite_strain(&m, &at(0.57), &[0.0; 9], &mut [0.0; 9], false).unwrap();
            let k = full.finite_strain(&m, &at(0.57), &[], &mut [], true).unwrap().tangent.unwrap();
            assert!(!positive_on(&k, &hex8_hourglass(&x, &at(0.57)).unwrap()), "full integration has lost it too");
        }
    }

    /// On a distorted hex8 the parent patterns carry a linear part (a mean gradient), the
    /// hourglass modes do not; in homogeneous compression the condensed tangent loses its
    /// stiffness on the modes at a stretch of 0.616, while on the parent patterns it stays
    /// positive to 0.612. Between the two, the check on the modes stops the element and one on
    /// the parent patterns would have let it through.
    #[test]
    fn hourglass_check_on_a_distorted_element_uses_the_true_modes() {
        let cube = unit_cube();
        // On the cube the modes are the parent patterns.
        let parents = hex8_hourglass(&cube, &[[0.0; 3]; 8]).unwrap();
        let mut x = cube.clone();
        x[1] = [1.2, -0.1, 0.1];
        x[6] = [1.4, 1.3, 1.2];
        let uniform = crate::uniform::Hex8::hex8_deformed(&std::array::from_fn(|i| x[i])).unwrap();
        let mean_gradient = |u: &[f64]| {
            uniform
                .gradient(&std::array::from_fn(|n| [u[3 * n], u[3 * n + 1], u[3 * n + 2]]))
                .as_flattened()
                .iter()
                .fold(0.0f64, |a, v| a.max(v.abs()))
        };
        let modes = hex8_hourglass(&x, &[[0.0; 3]; 8]).unwrap();
        assert!(modes.iter().all(|v| mean_gradient(v) < 1e-12), "no mean gradient");
        assert!(modes.iter().all(|v| v.iter().sum::<f64>().abs() < 1e-12), "no mean displacement");
        assert!(parents.iter().any(|v| mean_gradient(v) > 0.05), "the parent patterns have a linear part");

        let im = Solid::new(ElementKind::Hex8, &x, Formulation::IncompatibleModes).unwrap();
        let m = MaterialModel::NeoHookean { shear: 1.0, bulk: 10.0 };
        let at = |stretch: f64| -> Vec<[f64; 3]> {
            let lateral = 1.0 / stretch.sqrt();
            x.iter().map(|p| [(lateral - 1.0) * p[0], (lateral - 1.0) * p[1], (stretch - 1.0) * p[2]]).collect()
        };
        // K* without the check.
        let tangent = |u: &[[f64; 3]]| -> Matrix {
            let none = || Ok((Vec::new(), String::new()));
            let (finite, _) = condensed(&m, u, &x, Kinematics::Finite, 9, &[0.0; 9], &mut [0.0; 9], true, none, |u, material, updated| {
                im.evaluate_points(&m, u, Kinematics::Finite, None, material, updated, true)
            })
            .unwrap();
            finite.tangent.unwrap()
        };
        im.finite_strain(&m, &at(0.62), &[0.0; 9], &mut [0.0; 9], true).unwrap();
        let u = at(0.614);
        let k = tangent(&u);
        assert!(positive_on(&k, &parents), "the parent patterns still see stiffness");
        assert!(!positive_on(&k, &hex8_hourglass(&x, &u).unwrap()), "the hourglass modes do not");
        let err = im.finite_strain(&m, &u, &[0.0; 9], &mut [0.0; 9], true).unwrap_err();
        assert!(matches!(&err, Error::Degenerate(s) if s.contains("hourglass")), "{err}");
    }

    /// Incompatible-mode quad4, plane strain and axisymmetric, of a neo-Hookean material (bulk
    /// modulus ten times the shear modulus) in homogeneous compression, free laterally: it keeps its hourglass
    /// stiffness at a stretch of 0.6 and stops with an error at 0.5, a distorted element as well
    /// as a square one. The fully integrated quad4 keeps its stiffness on the same modes there:
    /// in 2D the check catches the enhanced element's own instability.
    #[test]
    fn hourglass_instability_of_quad4_in_compression_is_an_error() {
        let square = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        let distorted = [[1.0, 0.0, 0.0], [2.0, 0.1, 0.0], [2.2, 1.0, 0.0], [0.9, 1.1, 0.0]];
        let m = MaterialModel::NeoHookean { shear: 1.0, bulk: 10.0 };
        for theory in [Theory::PlaneStrain { thickness: 1.0 }, Theory::Axisymmetric] {
            for x in [&square, &distorted] {
                let im = Plane::new(ElementKind::Quad4, x, Formulation::IncompatibleModes, theory).unwrap();
                let full = Plane::new(ElementKind::Quad4, x, Formulation::Full, theory).unwrap();
                let at = |stretch: f64| -> Vec<[f64; 2]> {
                    let lateral = if theory == Theory::Axisymmetric { 1.0 / stretch.sqrt() } else { 1.0 / stretch };
                    x.iter().map(|p| [(lateral - 1.0) * p[0], (stretch - 1.0) * p[1]]).collect()
                };
                im.finite_strain(&m, &at(0.6), &[0.0; 4], &mut [0.0; 4], true).unwrap();
                let err = im.finite_strain(&m, &at(0.5), &[0.0; 4], &mut [0.0; 4], true).unwrap_err();
                assert!(matches!(&err, Error::Degenerate(s) if s.contains("hourglass")), "{theory:?} {x:?}: {err}");
                im.finite_strain(&m, &at(0.5), &[0.0; 4], &mut [0.0; 4], false).unwrap();
                let k = full.finite_strain(&m, &at(0.5), &[], &mut [], true).unwrap().tangent.unwrap();
                assert!(positive_on(&k, &quad4_hourglass(x, &at(0.5)).unwrap()), "{theory:?} {x:?}: full integration keeps it");
            }
        }
    }

    /// An element that barely moves still solves its modes: the tolerance stays above the
    /// round-off of the Newton correction (it used to fail below 1e-6 of the element's size).
    #[test]
    fn incompatible_modes_converge_at_tiny_displacements() {
        let x = unit_cube();
        let im = Solid::new(ElementKind::Hex8, &x, Formulation::IncompatibleModes).unwrap();
        for bulk in [10.0, 1000.0, 1e6] {
            let m = MaterialModel::NeoHookean { shear: 1.0, bulk };
            for e in [6, 9, 12] {
                let scale = 10f64.powi(-e);
                let mut one = vec![[0.0; 3]; 8];
                one[3][1] = scale;
                for u in [displacement(8, scale), one] {
                    im.finite_strain(&m, &u, &[0.0; 9], &mut [0.0; 9], true).unwrap();
                }
            }
        }
    }

    /// uᵀ K_σ u = ∫ σ_xx (∂u_y/∂x)² dV = σ_xx V for u_y = x under a uniform σ_xx, in a distorted
    /// hex20 and a plane strain quad8, and K_σ annihilates rigid translations.
    #[test]
    fn geometric_stiffness_of_a_uniform_stress() {
        let sigma = 2.5;
        let hex: Vec<[f64; 3]> = crate::basis::Basis::Hex20
            .nodes()
            .iter()
            .map(|p| [1.0 + 0.6 * p[0] + 0.05 * p[1] * p[2], 0.8 * p[1], 0.7 * p[2] + 0.03 * p[0]])
            .collect();
        let solid = Solid::new(ElementKind::Hex20, &hex, Formulation::Full).unwrap();
        let points = solid.rule().points().len();
        let k = solid.geometric_stiffness(&vec![[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]; points]).unwrap();
        let u: Vec<f64> = hex.iter().flat_map(|x| [0.0, x[0], 0.0]).collect();
        let energy: f64 = (0..u.len()).map(|i| u[i] * (0..u.len()).map(|j| k.get(i, j) * u[j]).sum::<f64>()).sum();
        assert!((energy - sigma * solid.volume().unwrap()).abs() < 1e-12, "{energy}");
        let t: Vec<f64> = hex.iter().flat_map(|_| [0.3, -0.2, 1.0]).collect();
        assert!(k.mul_vec(&t).iter().all(|v| v.abs() < 1e-12));
        assert!(solid.geometric_stiffness(&[[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]]).is_err());
        let quad: Vec<[f64; 3]> = crate::basis::Basis::Quad8.nodes().iter().map(|p| [1.0 + 0.5 * p[0] + 0.04 * p[1], 0.6 * p[1], 0.0]).collect();
        let plane =
            crate::plane::Plane::new(ElementKind::Quad8, &quad, Formulation::Full, crate::plane::Theory::PlaneStrain { thickness: 0.3 }).unwrap();
        let k = plane.geometric_stiffness(&vec![[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]; plane.rule().points().len()]).unwrap();
        let u: Vec<f64> = quad.iter().flat_map(|x| [0.0, x[0]]).collect();
        let energy: f64 = (0..u.len()).map(|i| u[i] * (0..u.len()).map(|j| k.get(i, j) * u[j]).sum::<f64>()).sum();
        assert!((energy - sigma * plane.volume().unwrap()).abs() < 1e-12, "{energy}");
    }
}
