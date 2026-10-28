//! Vreman's subgrid eddy viscosity (`turbulence: {"type": "vreman"}`).
//!
//! A. W. Vreman, *An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory
//! and applications*, Phys. Fluids 16(10):3670–3681, 2004, doi:10.1063/1.1785131 (§I, eqs. 5–8):
//!
//! `ν_t = c √(B_β / (α_ij α_ij))`, `α_ij = ∂_i u_j`, `β_ij = Σ_m Δ_m² α_mi α_mj`,
//! `B_β = β₁₁β₂₂ − β₁₂² + β₁₁β₃₃ − β₁₃² + β₂₂β₃₃ − β₂₃²`,
//!
//! and `ν_t = 0` where `α_ij α_ij = 0`. The constant is `c ≈ 2.5 C_s²` with `C_s` the Smagorinsky
//! constant (Vreman: `C_s` = 0.17 gives `c` = 0.07; J. Jacob, O. Malaspinas, P. Sagaut, J. Turbul.
//! 19:1051, 2018, §5, use `C_s` = 0.18). The step gives `C_s` (`constant`); the solver uses
//! `c = 2.5 C_s²`.
//!
//! **Where it vanishes.** `β` is positive semidefinite and `B_β` its second invariant, so
//! `B_β ≥ 0`. On the uniform cubic cells of the lattice (`Δ_m = Δx`, 1 in lattice units)
//! `β = Δx² αᵀα`, and by the Cauchy–Binet formula `B_β = Δx⁴ Σ M²` over the nine 2 × 2 minors `M`
//! of `α`: `B_β` is zero exactly when `α` has rank at most one. That covers the thirteen locally
//! laminar flow types of Vreman's §II whose theoretical subgrid dissipation is zero (one nonzero
//! off-diagonal derivative, or two in one row or one column of `α`, eqs. 14–15) and all their
//! rotations, and a no-slip wall in pure shear. It is computed here as that sum of squared minors:
//! the same number as eq. 8, but never negative from round-off, and exactly 0 when `α` has a
//! single nonzero row or column (every minor is then a product with an exact zero). A rigid
//! rotation, `u = Ω (−y, x, 0)`, has rank two: `B_β = Ω⁴ Δx⁴`, `α:α = 2 Ω²`, so
//! `ν_t = c Δx² Ω / √2`, not zero (Vreman's zero set is defined by the theoretical subgrid
//! dissipation, and its "rotations" are rotations of the coordinate axes, §II).
//!
//! (Jacob et al. write `α_ij = ∂u_i/∂x_j`, the transpose; `αᵀα` and `ααᵀ` have the same
//! invariants, so `B_β` is the same.)
//!
//! **Gradients.** `α` comes from the same finite differences of the neighbours' velocity as the
//! hybrid collision's strain rate ([`crate::hrr::gradient`]): second-order centred, first-order
//! one-sided where one neighbour along an axis is not a fluid cell, 0 where neither is. The
//! velocities are those of a pass over every fluid cell before the collision, so the eddy
//! viscosity is **explicit**: from this step's neighbours, unlike Hou's implicit Smagorinsky
//! ([`crate::turbulence`]), which solves for the rate from the cell's own non-equilibrium flux.
//! The cell then relaxes with `τ̄ = 3 (ν + ν_t) + ½` (every collision; for two relaxation times
//! τ⁻ follows keeping Λ, for multiple relaxation times the five stresses take it). A wall cell
//! under the wall model keeps the wall model's rate.
//!
//! No code of any lattice Boltzmann or large-eddy simulation package was read; this is written
//! from the papers.

/// The Smagorinsky-equivalent constant `C_s` the step gives, as Vreman's `c = 2.5 C_s²`.
pub fn model_constant(smagorinsky_equivalent: f64) -> f64 {
    2.5 * smagorinsky_equivalent * smagorinsky_equivalent
}

/// Vreman's eddy viscosity (lattice units, `Δx = 1`) of a cell whose velocity gradient is `g`
/// (`g[i][j] = ∂_i u_j`), with model constant `c` (`2.5 C_s²`). 0 where the gradient is zero or
/// not a usable number.
#[inline(always)]
pub fn eddy_viscosity(g: &[[f64; 3]; 3], c: f64) -> f64 {
    let aa = crate::hrr::norm2(g);
    if !(aa > 0.0 && aa.is_finite()) {
        return 0.0;
    }
    let b = b_beta(g);
    let nu = c * (b / aa).sqrt();
    if nu.is_finite() && nu > 0.0 { nu } else { 0.0 }
}

/// `B_β` for `β = αᵀα` (`Δ = 1`): the sum of the squared 2 × 2 minors of `α` (Cauchy–Binet).
#[inline(always)]
fn b_beta(g: &[[f64; 3]; 3]) -> f64 {
    // Rows i < j and columns k < l of α, (0, 1), (0, 2), (1, 2) each.
    let pairs = |r: &[f64; 3]| [(r[0], r[1]), (r[0], r[2]), (r[1], r[2])];
    let [r0, r1, r2] = g;
    let mut b = 0.0;
    for (ri, rj) in [(r0, r1), (r0, r2), (r1, r2)] {
        for ((ik, il), (jk, jl)) in pairs(ri).into_iter().zip(pairs(rj)) {
            let m = ik * jl - il * jk;
            b += m * m;
        }
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Eq. 8 as Vreman writes it, with `β_ij = Σ_m α_mi α_mj`.
    fn b_beta_eq8(g: &[[f64; 3]; 3]) -> f64 {
        let beta = |i: usize, j: usize| (0..3).map(|m| g[m][i] * g[m][j]).sum::<f64>();
        beta(0, 0) * beta(1, 1) - beta(0, 1).powi(2) + beta(0, 0) * beta(2, 2) - beta(0, 2).powi(2) + beta(1, 1) * beta(2, 2) - beta(1, 2).powi(2)
    }

    #[test]
    fn the_minors_are_eq_8() {
        let mut x = 0.123_f64;
        for _ in 0..200 {
            let g = [0, 1, 2].map(|_| {
                [0, 1, 2].map(|_| {
                    x = (x * 9301.0 + 0.4927).fract();
                    x - 0.5
                })
            });
            let (a, b) = (b_beta(&g), b_beta_eq8(&g));
            assert!((a - b).abs() <= 1e-14 * b.abs().max(1e-3), "{a} vs {b}");
            assert!(a >= 0.0);
        }
    }

    #[test]
    fn zero_and_rank_one_gradients_give_zero_exactly() {
        assert_eq!(eddy_viscosity(&[[0.0; 3]; 3], 0.081), 0.0);
        // Pure shear ∂_y u_x = γ, each single component, and the two-component types of Vreman's
        // eqs. 14 and 15 (one row, one column of α).
        let gamma = 0.0137;
        for i in 0..3 {
            for j in 0..3 {
                let mut g = [[0.0; 3]; 3];
                g[i][j] = gamma;
                assert_eq!(eddy_viscosity(&g, 0.081), 0.0, "{i}{j}");
            }
        }
        let row = [[0.0, 0.3, -0.7], [0.0; 3], [0.0; 3]];
        let column = [[0.0, 0.0, 0.3], [0.0, 0.0, -0.7], [0.0; 3]];
        assert_eq!(eddy_viscosity(&row, 0.081), 0.0);
        assert_eq!(eddy_viscosity(&column, 0.081), 0.0);
        // A general rank-one gradient a bᵀ: zero up to round-off.
        let (a, b) = ([0.3, -1.1, 0.7], [0.2, 0.5, -0.4]);
        let g = a.map(|ai| b.map(|bj| ai * bj));
        assert!(eddy_viscosity(&g, 0.081) < 1e-9, "{}", eddy_viscosity(&g, 0.081));
        // Not a number: no eddy viscosity.
        assert_eq!(eddy_viscosity(&[[f64::NAN; 3]; 3], 0.081), 0.0);
        assert_eq!(eddy_viscosity(&[[f64::INFINITY; 3]; 3], 0.081), 0.0);
    }

    /// A rigid rotation u = Ω (−y, x, 0): α_yx = −Ω, α_xy = Ω; β₁₁ = β₂₂ = Ω², B_β = Ω⁴,
    /// α:α = 2Ω², so ν_t = c Ω / √2 (Δ = 1). Not zero.
    #[test]
    fn a_rigid_rotation_has_c_omega_over_root_two() {
        for omega in [1e-3, 0.02, -0.3] {
            let g = [[0.0, omega, 0.0], [-omega, 0.0, 0.0], [0.0; 3]];
            let c = model_constant(0.18);
            let want = c * f64::abs(omega) / 2f64.sqrt();
            assert!((eddy_viscosity(&g, c) - want).abs() <= 1e-15 * want, "{omega}");
        }
    }

    /// Axisymmetric strain u = s (x, y, −2z) (rank three): β = s² diag(1, 1, 4),
    /// B_β = s⁴ (1 + 4 + 4) = 9 s⁴, α:α = 6 s², so ν_t = c √(3/2) |s|; and plane strain
    /// u = s (x, −y, 0): B_β = s⁴, α:α = 2 s², ν_t = c |s| / √2, as the rotation.
    #[test]
    fn strains_have_their_closed_forms() {
        let (c, s) = (model_constant(0.17), 0.04);
        let axisymmetric = [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, -2.0 * s]];
        let want = c * 1.5f64.sqrt() * s;
        assert!((eddy_viscosity(&axisymmetric, c) - want).abs() <= 1e-15 * want);
        let plane = [[s, 0.0, 0.0], [0.0, -s, 0.0], [0.0; 3]];
        let want = c * s / 2f64.sqrt();
        assert!((eddy_viscosity(&plane, c) - want).abs() <= 1e-15 * want);
        assert!((model_constant(0.17) - 0.07225).abs() < 1e-15);
    }
}
