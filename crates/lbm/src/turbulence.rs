//! Smagorinsky's subgrid eddy viscosity for the flow collision.
//!
//! The eddy viscosity is `ν_t = (C Δx)² |S|` with `|S| = √(2 S:S)` the local strain rate
//! (J. Smagorinsky, *General circulation experiments with the primitive equations*, Mon. Weather
//! Rev. 91(3):99–164, 1963). In the lattice Boltzmann method the strain rate is local: the
//! non-equilibrium momentum flux of a cell is `P ≈ −2 ρ c_s² τ S` (in lattice units, Δx = Δt = 1),
//! so the cell's own relaxation time `τ = τ₀ + 3 ν_t` solves a quadratic with no finite
//! differences (S. Hou, J. Sterling, S. Chen, G. D. Doolen, *A lattice Boltzmann subgrid model
//! for high Reynolds number flows*, Fields Inst. Commun. 6:151–166, 1996):
//!
//! `τ = ½ (τ₀ + √(τ₀² + 18 √2 C² |P| / ρ))`, with `|P| = √(P:P)`.
//!
//! (From `c_s² (τ − ½) = c_s² (τ₀ − ½) + C² √(2 P:P) / (2 ρ c_s² τ)`, with `c_s² = ⅓`.)
//! The eddy viscosity only adds: `τ ≥ τ₀`, and `τ = τ₀` where the flow does not shear.
//!
//! Built on the regularised collision's flux, as in O. Malaspinas, P. Sagaut, *Consistent subgrid
//! scale modelling for lattice Boltzmann methods*, J. Fluid Mech. 700:514–542, 2012. With heat, the
//! turbulent diffusivity α_t = ν_t / Pr_t (the Reynolds analogy; Pr_t 0.85 to 0.9 for air and
//! water) adds to the scalar lattice's relaxation per cell: τ_T = 4 (α + α_t) + ½. No wall
//! damping yet: Smagorinsky over-predicts the eddy viscosity in the first cells off a wall.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.

use crate::lattice::{C, CS2, Q};

/// The non-equilibrium momentum flux of a cell (xx, yy, zz, xy, xz, yz): `Σ c c f − ρ (c_s² I +
/// u u) + ½ (u F + F u)`, the flux Guo's forcing ties to the strain rate (see
/// [`crate::regularized`]).
#[inline(always)]
pub fn flux(fq: &[f64; Q], rho: f64, u: [f64; 3], force: [f64; 3], forced: bool) -> [f64; 6] {
    let mut pi = [0.0f64; 6];
    for q in (1..Q).step_by(2) {
        let c = C[q];
        let (x, y, z) = (f64::from(c[0]), f64::from(c[1]), f64::from(c[2]));
        let f = fq[q] + fq[q + 1];
        pi[0] += f * x * x;
        pi[1] += f * y * y;
        pi[2] += f * z * z;
        pi[3] += f * x * y;
        pi[4] += f * x * z;
        pi[5] += f * y * z;
    }
    let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    let mut p = [0.0f64; 6];
    for (k, &(a, b)) in pairs.iter().enumerate() {
        let iso = if a == b { CS2 } else { 0.0 };
        p[k] = pi[k] - rho * (iso + u[a] * u[b]) + if forced { 0.5 * (u[a] * force[b] + force[a] * u[b]) } else { 0.0 };
    }
    p
}

/// The relaxation rate of a cell with non-equilibrium flux `p` and density `rho`, where the
/// fluid's own rate is `omega0` and Smagorinsky's constant is `c` (lattice units). Never faster
/// than `omega0`, and `omega0` itself where the flux is not a usable number.
#[inline(always)]
pub fn omega(p: &[f64; 6], rho: f64, omega0: f64, c: f64) -> f64 {
    let tau0 = 1.0 / omega0;
    let pp = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + 2.0 * (p[3] * p[3] + p[4] * p[4] + p[5] * p[5]);
    let tau = 0.5 * (tau0 + (tau0 * tau0 + 18.0 * std::f64::consts::SQRT_2 * c * c * pp.sqrt() / rho).sqrt());
    if tau.is_finite() && tau >= tau0 { 1.0 / tau } else { omega0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// In simple shear `S_xy = γ/2`, the flux is `P_xy = −2 ρ c_s² τ S_xy`; the rate returned
    /// gives the viscosity `ν₀ + C² γ` exactly.
    #[test]
    fn simple_shear_gets_the_smagorinsky_viscosity() {
        for (tau0, c, gamma, rho) in [(0.51, 0.1, 1e-3, 1.0), (0.8, 0.17, 5e-2, 1.3), (1.2, 0.05, 1e-4, 0.9)] {
            // τ such that c_s² (τ − ½) = c_s² (τ₀ − ½) + C² γ.
            let tau = tau0 + c * c * gamma / CS2;
            let p = [0.0, 0.0, 0.0, -2.0 * rho * CS2 * tau * gamma / 2.0, 0.0, 0.0];
            let got = 1.0 / omega(&p, rho, 1.0 / tau0, c);
            assert!((got - tau).abs() < 1e-12 * tau, "τ₀ {tau0}: {got} vs {tau}");
        }
    }

    #[test]
    fn no_shear_no_eddy_viscosity_and_hostile_flux_keeps_the_rate() {
        assert_eq!(omega(&[0.0; 6], 1.0, 1.5, 0.17), 1.5);
        assert_eq!(omega(&[f64::NAN; 6], 1.0, 1.5, 0.17), 1.5);
        assert_eq!(omega(&[1.0; 6], 0.0, 1.5, 0.17), 1.5);
        assert!(omega(&[0.0, 0.0, 0.0, 1e-3, 0.0, 0.0], 1.0, 1.5, 0.17) < 1.5);
    }

    /// The flux of a sheared equilibrium-plus-offset population set is what the regularised
    /// collision computes.
    #[test]
    fn the_flux_is_the_regularised_collisions() {
        let mut f = [0.0; Q];
        for (q, v) in f.iter_mut().enumerate() {
            *v = crate::lattice::equilibrium(q, 1.01, [0.02, -0.01, 0.0]) * (1.0 + 0.003 * ((q * 5 % 7) as f64 - 3.0));
        }
        let rho: f64 = f.iter().sum();
        let p = flux(&f, rho, [0.02, -0.01, 0.0], [0.0; 3], false);
        let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let raw: f64 = (0..Q).map(|q| f[q] * crate::lattice::c(q)[a] * crate::lattice::c(q)[b]).sum();
            let iso = if a == b { CS2 } else { 0.0 };
            let want = raw - rho * (iso + [0.02, -0.01, 0.0][a] * [0.02, -0.01, 0.0][b]);
            assert!((p[k] - want).abs() < 1e-15, "{k}");
        }
    }
}
