//! The hybrid recursive regularised collision (`collision: "hrr"`) on D3Q19.
//!
//! # Scheme
//! J. Jacob, O. Malaspinas, P. Sagaut, *A new hybrid recursive regularised Bhatnagar–Gross–Krook
//! collision model for Lattice Boltzmann method-based large eddy simulation*, J. Turbul.
//! 19:1051–1076, 2018, doi:10.1080/14685248.2018.1540879. It is the recursive regularised
//! collision of [`crate::regularized`] (their RR-BGK: the same equilibrium, the same third-order
//! Hermite pairs, the same Guo force term and relaxation) with one change: the second-order
//! non-equilibrium flux the populations are rebuilt from is a blend (their eq. 31, 37)
//!
//! `P̃ = σ P + (1 − σ)(−2 ρ c_s² τ̄ S^FD)`, `0 < σ ≤ 1`,
//!
//! of the populations' own flux `P` (with Guo's `½ (u F + F u)`) and the flux the strain rate of
//! finite differences of the velocity gives. The third-order non-equilibrium is rebuilt
//! recursively from `P̃` (`u_a P̃_bc + u_b P̃_ac + u_c P̃_ab`). σ = 1 is the recursive regularised
//! collision, bit for bit; σ is 0.99 unless the step gives another (`sigma`).
//!
//! **τ̄.** The paper writes the blend with the continuous relaxation time τ and discretises with
//! `τ̄ = τ + ½` (their eq. 46). Here `τ̄ = 1/ω` is the cell's discrete relaxation time, the one its
//! collision actually uses: with Smagorinsky's eddy viscosity it includes it, and a wall cell
//! under the wall model takes the wall model's rate. The reason is the discrete flux: the
//! Chapman–Enskog analysis of the discrete scheme ties the populations' flux before the
//! collision to the strain rate by `P ≈ −2 ρ c_s² τ̄ S` (with τ̄, not τ̄ − ½), so with τ̄ the
//! finite-difference flux estimates the same quantity, and the scheme keeps the viscosity
//! `ν = c_s² (τ̄ − ½)` exactly as the wavenumber goes to zero (the linear analysis of a shear
//! wave in `tests/hrr.rs`, checked against the measured decay).
//!
//! **Stencil.** `S^FD_ab = ½ (∂_a u_b + ∂_b u_a)` with second-order centred differences of the
//! neighbouring cells' velocity `u` (with Guo's half force, the velocity the collision uses):
//! `∂_y u_x (i, j, k) = (u_x(i, j+1, k) − u_x(i, j−1, k)) / 2` in lattice units. The neighbours
//! are those streaming uses: across a periodic face the periodic image. Y. Feng, P. Boivin,
//! J. Jacob, P. Sagaut, *Hybrid recursive regularized thermal lattice Boltzmann model for high
//! subsonic compressible flows*, J. Comput. Phys. 394:82–99, 2019, doi:10.1016/j.jcp.2019.05.031,
//! use an isotropic stencil weighted by the lattice weights over all the directions instead; that
//! variant is not implemented.
//!
//! **Walls and missing neighbours.** Where one neighbour along an axis is not a fluid cell (an
//! obstacle cell, behind a wall face or outside an inlet or outlet face of the box), that
//! derivative is the first-order one-sided difference with the other neighbour,
//! `u(i+1) − u(i)` or `u(i) − u(i−1)`; with neither, it is 0. This is our choice for every kind
//! of boundary. Jacob et al. use one-sided first-order differences where neighbours are missing
//! at their grid-refinement interfaces (§3) and state no wall treatment.
//!
//! **Hyperviscosity.** The centred difference's leading error, `∂ g + Δx²/6 ∂³ g`, makes the
//! blend add a bi-Laplacian dissipation; at leading order (their eq. 47)
//! `ν_σ = (1 − σ)/(6σ) Δx² c_s² τ`, so a shear wave of wavenumber k decays at
//! `ν k² [1 + ((1 − σ)/σ)(1 − sin(k Δx)/(k Δx))]` for small k Δx: the shortest waves, which an
//! under-resolved large-eddy simulation piles energy into, are damped, the resolved ones barely.
//!
//! The velocities of the neighbours come from a pass over every fluid cell before the collision
//! ([`crate::flow`]), by fixed blocks, so the result is the same on any number of threads.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.

/// The strain rate of finite differences at a cell of velocity `u`: `neighbour(axis, up)` is
/// the velocity of the fluid cell next to it along `axis` (towards `+axis` with `up`), `None`
/// when that cell is not a fluid cell. Centred differences, one-sided where one neighbour is
/// missing, 0 where both are. Returns `xx, yy, zz, xy, xz, yz` (lattice units).
#[inline(always)]
pub fn strain(u: [f64; 3], neighbour: impl Fn(usize, bool) -> Option<[f64; 3]>) -> [f64; 6] {
    // g[a][b] = ∂_a u_b
    let g: [[f64; 3]; 3] = [0, 1, 2].map(|a| match (neighbour(a, false), neighbour(a, true)) {
        (Some(m), Some(p)) => [0, 1, 2].map(|b| 0.5 * (p[b] - m[b])),
        (Some(m), None) => [0, 1, 2].map(|b| u[b] - m[b]),
        (None, Some(p)) => [0, 1, 2].map(|b| p[b] - u[b]),
        (None, None) => [0.0; 3],
    });
    [g[0][0], g[1][1], g[2][2], 0.5 * (g[0][1] + g[1][0]), 0.5 * (g[0][2] + g[2][0]), 0.5 * (g[1][2] + g[2][1])]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::{CS2, Q, c};

    /// The hybrid collision of one cell, as the flow calls it.
    #[allow(clippy::too_many_arguments)]
    fn collide(
        f: &[f64; Q],
        rho: f64,
        u: [f64; 3],
        force: [f64; 3],
        omega: f64,
        smagorinsky: f64,
        forced: bool,
        sigma: f64,
        strain: &[f64; 6],
        out: &mut [f64; Q],
    ) -> f64 {
        crate::regularized::rebuild(f, rho, u, force, omega, smagorinsky, forced, Some((sigma, strain)), out)
    }

    fn state() -> ([f64; Q], f64, [f64; 3], [f64; 3]) {
        let mut f = [0.0; Q];
        for (q, x) in f.iter_mut().enumerate() {
            *x = crate::lattice::equilibrium(q, 1.01, [0.02, -0.01, 0.03]) * (1.0 + 0.02 * ((q * 5 % 7) as f64 - 3.0));
        }
        let force = [2e-4, 1e-4, -3e-5];
        let (mut rho, mut j) = (0.0, [0.0; 3]);
        for (q, x) in f.iter().enumerate() {
            rho += x;
            for k in 0..3 {
                j[k] += x * c(q)[k];
            }
        }
        (f, rho, [0, 1, 2].map(|k| (j[k] + 0.5 * force[k]) / rho), force)
    }

    #[test]
    fn sigma_one_is_the_regularised_collision_bit_for_bit() {
        let (f, rho, u, force) = state();
        let s = [1e-3, -2e-3, 5e-4, 3e-3, -1e-3, 2e-3];
        for smagorinsky in [0.0, 0.17] {
            let (mut a, mut b) = ([0.0; Q], [0.0; Q]);
            let wa = crate::regularized::collide(&f, rho, u, force, 1.7, smagorinsky, true, &mut a);
            let wb = collide(&f, rho, u, force, 1.7, smagorinsky, true, 1.0, &s, &mut b);
            assert_eq!(wa.to_bits(), wb.to_bits());
            for q in 0..Q {
                assert_eq!(a[q].to_bits(), b[q].to_bits(), "{q}");
            }
        }
    }

    #[test]
    fn the_rebuilt_flux_is_the_blend() {
        let (f, rho, u, force) = state();
        let s = [1e-3, -2e-3, 5e-4, 3e-3, -1e-3, 2e-3];
        let (omega, sigma) = (1.6, 0.7);
        let mut out = [0.0; Q];
        collide(&f, rho, u, force, omega, 0.0, true, sigma, &s, &mut out);
        // The post-collision flux is ρ(c_s² I + u u) + (1 − ω) P̃ + ½ (u F + F u).
        let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let (mut pi, mut pi_out) = (0.0, 0.0);
            for q in 0..Q {
                pi += f[q] * c(q)[a] * c(q)[b];
                pi_out += out[q] * c(q)[a] * c(q)[b];
            }
            let iso = if a == b { CS2 } else { 0.0 };
            let uf = 0.5 * (u[a] * force[b] + force[a] * u[b]);
            let p = pi - rho * (iso + u[a] * u[b]) + uf;
            let blend = sigma * p + (1.0 - sigma) * (-2.0 * CS2 * rho / omega) * s[k];
            let want = rho * (iso + u[a] * u[b]) + (1.0 - omega) * blend + uf;
            assert!((pi_out - want).abs() < 1e-15, "{k}: {pi_out} vs {want}");
        }
        // Mass and momentum as before.
        let (r, j): (f64, [f64; 3]) = (out.iter().sum(), [0, 1, 2].map(|k| (0..Q).map(|q| out[q] * c(q)[k]).sum()));
        assert!((r - rho).abs() < 1e-15);
        for k in 0..3 {
            assert!((j[k] - (rho * u[k] + 0.5 * force[k])).abs() < 1e-15);
        }
    }

    #[test]
    fn the_strain_rate_is_centred_one_sided_or_zero() {
        let u = [1.0, 2.0, 3.0];
        // u_x = y, u_y = 2 x, u_z = 0.5 z (+ constants), every neighbour present.
        let full = |a: usize, up: bool| {
            let d = if up { 1.0 } else { -1.0 };
            Some(match a {
                0 => [u[0], u[1] + 2.0 * d, u[2]],
                1 => [u[0] + d, u[1], u[2]],
                _ => [u[0], u[1], u[2] + 0.5 * d],
            })
        };
        let s = strain(u, full);
        assert_eq!(s, [0.0, 0.0, 0.5, 1.5, 0.0, 0.0]);
        // The +y neighbour missing: the one-sided difference with −y gives the same here; both
        // y neighbours missing: ∂_y = 0.
        let s = strain(u, |a, up| if a == 1 && up { None } else { full(a, up) });
        assert_eq!(s, [0.0, 0.0, 0.5, 1.5, 0.0, 0.0]);
        let s = strain(u, |a, up| if a == 1 { None } else { full(a, up) });
        assert_eq!(s, [0.0, 0.0, 0.5, 1.0, 0.0, 0.0]);
        assert_eq!(strain(u, |_, _| None), [0.0; 6]);
    }
}
