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
//! **Dynamic σ** (`sigma: "dynamic"`, their §4, eqs. 48–50). With a subgrid model's eddy viscosity
//! ν_t, σ is chosen per cell and step so that the blend's hyperviscous dissipation
//! `ε_σ = ν_σ |∇²u|²` equals the subgrid one `ε_sgs = ν_t |∇u|²`: `ν_σ = L_VK² ν_t` with
//! `L_VK = |∇u| / |∇²u|` (a von Kármán length), `|∇u|² = α_ij α_ij` from [`gradient`] and
//! `|∇²u|² = Σ_j (∇²u_j)²` from [`laplacian`]. With `ν_σ = (1 − σ)/(6σ) Δx² c_s² τ` and the paper's
//! continuous τ read as `c_s² τ = ν`, the molecular viscosity (the hyperviscosity the scheme's
//! linear analysis shows, `ν (1 − σ)/σ (1 − sin(kΔx)/(kΔx)) ≈ ν (1 − σ)/(6σ) (kΔx)²`, carries
//! ν, not ν + ν_t), so
//!
//! `σ = 1 / (1 + 6 L_VK² ν_t / (Δx² ν))`,
//!
//! computed as `|∇²u|² / (|∇²u|² + 6 |∇u|² ν_t / ν)` (Δx = 1), which lies in [0, 1] by
//! construction: ν_t = 0 gives σ = 1 exactly (the regularised collision, bit for bit);
//! `|∇²u| = 0` with ν_t > 0 gives σ = 0 (the finite-difference flux alone, which the blend
//! allows; only a fixed σ = 0 is refused). The subgrid viscosity is then **not** added to τ (the
//! paper's case 4, implicit LES): the cell relaxes with the molecular τ̄, and the model only sets σ.
//! Vreman's ν_t ([`crate::vreman`]) is the paper's; with Smagorinsky's model ν_t is Hou's, the
//! `(1/ω − 1/ω₀)/3` the implicit rate [`crate::turbulence::omega`] would give the cell, computed
//! without that cancellation and 0 where the flux is below its round-off
//! ([`crate::turbulence::eddy_viscosity`]): σ multiplies ν_t by `L_VK²`, unbounded where ∇²u
//! vanishes, so a rounding ν_t would set σ there (towards 0). A
//! wall cell under the wall model keeps σ = 1 (its rate is the wall model's).
//!
//! **Laplacian.** `∇²u_j` is the 7-point second-order centred Laplacian of each component from the
//! same neighbours, `Σ_a (u(i + e_a) − 2 u(i) + u(i − e_a))`; along an axis where either neighbour
//! is missing that axis's second difference is 0 (our choice, as one-sided first differences carry
//! no second derivative).
//!
//! The velocities of the neighbours come from a pass over every fluid cell before the collision
//! ([`crate::flow`]), by fixed blocks, so the result is the same on any number of threads.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.

/// The velocities of a cell's six neighbours, `[axis][up]` (`up` towards `+axis`): `None` where
/// that cell is not a fluid cell (or lies outside a face that is not periodic).
pub type Neighbours = [[Option<[f64; 3]>; 2]; 3];

/// The velocity gradient of finite differences at a cell of velocity `u` with neighbours `nb`:
/// `g[a][b] = ∂_a u_b` (lattice units). Centred differences, first-order one-sided where one
/// neighbour along an axis is missing, 0 where both are. Vreman's eddy viscosity
/// ([`crate::vreman`]) uses the same gradient.
#[inline(always)]
pub fn gradient(u: [f64; 3], nb: &Neighbours) -> [[f64; 3]; 3] {
    nb.map(|[m, p]| match (m, p) {
        (Some(m), Some(p)) => [0, 1, 2].map(|b| 0.5 * (p[b] - m[b])),
        (Some(m), None) => [0, 1, 2].map(|b| u[b] - m[b]),
        (None, Some(p)) => [0, 1, 2].map(|b| p[b] - u[b]),
        (None, None) => [0.0; 3],
    })
}

/// The squared norm `α_ij α_ij` of gradient `g`.
#[inline(always)]
pub fn norm2(g: &[[f64; 3]; 3]) -> f64 {
    g.iter().flatten().map(|x| x * x).sum()
}

/// The Laplacian of finite differences of each velocity component at a cell of velocity `u` with
/// neighbours `nb`: along each axis with both neighbours `u(i+1) − 2 u(i) + u(i−1)`, nothing along
/// an axis where either is missing.
#[inline(always)]
pub fn laplacian(u: [f64; 3], nb: &Neighbours) -> [f64; 3] {
    let mut l = [0.0; 3];
    for [m, p] in nb {
        if let (Some(m), Some(p)) = (m, p) {
            for (((lb, pb), mb), ub) in l.iter_mut().zip(p).zip(m).zip(u) {
                *lb += pb - 2.0 * ub + mb;
            }
        }
    }
    l
}

/// The dynamic σ of a cell (Jacob et al. eqs. 48–50, see the module docs) with gradient `g`,
/// Laplacian `lap`, subgrid eddy viscosity `nu_t` and molecular viscosity `nu` (lattice units):
/// `|∇²u|² / (|∇²u|² + 6 |∇u|² ν_t / ν)`, in [0, 1]; exactly 1 where `nu_t` is not above 0, and 1
/// where the result is not a number (no gradient and no Laplacian).
#[inline(always)]
pub fn dynamic_sigma(g: &[[f64; 3]; 3], lap: &[f64; 3], nu_t: f64, nu: f64) -> f64 {
    if !(nu_t > 0.0 && nu > 0.0) {
        return 1.0;
    }
    let grad2 = norm2(g);
    let lap2: f64 = lap.iter().map(|x| x * x).sum();
    let sigma = lap2 / (lap2 + 6.0 * grad2 * nu_t / nu);
    if sigma.is_finite() { sigma.clamp(0.0, 1.0) } else { 1.0 }
}

/// The strain rate of gradient `g` ([`gradient`]): `xx, yy, zz, xy, xz, yz`.
#[inline(always)]
pub fn strain(g: &[[f64; 3]; 3]) -> [f64; 6] {
    let [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = *g;
    [xx, yy, zz, 0.5 * (xy + yx), 0.5 * (xz + zx), 0.5 * (yz + zy)]
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

    /// The neighbours of `u` given by `at(axis, up)`.
    fn neighbours(at: impl Fn(usize, bool) -> Option<[f64; 3]>) -> Neighbours {
        [0, 1, 2].map(|a| [false, true].map(|up| at(a, up)))
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
        let g = gradient(u, &neighbours(full));
        assert_eq!(g, [[0.0, 2.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.5]]);
        assert_eq!(strain(&g), [0.0, 0.0, 0.5, 1.5, 0.0, 0.0]);
        // The +y neighbour missing: the one-sided difference with −y gives the same here; both
        // y neighbours missing: ∂_y = 0.
        let s = strain(&gradient(u, &neighbours(|a, up| if a == 1 && up { None } else { full(a, up) })));
        assert_eq!(s, [0.0, 0.0, 0.5, 1.5, 0.0, 0.0]);
        let s = strain(&gradient(u, &neighbours(|a, up| if a == 1 { None } else { full(a, up) })));
        assert_eq!(s, [0.0, 0.0, 0.5, 1.0, 0.0, 0.0]);
        assert_eq!(strain(&gradient(u, &neighbours(|_, _| None))), [0.0; 6]);
    }

    /// u_x = sin(k y) on the lattice: the centred differences give α_yx = cos(ky) sin k and
    /// ∇²u_x = 2 sin(ky)(cos k − 1), so L_VK² = cos²(ky) sin²k / (4 sin²(ky)(1 − cos k)²) and
    /// σ = 1 / (1 + 6 L_VK² ν_t / ν) (the continuous limit: L_VK = |cot(ky)| / k).
    #[test]
    fn dynamic_sigma_matches_the_closed_form_of_a_sine() {
        let (k, nu, nu_t) = (2.0 * std::f64::consts::PI / 16.0, 0.01, 0.003);
        for y in [0.3, 1.0, 2.5, 5.0] {
            let at = |y: f64| [(k * y).sin(), 0.0, 0.0];
            let u = at(y);
            let nb = neighbours(|a, up| Some(if a == 1 { at(y + if up { 1.0 } else { -1.0 }) } else { u }));
            let (g, lap) = (gradient(u, &nb), laplacian(u, &nb));
            let (alpha, lx) = ((k * y).cos() * k.sin(), 2.0 * (k * y).sin() * (k.cos() - 1.0));
            assert!((g[1][0] - alpha).abs() < 1e-15 && (lap[0] - lx).abs() < 1e-15, "{g:?} {lap:?}");
            let l2 = alpha * alpha / (lx * lx);
            let want = 1.0 / (1.0 + 6.0 * l2 * nu_t / nu);
            let got = dynamic_sigma(&g, &lap, nu_t, nu);
            assert!((got - want).abs() < 1e-14, "y {y}: {got} vs {want}");
            assert!((0.0..=1.0).contains(&got));
        }
        // The limits: no eddy viscosity is σ = 1 exactly; no Laplacian with a gradient is 0;
        // neither is 1; hostile numbers give 1.
        let g = [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0; 3]];
        assert_eq!(dynamic_sigma(&g, &[0.1, 0.0, 0.0], 0.0, 0.01), 1.0);
        assert_eq!(dynamic_sigma(&g, &[0.0; 3], 0.003, 0.01), 0.0);
        assert_eq!(dynamic_sigma(&[[0.0; 3]; 3], &[0.0; 3], 0.003, 0.01), 1.0);
        assert_eq!(dynamic_sigma(&[[f64::NAN; 3]; 3], &[f64::NAN; 3], 0.003, 0.01), 1.0);
        assert_eq!(dynamic_sigma(&g, &[0.1, 0.0, 0.0], f64::NAN, 0.01), 1.0);
        // The Laplacian of a quadratic is exact; an axis with a missing neighbour adds nothing.
        let quad = |x: f64, y: f64| [x * x + 3.0 * y * y, 0.0, 0.0];
        let nb = neighbours(|a, up| {
            let d = if up { 1.0 } else { -1.0 };
            Some(match a {
                0 => quad(1.0 + d, 2.0),
                1 => quad(1.0, 2.0 + d),
                _ => quad(1.0, 2.0),
            })
        });
        assert_eq!(laplacian(quad(1.0, 2.0), &nb), [8.0, 0.0, 0.0]);
        let mut half = nb;
        half[1][1] = None;
        assert_eq!(laplacian(quad(1.0, 2.0), &half), [2.0, 0.0, 0.0]);
    }
}
