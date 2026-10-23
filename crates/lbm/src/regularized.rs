//! The recursive regularised collision on D3Q19 with Guo forcing.
//!
//! # Scheme
//! Every population is rebuilt from the cell's hydrodynamic moments before it relaxes
//! (J. Latt, B. Chopard, *Lattice Boltzmann method with regularized pre-collision distribution
//! functions*, Math. Comput. Simul. 72(2–6):165–168, 2006, doi:10.1016/j.matcom.2006.05.017):
//! the density ρ, the velocity `u = (Σ f c + F/2) / ρ` of Guo's scheme, and the non-equilibrium
//! momentum flux. What the
//! populations carry beyond those moments, the non-hydrodynamic ("ghost") modes, is dropped
//! instead of being relaxed with `1 − ω`. As τ → ½ that factor approaches −1, so the ghost modes
//! of the two-relaxation-time scheme barely decay and grow at sharp corners and open faces; here
//! they are set to their equilibrium every step.
//!
//! The rebuilt populations are a Hermite expansion. Second order: the equilibrium of
//! [`crate::lattice::equilibrium`] and the non-equilibrium flux `P`. Third order: the moments
//! D3Q19 can hold, `xxy, yzz, xyy, xzz, xxz, yyz` (not `xyz`, and `xxx` is the same as `x`), with
//! the equilibrium `ρ u_a u_b u_c` and, recursively from `P`, the non-equilibrium
//! `u_a P_bc + u_b P_ac + u_c P_ab` (C. Coreixas, G. Wissocq, G. Puigt, J.-F. Boussuge,
//! P. Sagaut, *Recursive regularization step for high-order lattice Boltzmann methods*, Phys.
//! Rev. E 96:033306, 2017, doi:10.1103/PhysRevE.96.033306, with the D3Q19 form). On D3Q19 the
//! third-order Hermite polynomials of a pair such as `H_xxy` and `H_yzz` are not orthogonal under
//! the weights, so the pair enters through its sum and difference, which are: with
//! `Σ w (H_xxy + H_yzz)² = 2/27` and `Σ w (H_xxy − H_yzz)² = 6/27`, the term is
//! `w_q [27/2 (H_xxy + H_yzz)(a_xxy + a_yzz) + 9/2 (H_xxy − H_yzz)(a_xxy − a_yzz)]`, which gives
//! the populations exactly the moments `a_xxy` and `a_yzz`.
//!
//! # Body force
//! Guo, Zheng and Shi's update (Phys. Rev. E 65:046308, 2002, doi:10.1103/PhysRevE.65.046308) is
//! `f* = f^eq + (1 − ω) f^neq + (1 − ω/2) F_q` with `F_q = w_q [3 c·F + 9 (c·u)(c·F) − 3 u·F]`,
//! whose first and second Hermite moments are `F` and `u F + F u`. With `u` defined with the half
//! force, the populations' own first moment off equilibrium is `Σ c f^neq = −F/2`, and their
//! non-equilibrium flux is `Π^neq = Σ c c f − ρ c_s² I − ρ u u`. Regularising `f^neq` keeps its
//! first- and second-order Hermite parts; collecting terms,
//!
//! `f* = f^eq + (1 − ω) f^(1)[P] + ½ F_q`, with `P = Π^neq + ½ (u F + F u)`
//!
//! (the first-order parts give `w_q c·F (3/2)` for any ω, and `(1 − ω) Π^neq + (1 − ω/2)(uF + Fu)
//! = (1 − ω) P + ½ (uF + Fu)`). `P` is the flux the Chapman–Enskog analysis of Guo's scheme ties
//! to the strain rate, `P ≈ −2 ρ c_s² τ S`, so it is also what the third-order recursion takes. The
//! update conserves mass, gives the momentum `ρ u + F/2` (so the next step's velocity is right),
//! and has the viscosity `ν = c_s² (τ − ½)` of the other schemes; the fans', porous zones' and
//! buoyancy's forces enter through `F` and `u` exactly as in the two-relaxation-time collision.
//!
//! # Why this operator
//! On D3Q19 the cumulant collision (M. Geier, M. Schönherr, A. Pasquali, M. Krafczyk, Comput.
//! Math. Appl. 70:507–547, 2015, doi:10.1016/j.camwa.2015.05.001) is not available in its proper
//! form: it is built on the 27 moments of D3Q27, and D3Q19 lacks the corner directions its higher
//! cumulants need. The projected regularisation (second order only) also filters the ghost
//! modes, but sets the third-order moments to zero, which loses Galilean invariance of the
//! viscous stress as the velocity grows. The recursive form keeps the third-order moments D3Q19
//! supports and is the more stable of the two in the comparisons of Coreixas et al. (2017) and
//! of C. Coreixas, B. Chopard, J. Latt, *Comprehensive comparison of collision models in the
//! lattice Boltzmann framework*, Phys. Rev. E 100:033305, 2019, doi:10.1103/PhysRevE.100.033305,
//! and in our own runs (`docs/lbm-collision.md`): with the projected form a lid-driven cavity at
//! τ⁺ = 0.51 and a shear layer at lattice speed 0.2 failed where the recursive one ran.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.
//!
//! # Open boundaries
//! A pressure outlet under this collision holds the equilibrium of its density and normal
//! velocity ([`crate::flow`]); velocity inlets keep Zou and He's rule.

use crate::lattice::{C, Q, W};

/// The third-order Hermite polynomials D3Q19 holds, per direction, combined in the pairs that are
/// orthogonal under the weights: `[x+, x−, y+, y−, z+, z−]` with
/// `x± = H_xyy ± H_xzz`, `y± = H_xxy ± H_yzz`, `z± = H_xxz ± H_yyz`, `H_abb = c_a c_b² − c_a/3`.
const fn hermite3(q: usize) -> [f64; 6] {
    let c = C[q];
    let (x, y, z) = (c[0] as f64, c[1] as f64, c[2] as f64);
    let third = 1.0 / 3.0;
    let (xyy, xzz) = (x * y * y - x * third, x * z * z - x * third);
    let (xxy, yzz) = (x * x * y - y * third, y * z * z - y * third);
    let (xxz, yyz) = (x * x * z - z * third, y * y * z - z * third);
    [xyy + xzz, xyy - xzz, xxy + yzz, xxy - yzz, xxz + yyz, xxz - yyz]
}

/// [`hermite3`] of every direction.
const HERMITE3: [[f64; 6]; Q] = {
    let mut h = [[0.0; 6]; Q];
    let mut q = 0;
    while q < Q {
        h[q] = hermite3(q);
        q += 1;
    }
    h
};

/// The weights' inverse norms of the pairs in [`hermite3`]: 27/2 for the sums, 9/2 for the
/// differences.
const NORM3: [f64; 2] = [13.5, 4.5];

/// Post-collision populations of one cell: `fq` before the collision, its density `rho`, its
/// velocity `u` (with the half force) and its force per volume `force` (`forced`: whether there is
/// one anywhere in the flow; without, the force terms are skipped), relaxation rate `omega`.
/// With a Smagorinsky constant `smagorinsky` above zero the cell relaxes with
/// [`crate::turbulence::omega`] of its own flux instead. Returns the rate it relaxed with.
#[inline(always)]
pub fn collide(fq: &[f64; Q], rho: f64, u: [f64; 3], force: [f64; 3], omega: f64, smagorinsky: f64, forced: bool, out: &mut [f64; Q]) -> f64 {
    // The momentum flux Σ c c f: xx, yy, zz, xy, xz, yz (opposite directions share c c, so by
    // pairs).
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
    let cs2 = crate::lattice::CS2;
    let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    // P = Π − ρ (c_s² I + u u) + ½ (u F + F u).
    let mut p = [0.0f64; 6];
    for (k, &(a, b)) in pairs.iter().enumerate() {
        let iso = if a == b { cs2 } else { 0.0 };
        p[k] = pi[k] - rho * (iso + u[a] * u[b]) + if forced { 0.5 * (u[a] * force[b] + force[a] * u[b]) } else { 0.0 };
    }
    let omega = if smagorinsky > 0.0 { crate::turbulence::omega(&p, rho, omega, smagorinsky) } else { omega };
    let [pxx, pyy, pzz, pxy, pxz, pyz] = p;
    let trace = pxx + pyy + pzz;
    let relax = 1.0 - omega;
    let [ux, uy, uz] = u;
    // Third-order coefficients: equilibrium ρ u_a u_b u_c plus (1 − ω) times the recursive
    // non-equilibrium u_a P_bc + u_b P_ac + u_c P_ab.
    let a3 = |e: f64, n: f64| e + relax * n;
    let xyy = a3(rho * ux * uy * uy, 2.0 * uy * pxy + ux * pyy);
    let xzz = a3(rho * ux * uz * uz, 2.0 * uz * pxz + ux * pzz);
    let xxy = a3(rho * ux * ux * uy, 2.0 * ux * pxy + uy * pxx);
    let yzz = a3(rho * uy * uz * uz, 2.0 * uz * pyz + uy * pzz);
    let xxz = a3(rho * ux * ux * uz, 2.0 * ux * pxz + uz * pxx);
    let yyz = a3(rho * uy * uy * uz, 2.0 * uy * pyz + uz * pyy);
    let b = [
        NORM3[0] * (xyy + xzz),
        NORM3[1] * (xyy - xzz),
        NORM3[0] * (xxy + yzz),
        NORM3[1] * (xxy - yzz),
        NORM3[0] * (xxz + yyz),
        NORM3[1] * (xxz - yyz),
    ];
    let uu = ux * ux + uy * uy + uz * uz;
    let uf = ux * force[0] + uy * force[1] + uz * force[2];
    // The rest population, then each pair of opposite directions: the even terms (orders 0, 2
    // and the force's second) are the same for both, the odd ones (orders 1, 3 and the force's
    // first) change sign.
    let second0 = -1.5 * trace;
    out[0] = W[0] * (rho * (1.0 - 1.5 * uu) + relax * second0) - if forced { 1.5 * W[0] * uf } else { 0.0 };
    for q in (1..Q).step_by(2) {
        let c = C[q];
        let (x, y, z) = (f64::from(c[0]), f64::from(c[1]), f64::from(c[2]));
        let w = W[q];
        let cu = x * ux + y * uy + z * uz;
        let h = &HERMITE3[q];
        let third = h[0] * b[0] + h[1] * b[1] + h[2] * b[2] + h[3] * b[3] + h[4] * b[4] + h[5] * b[5];
        // 1/(2 c_s⁴) Q:P = 9/2 (c·P·c − tr P / 3).
        let cpc = x * x * pxx + y * y * pyy + z * z * pzz + 2.0 * (x * y * pxy + x * z * pxz + y * z * pyz);
        let second = 4.5 * (cpc - trace / 3.0);
        let mut even = w * (rho * (1.0 + 4.5 * cu * cu - 1.5 * uu) + relax * second);
        let mut odd = w * (rho * 3.0 * cu + third);
        if forced {
            let cf = x * force[0] + y * force[1] + z * force[2];
            even += 0.5 * w * (9.0 * cu * cf - 3.0 * uf);
            odd += 1.5 * w * cf;
        }
        out[q] = even + odd;
        out[q + 1] = even - odd;
    }
    omega
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::{CS2, c};

    /// Σ f, Σ f c, Σ f c c (xx yy zz xy xz yz) and the six supported third moments.
    fn moments(f: &[f64; Q]) -> (f64, [f64; 3], [f64; 6], [f64; 6]) {
        let (mut r, mut j, mut pi, mut m3) = (0.0, [0.0; 3], [0.0; 6], [0.0; 6]);
        for q in 0..Q {
            let [x, y, z] = c(q);
            r += f[q];
            for (k, v) in [x, y, z].iter().enumerate() {
                j[k] += f[q] * v;
            }
            for (k, v) in [x * x, y * y, z * z, x * y, x * z, y * z].iter().enumerate() {
                pi[k] += f[q] * v;
            }
            for (k, v) in [x * y * y, x * z * z, x * x * y, y * z * z, x * x * z, y * y * z].iter().enumerate() {
                m3[k] += f[q] * v;
            }
        }
        (r, j, pi, m3)
    }

    #[test]
    fn the_third_order_hermite_pairs_are_orthogonal_with_the_stated_norms() {
        for a in 0..6 {
            for b in 0..6 {
                let s: f64 = (0..Q).map(|q| W[q] * HERMITE3[q][a] * HERMITE3[q][b]).sum();
                let want = if a != b {
                    0.0
                } else if a % 2 == 0 {
                    2.0 / 27.0
                } else {
                    6.0 / 27.0
                };
                assert!((s - want).abs() < 1e-15, "{a} {b}: {s}");
            }
            // Orthogonal to the lower orders too.
            let s0: f64 = (0..Q).map(|q| W[q] * HERMITE3[q][a]).sum();
            assert!(s0.abs() < 1e-15);
            for k in 0..3 {
                let s1: f64 = (0..Q).map(|q| W[q] * HERMITE3[q][a] * c(q)[k]).sum();
                assert!(s1.abs() < 1e-15);
            }
        }
    }

    #[test]
    fn the_collision_keeps_the_hydrodynamic_moments_and_relaxes_the_flux() {
        // An arbitrary state off equilibrium, with a force.
        let mut f = [0.0; Q];
        for q in 0..Q {
            f[q] = crate::lattice::equilibrium(q, 1.02, [0.03, -0.02, 0.01]) * (1.0 + 0.01 * ((q * 7 % 5) as f64 - 2.0));
        }
        let force = [1e-4, -2e-4, 3e-5];
        let omega = 1.9;
        let (rho, j, pi, _) = moments(&f);
        let u = [0, 1, 2].map(|k| (j[k] + 0.5 * force[k]) / rho);
        let mut out = [0.0; Q];
        collide(&f, rho, u, force, omega, 0.0, true, &mut out);
        let (r2, j2, pi2, m3) = moments(&out);
        assert!((r2 - rho).abs() < 1e-15, "mass");
        for k in 0..3 {
            assert!((j2[k] - (j[k] + force[k])).abs() < 1e-15, "momentum gains the force");
        }
        let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let iso = if a == b { CS2 } else { 0.0 };
            let eq = rho * (iso + u[a] * u[b]);
            let uf = 0.5 * (u[a] * force[b] + force[a] * u[b]);
            let p = pi[k] - eq + uf;
            let want = eq + (1.0 - omega) * p + uf;
            assert!((pi2[k] - want).abs() < 1e-15, "flux {k}: {} vs {want}", pi2[k]);
        }
        // Third moments: ρ u_a u_b u_c plus (1 − ω) times the recursive part.
        let p = |a: usize, b: usize| {
            let k = pairs.iter().position(|&(x, y)| (x, y) == (a.min(b), a.max(b))).unwrap();
            let iso = if a == b { CS2 } else { 0.0 };
            pi[k] - rho * (iso + u[a] * u[b]) + 0.5 * (u[a] * force[b] + force[a] * u[b])
        };
        let third = |a: usize, b: usize, c: usize| {
            let e = rho * u[a] * u[b] * u[c];
            let n = u[a] * p(b, c) + u[b] * p(a, c) + u[c] * p(a, b);
            // the Hermite moment plus the c_s² parts of the raw moment c_a c_b c_c
            let h = e + (1.0 - omega) * n;
            let raw_extra = |i: usize, jj: usize, kk: usize| if jj == kk { CS2 * (j2[i]) } else { 0.0 };
            h + raw_extra(a, b, c)
        };
        let want = [third(0, 1, 1), third(0, 2, 2), third(1, 0, 0), third(1, 2, 2), third(2, 0, 0), third(2, 1, 1)];
        for k in 0..6 {
            assert!((m3[k] - want[k]).abs() < 1e-14, "third moment {k}: {} vs {}", m3[k], want[k]);
        }
    }

    #[test]
    fn an_equilibrium_without_force_is_a_fixed_point_up_to_its_third_order_part() {
        // A second-order equilibrium at rest is the scheme's own equilibrium.
        let mut f = [0.0; Q];
        for q in 0..Q {
            f[q] = crate::lattice::equilibrium(q, 0.97, [0.0; 3]);
        }
        let mut out = [0.0; Q];
        collide(&f, 0.97, [0.0; 3], [0.0; 3], 1.7, 0.0, false, &mut out);
        for q in 0..Q {
            assert!((out[q] - f[q]).abs() < 1e-16, "{q}");
        }
    }
}
