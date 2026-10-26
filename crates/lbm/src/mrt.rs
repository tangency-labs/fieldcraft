//! The D3Q19 multiple-relaxation-time collision.
//!
//! The populations are mapped to 19 orthogonal moments, each relaxed towards its equilibrium at
//! its own rate, and mapped back (D. d'Humières, I. Ginzburg, M. Krafczyk, P. Lallemand,
//! L.-S. Luo, *Multiple-relaxation-time lattice Boltzmann models in three dimensions*, Phil.
//! Trans. R. Soc. Lond. A 360:437–451, 2002, doi:10.1098/rsta.2001.0955, appendix A). The
//! moments are the integer polynomials of the lattice velocity c (with c² = |c|²): the density
//! 1; the energy 19 c² − 30 and its square (21 c⁴ − 53 c² + 24)/2; the momentum c_a and the
//! energy flux (5 c² − 9) c_a; the stresses 3 c_x² − c², c_y² − c_z², c_x c_y, c_y c_z, c_x c_z and
//! the two quartic modes (3 c² − 5)(3 c_x² − c²), (3 c² − 5)(c_y² − c_z²); and the three cubic
//! modes (c_y² − c_z²) c_x, (c_z² − c_x²) c_y, (c_x² − c_y²) c_z. They are orthogonal, so the inverse
//! map is the transpose over each moment's squared norm. Built from the polynomials, the map
//! does not depend on the order of the lattice's directions.
//!
//! The equilibria are the second-order ones equivalent to the single-relaxation-time model
//! (e = ρ (−11 + 19 u²), ε = 3 ρ − 11/2 ρ u², q = −2/3 ρ u, the stresses ρ u_a u_b in the same
//! combinations, the quartic modes −½ of the stresses', the cubic ones 0), and the rates those
//! of O. Malaspinas, P. Sagaut, *Wall model for large-eddy simulation based on the lattice
//! Boltzmann method*, J. Comput. Phys. 275:25–40, 2014 (doi:10.1016/j.jcp.2014.06.020,
//! appendix A, from d'Humières et al.'s linear stability analysis): 1.19 for the energy, 1.4
//! for its square and the quartic modes, 1.2 for the energy flux, 1.98 for the cubic modes, and
//! the viscous rate 1 / (3 ν + ½) for the five stresses, which Smagorinsky's model replaces per
//! cell (as in the paper's eqs. 47–48). Guo's force enters as `(I − S/2)` times its moments, as
//! in the paper's eq. 34, with the equilibrium at the velocity that includes half the force.
//!
//! Where this differs from the papers: Malaspinas and Sagaut's eq. 10 for the energy square's
//! equilibrium is printed as (11/2) ρ (3 − u²); this module uses d'Humières et al.'s
//! single-relaxation-time-equivalent 3 ρ − 11/2 ρ u² (w_ε = 3, w_εj = −11/2), the choice the
//! paper's other equilibria (the quartic modes at −½, its w_xx) correspond to. The equilibria use
//! j·j/ρ (the paper's ρ u u), not d'Humières et al.'s j·j/ρ₀.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.

use crate::lattice::{C, Q};

/// The value of moment `k`'s polynomial at velocity `c` (integers).
const fn polynomial(k: usize, c: [i32; 3]) -> i32 {
    let (x, y, z) = (c[0], c[1], c[2]);
    let c2 = x * x + y * y + z * z;
    match k {
        0 => 1,
        1 => 19 * c2 - 30,
        2 => (21 * c2 * c2 - 53 * c2 + 24) / 2,
        3 => x,
        4 => (5 * c2 - 9) * x,
        5 => y,
        6 => (5 * c2 - 9) * y,
        7 => z,
        8 => (5 * c2 - 9) * z,
        9 => 3 * x * x - c2,
        10 => (3 * c2 - 5) * (3 * x * x - c2),
        11 => y * y - z * z,
        12 => (3 * c2 - 5) * (y * y - z * z),
        13 => x * y,
        14 => y * z,
        15 => x * z,
        16 => (y * y - z * z) * x,
        17 => (z * z - x * x) * y,
        _ => (x * x - y * y) * z,
    }
}

/// The map from populations to moments: row `k` is moment `k`'s polynomial over the lattice's
/// directions.
const M: [[i32; Q]; Q] = {
    let mut m = [[0; Q]; Q];
    let mut k = 0;
    while k < Q {
        let mut q = 0;
        while q < Q {
            m[k][q] = polynomial(k, C[q]);
            q += 1;
        }
        k += 1;
    }
    m
};

/// Each moment's squared norm Σ_q M_kq², which the inverse map divides by.
const NORM: [i32; Q] = {
    let mut n = [0; Q];
    let mut k = 0;
    while k < Q {
        let mut q = 0;
        while q < Q {
            n[k] += M[k][q] * M[k][q];
            q += 1;
        }
        k += 1;
    }
    n
};

/// The moments relaxed at the viscous rate: the five stresses.
const VISCOUS: [usize; 5] = [9, 11, 13, 14, 15];

/// Each moment's relaxation rate, the viscous ones set per cell (the conserved moments' rates do
/// not act: their equilibria are themselves).
const RATES: [f64; Q] = [0.0, 1.19, 1.4, 0.0, 1.2, 0.0, 1.2, 0.0, 1.2, 0.0, 1.4, 0.0, 1.4, 0.0, 0.0, 0.0, 1.98, 1.98, 1.98];

/// The equilibrium moments at density `rho` and velocity `u`.
#[inline(always)]
fn equilibrium(rho: f64, u: [f64; 3]) -> [f64; Q] {
    let [ux, uy, uz] = u;
    let uu = ux * ux + uy * uy + uz * uz;
    let pxx = rho * (2.0 * ux * ux - uy * uy - uz * uz);
    let pww = rho * (uy * uy - uz * uz);
    [
        rho,
        rho * (-11.0 + 19.0 * uu),
        3.0 * rho - 5.5 * rho * uu,
        rho * ux,
        -2.0 / 3.0 * rho * ux,
        rho * uy,
        -2.0 / 3.0 * rho * uy,
        rho * uz,
        -2.0 / 3.0 * rho * uz,
        pxx,
        -0.5 * pxx,
        pww,
        -0.5 * pww,
        rho * ux * uy,
        rho * uy * uz,
        rho * ux * uz,
        0.0,
        0.0,
        0.0,
    ]
}

/// Post-collision populations of one cell: `fq` before the collision, its density `rho`, its
/// velocity `u` (with the half force) and its force per volume `force` (`forced`: whether there
/// is one anywhere in the flow), the stresses relaxing at `omega` (the fluid's own rate, or the
/// cell's under Smagorinsky's model or the wall model).
#[inline(always)]
pub fn collide(fq: &[f64; Q], rho: f64, u: [f64; 3], force: [f64; 3], omega: f64, forced: bool, out: &mut [f64; Q]) {
    let eq = equilibrium(rho, u);
    // Guo's force in the populations, w_q (3 (c − u) + 9 (c·u) c)·F, then its moments.
    let mut source = [0.0f64; Q];
    if forced {
        for (q, s) in source.iter_mut().enumerate() {
            let c = C[q];
            let cv = [f64::from(c[0]), f64::from(c[1]), f64::from(c[2])];
            let cu = cv[0] * u[0] + cv[1] * u[1] + cv[2] * u[2];
            *s = crate::lattice::W[q] * (0..3).map(|a| (3.0 * (cv[a] - u[a]) + 9.0 * cu * cv[a]) * force[a]).sum::<f64>();
        }
    }
    // The change of each moment: −s (m − m_eq) + (1 − s/2) M F, over its squared norm.
    let mut change = [0.0f64; Q];
    for k in 0..Q {
        let row = &M[k];
        let (mut m, mut mf) = (0.0, 0.0);
        for q in 0..Q {
            m += f64::from(row[q]) * fq[q];
            mf += f64::from(row[q]) * source[q];
        }
        let s = if VISCOUS.contains(&k) { omega } else { RATES[k] };
        change[k] = (-s * (m - eq[k]) + (1.0 - 0.5 * s) * mf) / f64::from(NORM[k]);
    }
    for q in 0..Q {
        let mut d = 0.0;
        for k in 0..Q {
            d += f64::from(M[k][q]) * change[k];
        }
        out[q] = fq[q] + d;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::W;

    /// The rows are orthogonal, and the conserved moments are the density and the momentum.
    #[test]
    fn the_moments_are_orthogonal() {
        for a in 0..Q {
            for b in 0..Q {
                let dot: i32 = (0..Q).map(|q| M[a][q] * M[b][q]).sum();
                assert_eq!(dot == 0, a != b, "{a} {b}: {dot}");
            }
        }
        assert_eq!(NORM[0], 19);
        for q in 0..Q {
            assert_eq!([M[3][q], M[5][q], M[7][q]], C[q]);
        }
    }

    /// The equilibrium populations have the equilibrium moments, so the collision leaves an
    /// equilibrium unforced state as it is, conserving mass and momentum.
    #[test]
    fn an_equilibrium_is_left_as_it_is() {
        let (rho, u) = (1.03, [0.04, -0.02, 0.01]);
        let fq: [f64; Q] = std::array::from_fn(|q| crate::lattice::equilibrium(q, rho, u));
        let mut out = [0.0; Q];
        collide(&fq, rho, u, [0.0; 3], 1.7, false, &mut out);
        for q in 0..Q {
            assert!((out[q] - fq[q]).abs() < 1e-15, "{q}: {} vs {}", out[q], fq[q]);
        }
        let eq = equilibrium(rho, u);
        for (k, e) in eq.iter().enumerate() {
            let m: f64 = (0..Q).map(|q| f64::from(M[k][q]) * fq[q]).sum();
            assert!((m - e).abs() < 1e-14 * (1.0 + e.abs()), "moment {k}: {m} vs {e}");
        }
        assert!((W.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    }

    /// With every rate equal to the viscous one the collision is the single-relaxation-time one:
    /// the non-equilibrium relaxes by (1 − ω) in every direction, only the rates differ here.
    /// The stresses' non-equilibrium (a simple shear) relaxes at the viscous rate.
    #[test]
    fn a_shear_relaxes_at_the_viscous_rate() {
        let (rho, u, omega) = (1.0, [0.02, 0.0, 0.0], 1.6);
        let feq: [f64; Q] = std::array::from_fn(|q| crate::lattice::equilibrium(q, rho, u));
        // A pure xy stress: w_q 9/2 (c_x c_y) P_xy.
        let pxy = 1e-3;
        let fq: [f64; Q] = std::array::from_fn(|q| feq[q] + W[q] * 4.5 * 2.0 * f64::from(C[q][0]) * f64::from(C[q][1]) * pxy);
        let mut out = [0.0; Q];
        collide(&fq, rho, u, [0.0; 3], omega, false, &mut out);
        for q in 0..Q {
            let want = feq[q] + (1.0 - omega) * (fq[q] - feq[q]);
            assert!((out[q] - want).abs() < 1e-15, "{q}: {} vs {want}", out[q]);
        }
    }

    /// Guo's force adds its momentum: Σ c f* − Σ c f = F.
    #[test]
    fn the_force_adds_its_momentum() {
        let (rho, force) = (1.0, [1e-4, -2e-4, 3e-5]);
        let u0 = [0.01, 0.02, 0.0];
        let fq: [f64; Q] = std::array::from_fn(|q| crate::lattice::equilibrium(q, rho, u0));
        let j: [f64; 3] = std::array::from_fn(|a| (0..Q).map(|q| f64::from(C[q][a]) * fq[q]).sum());
        let u = [0, 1, 2].map(|a| (j[a] + 0.5 * force[a]) / rho);
        let mut out = [0.0; Q];
        collide(&fq, rho, u, force, 1.5, true, &mut out);
        for a in 0..3 {
            let after: f64 = (0..Q).map(|q| f64::from(C[q][a]) * out[q]).sum();
            assert!((after - j[a] - force[a]).abs() < 1e-15, "{a}");
        }
        assert!((out.iter().sum::<f64>() - fq.iter().sum::<f64>()).abs() < 1e-15);
    }
}
