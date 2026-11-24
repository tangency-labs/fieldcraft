//! Linear (von Neumann) stability of the D3Q19 multiple-relaxation-time collision (`collision: "mrt"`,
//! `crates/lbm/src/mrt.rs`) about a uniform state: the origin of [`fieldcraft_lbm::mrt::MIN_VISCOSITY`].
//!
//! The analysis is Lallemand & Luo's (Phys. Rev. E 61:6546–6562, 2000, Sec. IV), the one d'Humières et al.
//! (Phil. Trans. R. Soc. Lond. A 360:437–451, 2002, appendix A) cite for their rates, in the form of
//! `trt_stability.rs`: about ρ = 1 and a velocity u, a plane wave of wavevector k advances by
//! G(k) = diag(e^{−ik·c_q}) · J, with J the Jacobian of the collision at the equilibrium (central differences
//! of `mrt::collide` itself, so the test follows the collision as it is implemented). A mode grows when the
//! spectral radius of G(k) exceeds one; the growth rate per step is ln ‖G^N‖ / N with N = 2²⁴ by repeated
//! squaring.
//!
//! What it shows: with the rates 1.19, 1.4, 1.2, 1.98 and the single-relaxation-time-equivalent equilibria the
//! state of rest is unstable below a lattice viscosity of 1.18008e-3 (τ = 0.50354). The last mode to become
//! stable is at k = π (0.2537, 0.4207, 0.4207) and its images under the lattice's symmetries; at Han, Ooka &
//! Kikumoto's A20 viscosity, 5e-5, the fastest mode grows by 7.9e-3 per step. The bound rises with the speed of
//! the stream: 1.2295e-3 at lattice Mach 0.1 (no wave grows above it, over 45 directions of the stream) and
//! 1.4016e-3 at 0.2, along (1, ½, ½); beyond Mach 0.1 slower waves grow above that bound too.
//! d'Humières et al. state ν > 2.54e-3 at speeds up to 0.19 for the same rates with their optimised equilibria.
use fieldcraft_lbm::lattice::{C, Q, equilibrium};
use fieldcraft_lbm::mrt::{MIN_VISCOSITY, MIN_VISCOSITY_AT_MACH_0_1, collide};
use std::f64::consts::PI;

/// A rate within this of zero is neutral: the floor of ln ‖G^N‖ / N at N = 2²⁴ (the norm of a power exceeds
/// the spectral radius's by a bounded factor, whose logarithm is divided by N).
const NEUTRAL: f64 = 1e-6;

/// The Jacobian of the collision at the equilibrium of density 1 and velocity `u`: row q, column p is
/// ∂f*_q/∂f_p, by central differences of `mrt::collide` (the collision is quadratic in the momentum over the
/// density, so the step's truncation error is of order h² ≈ 1e-10).
fn jacobian(nu: f64, u: [f64; 3]) -> [[f64; Q]; Q] {
    let omega = 1.0 / (3.0 * nu + 0.5);
    let post = |f: &[f64; Q]| {
        let rho: f64 = f.iter().sum();
        let v: [f64; 3] = std::array::from_fn(|a| (0..Q).map(|q| f64::from(C[q][a]) * f[q]).sum::<f64>() / rho);
        let mut out = [0.0; Q];
        collide(f, rho, v, [0.0; 3], omega, false, &mut out);
        out
    };
    let feq: [f64; Q] = std::array::from_fn(|q| equilibrium(q, 1.0, u));
    let h = 1e-5;
    let mut j = [[0.0; Q]; Q];
    for p in 0..Q {
        let (mut up, mut down) = (feq, feq);
        up[p] += h;
        down[p] -= h;
        let (a, b) = (post(&up), post(&down));
        for q in 0..Q {
            j[q][p] = (a[q] - b[q]) / (2.0 * h);
        }
    }
    j
}

/// The growth rate per step of the wave `k` under the collision with Jacobian `j`.
fn growth(j: &[[f64; Q]; Q], k: [f64; 3]) -> f64 {
    // (re, im) pairs, row-major.
    let mut g = vec![(0.0f64, 0.0f64); Q * Q];
    for q in 0..Q {
        let phase = -(0..3).map(|a| k[a] * f64::from(C[q][a])).sum::<f64>();
        for p in 0..Q {
            g[q * Q + p] = (phase.cos() * j[q][p], phase.sin() * j[q][p]);
        }
    }
    let norm = |m: &[(f64, f64)]| m.iter().map(|x| x.0 * x.0 + x.1 * x.1).sum::<f64>().sqrt();
    let first = norm(&g);
    g.iter_mut().for_each(|x| *x = (x.0 / first, x.1 / first));
    let (mut log, mut len) = (first.ln(), 1.0);
    for _ in 0..24 {
        let mut square = vec![(0.0f64, 0.0f64); Q * Q];
        for i in 0..Q {
            for l in 0..Q {
                let a = g[i * Q + l];
                for m in 0..Q {
                    let b = g[l * Q + m];
                    let s = &mut square[i * Q + m];
                    s.0 += a.0 * b.0 - a.1 * b.1;
                    s.1 += a.0 * b.1 + a.1 * b.0;
                }
            }
        }
        let s = norm(&square);
        square.iter_mut().for_each(|x| *x = (x.0 / s, x.1 / s));
        g = square;
        log = 2.0 * log + s.ln();
        len *= 2.0;
    }
    log / len
}

/// The wave that is the last to become stable at rest, and the fastest one at A20's viscosity.
const LAST: [f64; 3] = [0.2537 * PI, 0.4207 * PI, 0.4207 * PI];
const FASTEST_AT_A20: [f64; 3] = [0.2002 * PI, 0.5928 * PI, 0.5928 * PI];

/// The state of rest is unstable one per cent below `MIN_VISCOSITY` and stable one per cent above it, on the
/// wave that decides it; at A20's viscosity (τ = 0.50015) it grows by 7.9e-3 per step.
#[test]
fn the_state_of_rest_is_unstable_below_the_bound_and_stable_above_it() {
    let rest = [0.0; 3];
    let below = growth(&jacobian(0.99 * MIN_VISCOSITY, rest), LAST);
    let above = growth(&jacobian(1.01 * MIN_VISCOSITY, rest), LAST);
    let a20 = growth(&jacobian(5.0e-5, rest), FASTEST_AT_A20);
    println!("MRT at rest: growth per step {below:+.3e} at 0.99 ν_c, {above:+.3e} at 1.01 ν_c, {a20:+.3e} at ν 5e-5");
    assert!(below > NEUTRAL, "{below:e}");
    assert!(above < -NEUTRAL, "{above:e}");
    assert!((a20 - 7.9e-3).abs() < 1e-4, "{a20:e}");
}

/// Every wave of a 32³ box (k = 2π i / 32) at rest: none grows one per cent above the bound. With a stream
/// the bound rises: at lattice Mach 0.1 along (1, ½, ½), the direction that raises it most,
/// the deciding wave grows one per cent below `MIN_VISCOSITY_AT_MACH_0_1` and decays one per cent above it;
/// at Mach 0.2 the same at 1.4016e-3 (there slower waves, of 3e-5 per step, still grow up to ν ≈ 3e-3:
/// `docs/manual/fluid/collision.md`).
#[test]
fn no_wave_grows_above_the_bound_and_the_bound_rises_with_the_stream() {
    let n = 32usize;
    let j = jacobian(1.01 * MIN_VISCOSITY, [0.0; 3]);
    let mut worst = f64::NEG_INFINITY;
    // k and −k have the same rate: half the box.
    for m in 1..n * n * (n / 2 + 1) {
        let k = [m % n, m / n % n, m / (n * n)].map(|i| 2.0 * PI * i as f64 / n as f64);
        worst = worst.max(growth(&j, k));
    }
    println!("MRT at rest, 1.01 ν_c, {n}³ waves: largest growth {worst:+.3e} per step");
    assert!(worst <= NEUTRAL, "{worst:e}");
    let along = |speed: f64| {
        let d = [1.0f64, 0.5, 0.5];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        d.map(|c| speed * c / len)
    };
    // (lattice Mach number, the bound there, the deciding wave over π)
    for (mach, bound, wave) in [(0.1, MIN_VISCOSITY_AT_MACH_0_1, [0.2561, 0.4228, 0.4228]), (0.2, 1.4016e-3, [0.2644, 0.4266, 0.4266])] {
        let u = along(mach / 3f64.sqrt());
        let k = wave.map(|c| c * PI);
        let (below, above) = (growth(&jacobian(0.99 * bound, u), k), growth(&jacobian(1.01 * bound, u), k));
        println!("MRT at lattice Mach {mach}: growth {below:+.3e} at 0.99 × {bound:e}, {above:+.3e} at 1.01 ×");
        assert!(below > NEUTRAL && above < -NEUTRAL, "Mach {mach}: {below:e} {above:e}");
    }
}
