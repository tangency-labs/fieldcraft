//! #820: linear (von Neumann) stability of the D3Q19 two-relaxation-time collision about a uniform
//! stream, on every mode of the 16³ periodic box of the #820 reproduction.
//!
//! The same linearisation as Sterling & Chen 1996 (Eq. 31–33) and Wissocq, Sagaut & Boussuge 2019
//! (J. Comput. Phys. 380:311–333, doi:10.1016/j.jcp.2018.12.015, App. B.2): about ρ = 1, u = (U, 0, 0)
//! a plane wave of wavevector k advances by G(k) = diag(e^{−ik·c_q}) · (I − ω⁺(P⁺ − J⁺) − ω⁻(P⁻ − J⁻)),
//! with P± the symmetric and antisymmetric projections and J the equilibrium's Jacobian. A mode
//! grows when the spectral radius of G(k) exceeds one; the growth rate per step is taken from
//! ‖G^N‖^{1/N} with N = 2²⁴ by repeated squaring, which a non-normal transient cannot mimic. The
//! staggered k = π modes are neutral (Ginzburg, Verhaeghe & d'Humières 2008, Sec. 2 and 5.4), so
//! "stable" is a rate at round-off.
//!
//! What it shows at #820's lattice speed 0.05 (Mach 0.087) with Λ = 3/16: the stream itself grows
//! for τ⁺ ≤ 0.502 and is neutral from 0.5025; at τ⁺ = 0.5005 it is neutral only up to U = 0.01, or
//! with a near-BGK Λ (τ⁻ ≈ τ⁺). The growing modes are oblique, as Wissocq et al. find for BGK.
use fieldcraft_lbm::lattice::{C, Q, W, opposite};

#[derive(Clone, Copy)]
struct Complex(f64, f64);

impl Complex {
    fn mul(self, o: Complex) -> Complex {
        Complex(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn add(self, o: Complex) -> Complex {
        Complex(self.0 + o.0, self.1 + o.1)
    }
    fn scale(self, s: f64) -> Complex {
        Complex(self.0 * s, self.1 * s)
    }
    fn norm2(self) -> f64 {
        self.0 * self.0 + self.1 * self.1
    }
}

/// The largest growth rate per step over the modes of an n³ box (k = 0 left out).
fn growth(tau: f64, magic: f64, speed: f64, n: usize) -> f64 {
    let c = |q: usize| C[q].map(f64::from);
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let u = [speed, 0.0, 0.0];
    let (wp, wm) = (1.0 / tau, 1.0 / (0.5 + magic / (tau - 0.5)));
    // ∂f_eq,q/∂f_p at ρ = 1, j = u, with f_eq = w(ρ + 3c·j + 4.5(c·j)²/ρ − 1.5 j²/ρ).
    let jacobian = |q: usize, p: usize| {
        let cu = dot(c(q), u);
        let d_rho = W[q] * (1.0 - 4.5 * cu * cu + 1.5 * dot(u, u));
        let d_j: [f64; 3] = std::array::from_fn(|a| W[q] * (3.0 * c(q)[a] + 9.0 * cu * c(q)[a] - 3.0 * u[a]));
        d_rho + dot(c(p), d_j)
    };
    let mut collide = [[0.0; Q]; Q];
    for (q, row) in collide.iter_mut().enumerate() {
        let o = opposite(q);
        for (p, x) in row.iter_mut().enumerate() {
            let id = |i: usize| if i == p { 1.0 } else { 0.0 };
            let (plus, minus) = (0.5 * (id(q) + id(o)), 0.5 * (id(q) - id(o)));
            let (j_plus, j_minus) = (0.5 * (jacobian(q, p) + jacobian(o, p)), 0.5 * (jacobian(q, p) - jacobian(o, p)));
            *x = id(q) - wp * (plus - j_plus) - wm * (minus - j_minus);
        }
    }
    let norm = |m: &[Complex]| m.iter().map(|x| x.norm2()).sum::<f64>().sqrt();
    let mut worst = f64::NEG_INFINITY;
    for m in 1..n * n * n {
        let k = [m % n, m / n % n, m / (n * n)].map(|i| 2.0 * std::f64::consts::PI * i as f64 / n as f64);
        let mut g = vec![Complex(0.0, 0.0); Q * Q];
        for q in 0..Q {
            let phase = -dot(k, c(q));
            for p in 0..Q {
                g[q * Q + p] = Complex(phase.cos(), phase.sin()).scale(collide[q][p]);
            }
        }
        // log ‖G^(2^s)‖ / 2^s, renormalising after each squaring.
        let first = norm(&g);
        g.iter_mut().for_each(|x| *x = x.scale(1.0 / first));
        let (mut log, mut len) = (first.ln(), 1.0);
        for _ in 0..24 {
            let mut square = vec![Complex(0.0, 0.0); Q * Q];
            for i in 0..Q {
                for l in 0..Q {
                    let a = g[i * Q + l];
                    for j in 0..Q {
                        square[i * Q + j] = square[i * Q + j].add(a.mul(g[l * Q + j]));
                    }
                }
            }
            let s = norm(&square);
            square.iter_mut().for_each(|x| *x = x.scale(1.0 / s));
            g = square;
            log = 2.0 * log + s.ln();
            len *= 2.0;
        }
        worst = worst.max(log / len);
    }
    worst
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_trt_uniform_stream_is_linearly_unstable_close_to_one_half_at_the_820_mach_number() {
    const DEFAULT: f64 = 3.0 / 16.0;
    // Round-off of the neutral staggered modes stays below this; real growth is far above it.
    const NEUTRAL: f64 = 1e-7;
    for (tau, magic, speed, unstable) in [
        (0.5005, DEFAULT, 0.05, true),
        (0.501, DEFAULT, 0.05, true),
        (0.502, DEFAULT, 0.05, true),
        (0.5025, DEFAULT, 0.05, false),
        (0.505, DEFAULT, 0.05, false),
        (0.5005, DEFAULT, 0.01, false),
        (0.5005, DEFAULT, 0.02, true),
        (0.5005, 0.25, 0.05, true),
        (0.5005, 1e-4, 0.05, true),
        (0.5005, 1e-5, 0.05, false),
    ] {
        let rate = growth(tau, magic, speed, 16);
        println!("τ⁺ {tau} Λ {magic:e} U {speed}: largest growth {rate:.3e} per step");
        assert_eq!(rate > NEUTRAL, unstable, "τ⁺ {tau} Λ {magic} U {speed}: {rate:e}");
    }
}
