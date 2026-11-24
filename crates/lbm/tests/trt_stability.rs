//! #820, #1063: linear (von Neumann) stability of the D3Q19 two-relaxation-time collision about a uniform
//! stream. The 6.18 along (1, 1, 0) is checked here; [`fieldcraft_lbm::setup::TRT_STREAM_LIMIT`] (5.6) rests on
//! a search over the stream's direction (a κ onset, notes #1063 section 7), and the last test here brackets it
//! at one direction, it does not derive it.
//!
//! The same linearisation as Sterling & Chen 1996 (Eq. 31–33) and Wissocq, Sagaut & Boussuge 2019
//! (J. Comput. Phys. 380:311–333, doi:10.1016/j.jcp.2018.12.015, App. B.2): about ρ = 1 and a velocity u
//! a plane wave of wavevector k advances by G(k) = diag(e^{−ik·c_q}) · (I − ω⁺(P⁺ − J⁺) − ω⁻(P⁻ − J⁻)),
//! with P± the symmetric and antisymmetric projections and J the equilibrium's Jacobian. A mode
//! grows when the spectral radius of G(k) exceeds one; the growth rate per step is taken from
//! ‖G^N‖^{1/N} with N = 2²⁴ by repeated squaring, which a non-normal transient cannot mimic. The
//! staggered k = π modes are neutral (Ginzburg, Verhaeghe & d'Humières 2008, Sec. 2 and 5.4), so
//! "stable" is a rate at round-off.
//!
//! What it shows with Λ = 3/16. The stream is stable up to a lattice speed proportional to τ⁺ − ½: as
//! τ⁺ → ½ the symmetric modes relax by −(1 − 4(τ⁺ − ½)) and the antisymmetric ones by 1 − (τ⁺ − ½)/Λ, both
//! undamped in the limit, and the stream couples them at a rate proportional to its speed. Along (1, 1, 0)
//! the limit is 6.18 (τ⁺ − ½), a cell Reynolds number |u| Δx/ν of 18.5, on the standing wave
//! k = π (0.336, −0.336, 0.673), which is normal to the stream and changes sign every step; it is the same
//! from τ⁺ = 0.5005 to 0.52. At the laminar floor τ⁺ = 0.505 and the automatic time step's speed 0.05 the
//! stream grows by 6.2e-3 per step. Along a lattice axis the limit is 24.5 (τ⁺ − ½): the stream of #820
//! (U = 0.05 along x) grows for τ⁺ ≤ 0.502 and is neutral from 0.5025; at τ⁺ = 0.5005 it is neutral only
//! up to U = 0.01, or with a near-BGK Λ (τ⁻ ≈ τ⁺).
mod common;

use common::*;
use fieldcraft_lbm::lattice::{C, Q, W, opposite};
use fieldcraft_lbm::setup::TRT_STREAM_LIMIT;
use fieldcraft_lbm::{Limits, Setup};
use fieldcraft_model::{BoxFace, Collision, FaceCondition};
use std::f64::consts::PI;

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

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn c(q: usize) -> [f64; 3] {
    C[q].map(f64::from)
}

/// The collision linearised about ρ = 1 and the velocity `u`: row q, column p is ∂f*_q/∂f_p.
fn collision(tau: f64, magic: f64, u: [f64; 3]) -> [[f64; Q]; Q] {
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
    collide
}

/// The growth rate per step of the wave `k` under the linearised collision `collide`.
fn wave(collide: &[[f64; Q]; Q], k: [f64; 3]) -> f64 {
    let norm = |m: &[Complex]| m.iter().map(|x| x.norm2()).sum::<f64>().sqrt();
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
    log / len
}

/// The largest growth rate per step over the modes of an n³ box (k = 0 left out).
fn growth(tau: f64, magic: f64, u: [f64; 3], n: usize) -> f64 {
    let collide = collision(tau, magic, u);
    let mut worst = f64::NEG_INFINITY;
    for m in 1..n * n * n {
        let k = [m % n, m / n % n, m / (n * n)].map(|i| 2.0 * PI * i as f64 / n as f64);
        worst = worst.max(wave(&collide, k));
    }
    worst
}

const DEFAULT: f64 = 3.0 / 16.0;
/// A rate within this of zero is neutral: the floor of ln ‖G^N‖ / N at N = 2²⁴.
const NEUTRAL: f64 = 1e-6;

/// The stream limit along (1, 1, 0), over τ⁺ − ½, measured in this file's diagonal tests. It is not
/// `TRT_STREAM_LIMIT`: that is the smallest sampled direction's limit (5.6), so these tests keep the
/// (1, 1, 0) limit they were written for.
const DIAGONAL_LIMIT: f64 = 6.18;

/// A stream of lattice speed `speed` along (1, 1, 0).
fn diagonal(speed: f64) -> [f64; 3] {
    [speed / 2f64.sqrt(), speed / 2f64.sqrt(), 0.0]
}

/// Along (1, 1, 0) the stream is stable two per cent below `DIAGONAL_LIMIT` (τ⁺ − ½) and grows two per
/// cent above it, on the standing wave that decides it, from τ⁺ = 0.5005 to 0.52; at the laminar floor
/// and the automatic time step's speed 0.05 its fastest wave grows by 6.2e-3 per step.
#[test]
fn a_diagonal_stream_grows_past_the_stream_limit() {
    let deciding = [0.336 * PI, -0.336 * PI, 0.673 * PI];
    for tau in [0.5005, 0.505, 0.52] {
        let at = |factor: f64| wave(&collision(tau, DEFAULT, diagonal(factor * DIAGONAL_LIMIT * (tau - 0.5))), deciding);
        let (below, above) = (at(0.98), at(1.02));
        println!("τ⁺ {tau}: growth per step {below:+.3e} at 0.98 of the limit, {above:+.3e} at 1.02");
        assert!(below < -NEUTRAL && above > NEUTRAL, "τ⁺ {tau}: {below:e} {above:e}");
    }
    let floor = wave(&collision(0.505, DEFAULT, diagonal(0.05)), [0.3145 * PI, -0.3145 * PI, 0.63 * PI]);
    println!("τ⁺ 0.505, U 0.05 along (1, 1, 0): {floor:+.3e} per step");
    assert!((floor - 6.2e-3).abs() < 1e-4, "{floor:e}");
}

/// No wave of a 16³ box grows two per cent below the limit along (1, 1, 0) at the laminar floor.
#[test]
fn no_wave_of_a_box_grows_below_the_stream_limit() {
    let rate = growth(0.505, DEFAULT, diagonal(0.98 * DIAGONAL_LIMIT * 0.005), 16);
    println!("τ⁺ 0.505, 0.98 of the limit along (1, 1, 0), 16³ waves: largest growth {rate:+.3e} per step");
    assert!(rate <= NEUTRAL, "{rate:e}");
}

/// The smallest limit over the stream's direction is a long-wave onset beside a staggered mode: with
/// k = π e_y + q the growth is κ |q|² − β |q|⁴, and κ changes sign at 5.606 (τ⁺ − ½) for a stream along
/// (1, 0.437, 0.01) (5.603 at τ⁺ = 0.5005, 5.656 at 0.52). `TRT_STREAM_LIMIT` is that, rounded down. Two per
/// cent below it the wave decays and two per cent above it grows; the rates are of order 1e-6 (τ⁺ − ½)/0.005,
/// so they are compared with `WEAK`, not `NEUTRAL`.
#[test]
fn the_smallest_limit_is_a_long_wave_beside_a_staggered_mode() {
    /// About 16 times the error of ln ‖G^N‖ / N at N = 2²⁴ on these waves (6e-9 against the eigenvalues), and a
    /// sixth of the smallest rate asserted.
    const WEAK: f64 = 1e-7;
    let along = |speed: f64| {
        let d = [16.0f64, 7.0, 0.0];
        let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
        d.map(|c| speed * c / len)
    };
    // (τ⁺, the wave over π: π e_y + |q| (−0.68, 0.7332, 0) with |q| = 0.02 π and 0.015 π)
    for (tau, k) in [(0.505, [-0.0136, 1.01466, 0.0]), (0.52, [-0.0102, 1.011, 0.0])] {
        let k = k.map(|c: f64| c * PI);
        let at = |factor: f64| wave(&collision(tau, DEFAULT, along(factor * TRT_STREAM_LIMIT * (tau - 0.5))), k);
        let (below, above) = (at(0.98), at(1.02));
        println!("τ⁺ {tau}: growth per step {below:+.3e} at 0.98 of the limit, {above:+.3e} at 1.02");
        assert!(below < -WEAK && above > WEAK, "τ⁺ {tau}: {below:e} {above:e}");
    }
}

/// A periodic 8³ box (Δx = Δt = 1) in a uniform stream of lattice speed `speed` at `tau`, with the
/// two-relaxation-time parameter `magic` (`None` for the default Λ = 3/16).
fn stream_with(tau: f64, speed: f64, magic: Option<f64>, collision: Collision) -> Setup {
    let mut d = domain([0.0; 3], [8.0; 3], 1.0);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = FaceCondition::Periodic {};
    }
    d.initial_velocity = [speed, 0.0, 0.0];
    Setup::new(&model(1.0, (tau - 0.5) / 3.0, d, collide(10.0, 1.0, magic, collision)), None, &[], &Limits::default()).unwrap()
}

/// A periodic 8³ box in a uniform stream at the default Λ = 3/16.
fn stream(tau: f64, speed: f64, collision: Collision) -> Setup {
    stream_with(tau, speed, None, collision)
}

/// The setup warns where the characteristic lattice speed is past the limit at the fluid's τ⁺, and only there:
/// at the floor with the automatic time step's speed (cell Reynolds number 30), not at τ⁺ = 0.51 (15), not for
/// a slow flow close to ½, not for the regularized collision, and in the band between the smallest sampled
/// direction's limit (`TRT_STREAM_LIMIT`) and the (1, 1, 0) limit, where a stream is unstable in some direction.
#[test]
fn a_stream_past_the_limit_is_warned_about() {
    let warned = |s: &Setup| s.warnings.iter().any(|w| w.contains("linearly unstable") && w.contains("two-relaxation-time collision"));
    let floor = stream(0.505, 0.05, Collision::Trt);
    assert!(warned(&floor), "{:?}", floor.warnings);
    assert!(floor.warnings.iter().any(|w| w.contains("cell Reynolds number 30.0") && w.contains("16.8")), "{:?}", floor.warnings);
    assert!(!warned(&stream(0.51, 0.05, Collision::Trt)));
    assert!(!warned(&stream(0.505, 0.02, Collision::Trt)));
    assert!(warned(&stream(0.52, 0.15, Collision::Trt)));
    assert!(warned(&stream(0.505, 1.02 * TRT_STREAM_LIMIT * 0.005, Collision::Trt)));
    assert!(!warned(&stream(0.505, 0.05, Collision::Regularized)));
}

/// At a Λ other than the default the stream limit is not established: the setup warns once, and does not
/// apply the Λ = 3/16 limit to it. At Λ = 3/16 given explicitly it is the default, and there is no such warning.
#[test]
fn a_non_default_magic_is_warned_as_not_established() {
    let not_established = |s: &Setup| s.warnings.iter().filter(|w| w.contains("has not been established")).count();
    let half = stream_with(0.505, 0.05, Some(0.5), Collision::Trt);
    assert_eq!(not_established(&half), 1, "{:?}", half.warnings);
    assert!(half.warnings.iter().any(|w| w.contains("3.7") && w.contains("4.5")), "{:?}", half.warnings);
    assert!(!half.warnings.iter().any(|w| w.contains("linearly unstable")), "{:?}", half.warnings);
    assert_eq!(not_established(&stream_with(0.505, 0.001, Some(1.0 / 24.0), Collision::Trt)), 1);
    assert_eq!(not_established(&stream_with(0.505, 0.05, Some(3.0 / 16.0), Collision::Trt)), 0);
    assert_eq!(not_established(&stream(0.505, 0.05, Collision::Trt)), 0);
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_trt_uniform_stream_is_linearly_unstable_close_to_one_half_at_the_820_mach_number() {
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
        let rate = growth(tau, magic, [speed, 0.0, 0.0], 16);
        println!("τ⁺ {tau} Λ {magic:e} U {speed}: largest growth {rate:.3e} per step");
        assert_eq!(rate > NEUTRAL, unstable, "τ⁺ {tau} Λ {magic} U {speed}: {rate:e}");
    }
}
