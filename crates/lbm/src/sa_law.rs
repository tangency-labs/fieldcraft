//! The Spalart–Allmaras constant-stress wall law (`wall_model: "spalart_allmaras_bounce"`): the
//! mean velocity u⁺(y⁺) that the Spalart–Allmaras model itself produces in a layer of constant
//! total stress, in closed form, as a wall law beside Spalding's.
//!
//! # Sources
//! - P. R. Spalart, S. R. Allmaras, *A one-equation turbulence model for aerodynamic flows*,
//!   AIAA-92-0439, 1992 (doi:10.2514/6.1992-439); the same model in La Recherche Aérospatiale
//!   1994(1):5–21. Equation numbers are the 1992 paper's. Section "Near-wall region, high
//!   Reynolds number": in the log layer S = u_τ/(κ d) and ν_t = κ u_τ d, and production, diffusion
//!   and destruction balance provided c_w1 = c_b1/κ² + (1 + c_b2)/σ. Section "Near-wall region,
//!   finite Reynolds number": ν̃ is defined so that χ = κ y⁺ from the wall to the log layer, with
//!   f_v1 = χ³/(χ³ + c_v1³) (eq. 7) and f_v2 built so that S̃ = u_τ/(κ y) all the way to the wall
//!   (eq. 8); the transport equation (eq. 9) is then in equilibrium all the way to d = 0. The
//!   constants κ 0.41 and c_v1 7.1 (Appendix, with eqs. A1–A2).
//! - R. H. Nichols, C. C. Nelson, *Wall function boundary conditions including heat transfer and
//!   compressibility*, AIAA J. 42(6):1107–1114, 2004. The idea of a wall function that matches
//!   the first node's ν̃ to the model's own wall behaviour rather than to another model's.
//! - The closed form is the integral of a rational function by partial fractions (any calculus
//!   text). Its quartic is factored into two real quadratics by Ferrari's method, the resolvent
//!   cubic is solved by Viète's trigonometric formula (its three roots are real).
//!
//! # The law
//! In a layer of constant total stress u_τ² (wall units: ν = u_τ = 1), the shear is
//! du⁺/dy⁺ = 1 / (1 + ν_t⁺) with ν_t⁺ = χ f_v1(χ). The model's steady, one-dimensional equation
//!
//! ```text
//! 0 = c_b1 (1 − f_t2) S̃ χ − (c_w1 f_w − c_b1 f_t2 / κ²) (χ / y⁺)² + (1/σ) [((1 + χ) χ')' + c_b2 χ'²]
//! ```
//!
//! is solved at **every** height down to the wall by χ = κ y⁺: S̃ is then 1 / (κ y⁺) because
//! 1 − f_v2 = χ / (1 + χ f_v1), so r = 1 and f_w = 1; the f_t2 terms cancel (they enter P and D
//! with opposite signs and the same factor); and `c_w1 = c_b1 / κ² + (1 + c_b2) / σ` balances the
//! rest. So the wall law that goes with the model is the integral, with χ = κ y and
//! x = κ y, and the rational integrand (x³ + a) / (x⁴ + x³ + a), a = c_v1³:
//!
//! ```text
//! u⁺(y⁺) = ∫₀^y⁺ dy / (1 + κ y f_v1(κ y)) = (1/κ) F(κ y⁺),    F(X) = ∫₀^X (x³ + a) / (x⁴ + x³ + a) dx
//! ```
//!
//! with nothing fitted: the same κ as the model. The factoring x⁴ + x³ + a = (x² + p₁x + q₁)(x² +
//! p₂x + q₂) has p₁ + p₂ = 1, q₁q₂ = a and the partial fractions
//! (A₁x + B₁)/(x² + p₁x + q₁) + (A₂x + B₂)/(x² + p₂x + q₂), each integrated to
//! `½A log(x² + px + q) + ((2B − Ap)/d) atan((2x + p)/d)` with d = √(4q − p²) > 0.
//!
//! **Log-law intercept.** As X → ∞, F(X) = ln X + C + S(X) with S(X) = −∫_X^∞ (f − 1/x) dx, where
//! f − 1/x = −(x³ − a x + a) / (x (x⁴ + x³ + a)) ≤ 0 for x ≥ √a. So 0 ≤ S(X) ≤ ∫_X^∞ dx/x² = 1/X
//! (x³ − a x + a lies between 0 and (x⁴ + x³ + a)/x there). Hence u⁺ approaches ln(y⁺)/κ + B∞ from
//! above, with B∞ = (C + ln κ)/κ about 5.03, and the gap is at most 1/(κ² y⁺).
//!
//! The law differs from Spalding's (κ 0.41, B 5.5) in the buffer layer: at y⁺ 15 it gives 11.09
//! where Spalding's gives 10.40. A wall function that imposes ν̃ at the
//! first node and carries a stress at it should use the pair the model is in equilibrium with.
//!
//! **Domain.** A height whose κ y⁺ overflows the float range is not usable: the functions return
//! 0 for it, as for every non-finite or non-positive input.
//!
//! [`friction_velocity`] inverts the law: the friction velocity of a flow of speed U at distance y
//! from the wall with viscosity ν solves u_τ u⁺(y u_τ / ν) = U, i.e. y⁺ u⁺(y⁺) = U y / ν, whose
//! left side rises with y⁺. [`y_plus_for_reynolds`] finds the root by bisection at the geometric mean of a bracket, whose step count is bounded.

use crate::sa::{C_V1, KAPPA};
use std::f64::consts::FRAC_PI_2;

/// The factored law: p, q, A and B of the two quadratic factors of x⁴ + x³ + c_v1³ (index 0 is
/// factor 1, index 1 is factor 2), with the partial fractions (A x + B) of the integrand.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Factors {
    pub p: [f64; 2],
    pub q: [f64; 2],
    pub a: [f64; 2],
    pub b: [f64; 2],
}

/// Bits of the inversion's starting bracket in log space: ln(hi/lo) ≤ ln √re ≤ ln √f64::MAX = 354.9
/// < 2⁹ (u⁺(√re) ≥ 1 once √re ≥ 1, and hi/lo is about 1 below it; see [`y_plus_for_reynolds`]).
const LOG_RATIO_BITS: u32 = 9;

/// The bisections of the inversion: each halves ln(hi/lo), which must fall from below 2⁹ to the
/// stop width ε (relative), so 9 + 52 halvings suffice; one spare.
const INVERSION_STEPS: u32 = f64::MANTISSA_DIGITS + LOG_RATIO_BITS;

/// The factored law of x⁴ + x³ + c_v1³ (see the module docs), from the constant c_v1 only.
pub fn factors() -> Factors {
    let a = C_V1 * C_V1 * C_V1;
    // Resolvent cubic v³ − 64 a v − 64 a = 0 of the factoring; v = (p₁ − p₂)² − 1 with all three
    // roots real. Viète's trigonometric form gives the largest root.
    let k = 64.0 * a;
    let v = 2.0 * (k / 3.0).sqrt() * libm::cos(libm::acos(1.5 * (3.0 / k).sqrt()) / 3.0);
    let m = (v + 1.0).sqrt() / 2.0;
    let (p1, p2) = (0.5 + m, 0.5 - m);
    let s = m * m - 0.25;
    let q1 = 0.5 * s * (1.0 + 0.5 / m);
    let q2 = 0.5 * s * (1.0 - 0.5 / m);
    // Partial fractions of (x³ + a) / ((x² + p₁x + q₁)(x² + p₂x + q₂)), δ = q₂ − q₁.
    let delta = q2 - q1;
    let a1 = (p1 * p1 - q1 + 2.0 * m * (a + p1 * q1) / delta) / (delta + 2.0 * m * p1 + 4.0 * m * m * q1 / delta);
    let a2 = 1.0 - a1;
    let b1 = (a + p1 * q1 - 2.0 * m * q1 * a1) / delta;
    let b2 = -p1 + 2.0 * m * a1 - b1;
    Factors { p: [p1, p2], q: [q1, q2], a: [a1, a2], b: [b1, b2] }
}

/// The constant C of F(X) = ln X + C + S(X), S → 0: Σᵢ [ −½ Aᵢ ln qᵢ + ((2Bᵢ − Aᵢpᵢ)/dᵢ)(π/2 − atan(pᵢ/dᵢ)) ].
fn log_law_constant(f: &Factors) -> f64 {
    let mut c = 0.0;
    for (((p, q), a), b) in f.p.iter().zip(&f.q).zip(&f.a).zip(&f.b) {
        let d = (4.0 * q - p * p).sqrt();
        c += -0.5 * a * libm::log(*q) + ((2.0 * b - a * p) / d) * (FRAC_PI_2 - libm::atan(p / d));
    }
    c
}

/// B∞ = (C + ln κ) / κ: u⁺ − ln(y⁺)/κ tends to this from above as y⁺ grows (about 5.03).
pub fn log_law_intercept() -> f64 {
    (log_law_constant(&factors()) + libm::log(KAPPA)) / KAPPA
}

/// ln((x² + p x + q) / q) for a finite x ≥ 0, without the cancellation of two logarithms near
/// x = 0: below 1 it is log1p((x² + p x) / q); above 1 it is 2 ln x + log1p((p + q/x)/x) − ln q,
/// which never forms x² and so cannot overflow.
fn log_ratio(x: f64, p: f64, q: f64) -> f64 {
    if x < 1.0 { libm::log1p((x * x + p * x) / q) } else { 2.0 * libm::log(x) + libm::log1p((p + q / x) / x) - libm::log(q) }
}

/// One factor's integral from 0 to x: `½A ln((x² + px + q)/q) + ((2B − Ap)/d)(atan((2x + p)/d) − atan(p/d))`.
/// The difference of the two arctangents is atan2(2x/d, 1 + p(2x + p)/d²), the same angle with no
/// cancellation at small x. For x ≥ 0 the difference lies in [0, π) (atan rises and spans less than
/// π), and atan2 with a first argument 2x/d ≥ 0 returns exactly that range.
fn factor_integral(x: f64, p: f64, q: f64, a: f64, b: f64) -> f64 {
    let d = (4.0 * q - p * p).sqrt();
    let arctan_difference = libm::atan2(2.0 * x / d, 1.0 + p * (2.0 * x + p) / (d * d));
    0.5 * a * log_ratio(x, p, q) + ((2.0 * b - a * p) / d) * arctan_difference
}

/// F(x) = ∫₀^x (t³ + a) / (t⁴ + t³ + a) dt for a finite x ≥ 0.
fn integral(x: f64) -> f64 {
    let f = factors();
    factor_integral(x, f.p[0], f.q[0], f.a[0], f.b[0]) + factor_integral(x, f.p[1], f.q[1], f.a[1], f.b[1])
}

/// ν_t⁺ = κ y⁺ f_v1(κ y⁺) = x⁴ / (x³ + c_v1³) with x = κ y⁺, written x / (1 + c_v1³ / x³) so
/// that it does not form x⁴ (it stays finite up to the float range of x). 0 at y⁺ ≤ 0 or for a
/// height that is not usable (see the module docs).
pub fn eddy_viscosity_plus(y_plus: f64) -> f64 {
    let x = KAPPA * y_plus;
    if !(y_plus > 0.0 && x.is_finite()) {
        return 0.0;
    }
    let nt = x / (1.0 + C_V1 * C_V1 * C_V1 / (x * x * x));
    if nt.is_finite() { nt } else { 0.0 }
}

/// u⁺(y⁺): the mean velocity in wall units of the equilibrium layer at height `y_plus`. 0 at
/// y⁺ ≤ 0 or for a height that is not usable.
pub fn u_plus(y_plus: f64) -> f64 {
    let x = KAPPA * y_plus;
    if !(y_plus > 0.0 && x.is_finite()) {
        return 0.0;
    }
    let u = integral(x) / KAPPA;
    if u.is_finite() && u >= 0.0 { u } else { 0.0 }
}

/// The y⁺ with y⁺ u⁺(y⁺) = `re`, the inverse of the law's product (a Reynolds number U y / ν).
/// 0 for `re` ≤ 0, for a `re` that is not a number, and where the root is not usable.
///
/// The root is bracketed by lo = √re (u⁺ ≤ y⁺ gives y·u⁺ ≤ y²) and hi = re / u⁺(lo) (u⁺ rises, so
/// y·u⁺ ≥ hi·u⁺(lo) = re). Each step bisects at the geometric mean, halving ln(hi/lo), until the
/// bracket is ε·hi wide: at most [`INVERSION_STEPS`] steps, which the bound on the starting bracket
/// guarantees, so the loop never stops early.
pub fn y_plus_for_reynolds(re: f64) -> f64 {
    if !(re > 0.0 && re.is_finite()) {
        return 0.0;
    }
    let (mut lo, u_lo) = (re.sqrt(), u_plus(re.sqrt()));
    if !(u_lo > 0.0) {
        return 0.0;
    }
    let mut hi = re / u_lo;
    // A bracket end whose κy⁺ overflows is not usable (see the module docs).
    if !(hi.is_finite() && (KAPPA * hi).is_finite()) {
        return 0.0;
    }
    // Bisection at the geometric mean halves ln(hi/lo) each step: a bound on the steps that holds
    // whatever the function's shape ([`INVERSION_STEPS`]).
    for _ in 0..INVERSION_STEPS {
        if hi - lo <= f64::EPSILON * hi {
            break;
        }
        let mid = lo * (hi / lo).sqrt();
        let value = mid * u_plus(mid) - re;
        if value == 0.0 {
            return mid;
        }
        if value > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let y = lo * (hi / lo).sqrt();
    if y.is_finite() && y > 0.0 { y } else { 0.0 }
}

/// The friction velocity of a flow of speed `speed` parallel to a wall at distance `y` with
/// kinematic viscosity `nu`, from the equilibrium layer: u_τ with u_τ u⁺(y u_τ / ν) = U. 0 for a
/// speed, distance or viscosity that is not a positive number, or where no usable root exists.
pub fn friction_velocity(speed: f64, y: f64, nu: f64) -> f64 {
    if !(speed > 0.0 && speed.is_finite() && y > 0.0 && y.is_finite() && nu > 0.0 && nu.is_finite()) {
        return 0.0;
    }
    let re = speed * y / nu;
    let y_plus = y_plus_for_reynolds(re);
    let ut = nu * y_plus / y;
    if ut.is_finite() { ut } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa::{C_B1, C_B2, R_MAX, SIGMA, SaField, SaInput, Side, c_w1, f_t2, f_w, s_tilde};

    const EPS: f64 = f64::EPSILON;

    /// The equilibrium layer's shear du⁺/dy⁺ = 1 / (1 + ν_t⁺).
    fn shear(y: f64) -> f64 {
        1.0 / (1.0 + eddy_viscosity_plus(y))
    }

    /// 5-point Gauss–Legendre nodes and weights on [−1, 1] (Abramowitz & Stegun, Table 25.4).
    const GL_NODES: [f64; 5] = [0.0, 0.538_469_310_105_683_1, -0.538_469_310_105_683_1, 0.906_179_845_938_664, -0.906_179_845_938_664];
    const GL_WEIGHTS: [f64; 5] =
        [0.568_888_888_888_888_9, 0.478_628_670_499_366_5, 0.478_628_670_499_366_5, 0.236_926_885_056_189_1, 0.236_926_885_056_189_1];

    /// ∫ u⁺' dy from y⁺ 10⁻³ to `top` by composite Gauss–Legendre panels in ln y, `per_decade`
    /// panels per decade, plus the head [0, 10⁻³] from the sublayer series. Independent of the
    /// factoring: it uses only the shear.
    fn quadrature(top: f64, per_decade: f64) -> f64 {
        let a = C_V1 * C_V1 * C_V1;
        let (lo, hi) = (1e-3f64.ln(), top.ln());
        let panels = ((hi - lo) / std::f64::consts::LN_10 * per_decade).ceil().max(1.0) as usize;
        let h = (hi - lo) / panels as f64;
        let mut sum = 1e-3 - KAPPA.powi(4) * 1e-3f64.powi(5) / (5.0 * a);
        for p in 0..panels {
            let mid = lo + (p as f64 + 0.5) * h;
            for (x, w) in GL_NODES.iter().zip(GL_WEIGHTS) {
                let y = (mid + 0.5 * h * x).exp();
                sum += 0.5 * h * w * y * shear(y);
            }
        }
        sum
    }

    /// The sublayer series u⁺ = y⁺ − κ⁴ y⁺⁵ / (5 c_v1³) + O(y⁺⁸) and the log law, with the
    /// monotone bound u⁺ ≤ y⁺ and the intercept of the model's own log law.
    #[test]
    fn the_sublayer_and_the_log_law() {
        let a = C_V1 * C_V1 * C_V1;
        // The next series term is κ⁷ y⁸ / (8 a²) < 1e-20 here; the closed form's round-off is a
        // few ulps of terms of size ≲ 4, so 64 ε absolute holds.
        for y in [1e-4f64, 1e-3, 1e-2, 0.1] {
            let series = y - KAPPA.powi(4) * y.powi(5) / (5.0 * a);
            assert!((u_plus(y) - series).abs() <= 64.0 * EPS, "y⁺ {y}: {} vs {series}", u_plus(y));
        }
        assert_eq!(u_plus(0.0), 0.0);
        // u⁺ ≤ y⁺ since the integrand is at most 1; strictly increasing.
        let mut last = 0.0;
        for k in 0..=2000 {
            let y = 10f64.powf(-3.0 + 10.0 * k as f64 / 2000.0);
            let u = u_plus(y);
            assert!(u > last && u <= y + 64.0 * EPS * y, "y⁺ {y}: {u}");
            last = u;
        }
        // The intercept: u⁺ − ln(y⁺)/κ → B∞ from above, with the remainder S(X)/κ bounded in the module
        // docs: −1/(κ² y⁺) ≤ ln(y⁺)/κ + B∞ − u⁺ ≤ 0 for y⁺ ≥ √a/κ (round-off of ln y/κ ≈ 34 allows 1e-12).
        let b_inf = log_law_intercept();
        println!("B_inf = {b_inf}");
        for y in [a.sqrt() / KAPPA, 100.0, 1e3, 1e4, 1e5, 1e6] {
            let gap = y.ln() / KAPPA + b_inf - u_plus(y);
            assert!(-1.0 / (KAPPA * KAPPA * y) - 1e-12 <= gap && gap <= 1e-12, "y⁺ {y}: gap {gap}");
        }
    }

    /// The closed form against an independent quadrature of the same integral, at two panel
    /// densities. Gauss–Legendre's error falls by 2¹⁰ when the panels halve, so the coarser
    /// difference bounds the finer one's error; the closed form must agree within it.
    #[test]
    fn the_closed_form_matches_an_independent_quadrature() {
        for y in [0.01, 0.5, 3.0, 7.5, 15.0, 30.0, 60.0, 180.0, 1000.0, 5200.0, 1e5, 1e7 - 1.0] {
            let q32 = quadrature(y, 32.0);
            let q64 = quadrature(y, 64.0);
            let u = u_plus(y);
            assert!((u - q64).abs() <= (q32 - q64).abs() + 64.0 * EPS * u, "y⁺ {y}: closed form {u}, quadrature {q64}");
        }
    }

    /// χ = κ y⁺ solves the model's one-dimensional steady equation at every height, buffer layer
    /// included, with the crate's own S̃, f_w, f_t2 and constants. The terms are O(1) and cancel
    /// analytically, so the residual is bounded by round-off of the largest term.
    #[test]
    fn the_law_solves_the_models_equation() {
        for k in 0..=120 {
            let y = 0.05 * 10f64.powf(5.0 * k as f64 / 120.0);
            let chi = KAPPA * y;
            let s = s_tilde(shear(y), chi, 1.0, y);
            assert!((s * chi - 1.0).abs() <= 1e-12, "y⁺ {y}: S̃ κ y = {}", s * chi);
            let r = (chi / (s * KAPPA * KAPPA * y * y)).min(R_MAX);
            assert!((r - 1.0).abs() <= 1e-12, "y⁺ {y}: r {r}");
            let ft2 = f_t2(chi);
            let production = C_B1 * (1.0 - ft2) * s * chi;
            let destruction = (c_w1() * f_w(r) - C_B1 * ft2 / (KAPPA * KAPPA)) * (chi / y) * (chi / y);
            let diffusion = (1.0 + C_B2) * KAPPA * KAPPA / SIGMA;
            let residual = production - destruction + diffusion;
            // The production crosses zero where f_t2 = 1, so the bar is on the largest term.
            let scale = production.abs().max(destruction.abs()).max(diffusion.abs());
            assert!(residual.abs() <= 1e-11 * scale, "y⁺ {y}: residual {residual} against {scale}");
        }
    }

    /// The lattice transport is exact on the law's ν̃: on SaField's discrete operator a linear
    /// ν̃ = s y gives every face slope s (to a fluid cell, to the on-node wall's node one cell away,
    /// to an inflow face half a cell away), every face diffusivity ν + ν̃ at the face, and the c_b2
    /// term s², so one sub-step of transport alone adds exactly (1 + c_b2) s² / σ to every cell. With
    /// [`the_law_solves_the_models_equation`] (the sources take that much off) this is why the law is
    /// the steady state of the model at the wall-law node. (The model's update splits transport and
    /// sources, so its own fixed point carries a splitting error of order Δt; `sa_channel.rs`.)
    /// A small u_τ keeps the coefficient sum below 1, so one sub-step is the whole lattice step.
    #[test]
    fn the_lattice_transport_is_exact_on_the_law() -> crate::Result<()> {
        // u_τ 0.01 (lattice), h = 0.75 wall units per cell: ν = u_τ / h, ν̃ = κ u_τ y.
        let (n, h, ut) = (28usize, 0.75, 0.01);
        let (nu, s) = (ut / h, KAPPA * ut);
        let sides: Vec<[Side; 6]> = (0..n)
            .map(|j| {
                let below = if j == 0 { Side::Node } else { Side::Cell((j - 1) as u32) };
                let above = if j + 1 == n { Side::Inflow } else { Side::Cell((j + 1) as u32) };
                [Side::Open, Side::Open, below, above, Side::Open, Side::Open]
            })
            .collect();
        let mut field = SaField::new(sides, vec![true; n], s * (n as f64 + 0.5), 0.0)?;
        let y: Vec<f64> = (0..n).map(|j| (j + 1) as f64).collect();
        let law: Vec<f64> = y.iter().map(|y| s * y).collect();
        field.set_nu_tilde(&law)?;
        let u = vec![[0.0; 3]; n];
        let vorticity = vec![0.0; n];
        let input = SaInput { u: &u, vorticity: &vorticity, length: &y, nu, model: false, source: None, substeps: Some(1) };
        field.advance(&input)?;
        let gain = (1.0 + C_B2) * s * s / SIGMA;
        for (j, before) in law.iter().enumerate() {
            let got = field.nu_tilde().get(j).copied().unwrap_or(f64::NAN) - before;
            // One rounding of ν̃ + Δt·rate (ε ν̃) and the rate's sum of at most 16 terms, each at
            // most 2 (ν + ν̃ + s) s / σ and rounded once: 64 ε times their total.
            let bar = 64.0 * f64::EPSILON * (before + 32.0 * (nu + before + s) * s / SIGMA);
            assert!((got - gain).abs() <= bar, "cell {j}: gain {got} vs {gain} (bar {bar})");
        }
        Ok(())
    }

    /// The inversion: y⁺ with y⁺ u⁺(y⁺) = Re is bracketed to the round-off of the stop width, both
    /// for the Reynolds number of each height and for the friction velocity of each speed.
    #[test]
    fn the_inversion_round_trips() {
        for k in 0..=4000 {
            let y = 10f64.powf(-3.5 + 10.4 * k as f64 / 4000.0);
            let re = y * u_plus(y);
            let g = |t: f64| t * u_plus(t) - re;
            let back = y_plus_for_reynolds(re);
            // g rises with t and its log-derivative is at least 1, so g changes sign across the
            // round-off band of ±64 ε around the root.
            assert!(g(back * (1.0 - 64.0 * EPS)) <= 0.0 && 0.0 <= g(back * (1.0 + 64.0 * EPS)), "y⁺ {y}: back {back}");
        }
        let (y, nu) = (1.0, 1.57e-3);
        let mut last = 0.0;
        for k in 1..=200 {
            let speed = 1e-4 * 1.07f64.powi(k);
            let ut = friction_velocity(speed, y, nu);
            assert!(ut > last, "speed {speed}: u_τ {ut} not above {last}");
            last = ut;
            // h(t) = t u⁺(y t / ν) − U changes sign across the same ±64 ε band around u_τ.
            let h = |t: f64| t * u_plus(y * t / nu) - speed;
            assert!(h(ut * (1.0 - 64.0 * EPS)) <= 0.0 && 0.0 <= h(ut * (1.0 + 64.0 * EPS)), "speed {speed}: u_τ {ut}");
        }
    }

    /// Hostile input: zero, negative, infinite and not-a-number arguments give 0 or a finite
    /// non-negative value, never a panic or a non-finite number.
    #[test]
    fn hostile_input_gives_zero() {
        let bad = [0.0, -1.0, -0.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY, f64::MIN_POSITIVE, f64::MAX, 1e-300, 1e300];
        for x in bad {
            for v in [u_plus(x), y_plus_for_reynolds(x), eddy_viscosity_plus(x)] {
                assert!(v.is_finite() && v >= 0.0, "{x}: {v}");
            }
            for (speed, y, nu) in [(x, 1.0, 1e-3), (0.1, x, 1e-3), (0.1, 1.0, x)] {
                let ut = friction_velocity(speed, y, nu);
                assert!(ut.is_finite() && ut >= 0.0, "{speed} {y} {nu}: {ut}");
            }
        }
        assert_eq!(friction_velocity(f64::NAN, 1.0, 1e-3), 0.0);
        assert_eq!(friction_velocity(0.0, 1.0, 1e-3), 0.0);
        assert_eq!(friction_velocity(0.1, 1.0, 0.0), 0.0);
    }
}
