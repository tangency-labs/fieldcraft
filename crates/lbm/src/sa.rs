//! Spalart–Allmaras one-equation turbulence model (`turbulence: {"type": "spalart_allmaras"}`),
//! in its negative form (SA-neg).
//!
//! # Sources
//! - P. R. Spalart, S. R. Allmaras, *A one-equation turbulence model for aerodynamic flows*,
//!   La Recherche Aérospatiale 1994(1):5–21 (first as AIAA-92-0439, 1992). The transport equation,
//!   its constants, the near-wall functions f_v1, f_v2, f_w and the trip-suppression term f_t2.
//! - S. R. Allmaras, F. T. Johnson, P. R. Spalart, *Modifications and clarifications for the
//!   implementation of the Spalart–Allmaras turbulence model*, ICCFD7-1902, 2012. The limiter on
//!   the modified vorticity S̃ (their eqs. 11–12) and the negative branch, SA-neg (§2.4,
//!   eqs. 20–23), which keeps the equation well behaved where the discrete ν̃ goes below zero.
//!
//! # The model
//! The working variable ν̃ obeys
//!
//! ```text
//! ∂ν̃/∂t + u·∇ν̃ = P − D + (1/σ) [∇·((ν + ν̃) ∇ν̃) + c_b2 |∇ν̃|²]
//! P = c_b1 (1 − f_t2) S̃ ν̃,        D = (c_w1 f_w − c_b1 f_t2 / κ²) (ν̃ / d)²
//! ν_t = ν̃ f_v1,   χ = ν̃ / ν,   f_v1 = χ³ / (χ³ + c_v1³),   f_v2 = 1 − χ / (1 + χ f_v1)
//! S̃ = Ω + S̄,   S̄ = ν̃ f_v2 / (κ² d²),   f_t2 = c_t3 exp(−c_t4 χ²)
//! r = min(ν̃ / (S̃ κ² d²), 10),   g = r + c_w2 (r⁶ − r),   f_w = g ((1 + c_w3⁶) / (g⁶ + c_w3⁶))^(1/6)
//! ```
//!
//! with Ω the vorticity magnitude and d the distance to the nearest wall. Allmaras et al.
//! keep S̃ positive where S̄ would take it below 0.3 Ω:
//! `S̃ = Ω + Ω (c_v2² Ω + c_v3 S̄) / ((c_v3 − 2 c_v2) Ω − S̄)` when `S̄ < −c_v2 Ω`.
//!
//! **SA-neg.** Where ν̃ < 0 the eddy viscosity is 0, and the equation becomes
//!
//! ```text
//! P = c_b1 (1 − c_t3) Ω ν̃,   D = −c_w1 (ν̃ / d)²,   diffusion with ν + ν̃ f_n,
//! f_n = (c_n1 + χ³) / (c_n1 − χ³)
//! ```
//!
//! Both sources then drive ν̃ back towards 0, and `ν + ν̃ f_n` stays positive for every χ < 0
//! (its minimum over χ is about 0.0094 ν, at χ ≈ −1.37), so the diffusion never turns anti-diffusive.
//!
//! **Log layer.** ν̃ = κ u_τ y with the constant-stress shear `Ω = u_τ² / (ν + ν_t)` gives
//! `S̃ = u_τ / (κ y)`, r = 1, f_w = 1, and the balance `P − D + diffusion = 0` holds exactly once
//! f_t2 has died out (χ above about 10), since `c_w1 = c_b1 / κ² + (1 + c_b2) / σ` (the paper's §3).
//!
//! No code of any CFD package was read; this is written from the papers.

/// c_b1, the production constant.
pub const C_B1: f64 = 0.1355;
/// σ, the diffusion Prandtl number.
pub const SIGMA: f64 = 2.0 / 3.0;
/// c_b2, the gradient-squared diffusion constant.
pub const C_B2: f64 = 0.622;
/// κ, von Kármán's constant: 0.41, the model's own (Spalart & Allmaras 1994, the final form's
/// constants; c_w1 is built from it).
pub const KAPPA: f64 = 0.41;
pub const C_W2: f64 = 0.3;
pub const C_W3: f64 = 2.0;
pub const C_V1: f64 = 7.1;
pub const C_T3: f64 = 1.2;
pub const C_T4: f64 = 0.5;
/// c_v2 of the S̃ limiter (Allmaras et al. 2012).
pub const C_V2: f64 = 0.7;
/// c_v3 of the S̃ limiter (Allmaras et al. 2012).
pub const C_V3: f64 = 0.9;
/// c_n1 of SA-neg's diffusion (Allmaras et al. 2012).
pub const C_N1: f64 = 16.0;
/// The cap on r.
pub const R_MAX: f64 = 10.0;

#[inline(always)]
pub fn c_w1() -> f64 {
    C_B1 / (KAPPA * KAPPA) + (1.0 + C_B2) / SIGMA
}

#[inline(always)]
pub fn f_v1(chi: f64) -> f64 {
    let c3 = chi * chi * chi;
    c3 / (c3 + C_V1 * C_V1 * C_V1)
}

#[inline(always)]
pub fn f_v2(chi: f64) -> f64 {
    1.0 - chi / (1.0 + chi * f_v1(chi))
}

#[inline(always)]
pub fn f_t2(chi: f64) -> f64 {
    C_T3 * libm::exp(-C_T4 * chi * chi)
}

/// f_n(χ) = (c_n1 + χ³) / (c_n1 − χ³), SA-neg's diffusion factor (χ < 0 only; 1 otherwise).
#[inline(always)]
pub fn f_n(chi: f64) -> f64 {
    if chi >= 0.0 {
        return 1.0;
    }
    let c3 = chi * chi * chi;
    (C_N1 + c3) / (C_N1 - c3)
}

/// f_w(r) with `g = r + c_w2 (r⁶ − r)`.
#[inline(always)]
pub fn f_w(r: f64) -> f64 {
    let g = r + C_W2 * (r.powi(6) - r);
    let c6 = C_W3.powi(6);
    g * libm::pow((1.0 + c6) / (g.powi(6) + c6), 1.0 / 6.0)
}

/// The eddy viscosity `ν̃ f_v1(ν̃/ν)`, 0 where ν̃ is not above 0 (SA-neg) or a number is not
/// usable. Never negative, always finite.
#[inline(always)]
pub fn nu_t(nu_tilde: f64, nu: f64) -> f64 {
    if !(nu_tilde > 0.0 && nu > 0.0) {
        return 0.0;
    }
    let v = nu_tilde * f_v1(nu_tilde / nu);
    if v.is_finite() && v > 0.0 { v } else { 0.0 }
}

/// The diffusivity of ν̃ before the 1/σ: `ν + ν̃` for ν̃ ≥ 0, `ν + ν̃ f_n` below (positive for every
/// ν̃ when ν > 0).
#[inline(always)]
pub fn diffusivity(nu_tilde: f64, nu: f64) -> f64 {
    if nu_tilde >= 0.0 { nu + nu_tilde } else { nu + nu_tilde * f_n(nu_tilde / nu) }
}

/// The modified vorticity S̃ with Allmaras et al.'s limiter, for vorticity magnitude `omega`
/// and wall distance `d` (both positive).
#[inline(always)]
pub fn s_tilde(omega: f64, nu_tilde: f64, nu: f64, d: f64) -> f64 {
    let s_bar = nu_tilde * f_v2(nu_tilde / nu) / (KAPPA * KAPPA * d * d);
    if s_bar >= -C_V2 * omega { omega + s_bar } else { omega + omega * (C_V2 * C_V2 * omega + C_V3 * s_bar) / ((C_V3 - 2.0 * C_V2) * omega - s_bar) }
}

/// The source terms of a cell, split for a point-implicit update: the rate of ν̃ is
/// `production − destruction`. Both 0 where an input is not usable or `d` is not above 0.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Sources {
    pub production: f64,
    pub destruction: f64,
}

/// P and D of the cell (see the module docs) for ν̃ `nu_tilde`, molecular viscosity `nu`,
/// vorticity magnitude `omega` and length `d` (the wall distance, or DDES's d̃).
pub fn sources(nu_tilde: f64, nu: f64, omega: f64, d: f64) -> Sources {
    if !(nu > 0.0 && d > 0.0 && omega >= 0.0 && nu_tilde.is_finite() && omega.is_finite() && !d.is_nan()) {
        return Sources::default();
    }
    let chi = nu_tilde / nu;
    let ratio = nu_tilde / d;
    let out = if nu_tilde >= 0.0 {
        let ft2 = f_t2(chi);
        let st = s_tilde(omega, nu_tilde, nu, d);
        let denominator = st * KAPPA * KAPPA * d * d;
        let r = if denominator > 0.0 { (nu_tilde / denominator).min(R_MAX) } else { R_MAX };
        Sources { production: C_B1 * (1.0 - ft2) * st * nu_tilde, destruction: (c_w1() * f_w(r) - C_B1 / (KAPPA * KAPPA) * ft2) * ratio * ratio }
    } else {
        Sources { production: C_B1 * (1.0 - C_T3) * omega * nu_tilde, destruction: -c_w1() * ratio * ratio }
    };
    if out.production.is_finite() && out.destruction.is_finite() { out } else { Sources::default() }
}

/// The vorticity magnitude `Ω = √(2 W:W)` of a velocity gradient `g[i][j] = ∂_i u_j`.
#[inline(always)]
pub fn vorticity(g: &[[f64; 3]; 3]) -> f64 {
    let wx = g[1][2] - g[2][1];
    let wy = g[2][0] - g[0][2];
    let wz = g[0][1] - g[1][0];
    let w = (wx * wx + wy * wy + wz * wz).sqrt();
    if w.is_finite() { w } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The log layer: ν̃ = κ u_τ y, Ω = u_τ² / (ν + ν_t) (constant total stress). Then S̃ = u_τ/(κy),
    /// r = f_w = 1 and P − D + (1/σ)[d/dy((ν+ν̃) dν̃/dy) + c_b2 (dν̃/dy)²] = 0, with the diffusion
    /// (1 + c_b2) κ² u_τ² / σ for the linear ν̃. Holds to round-off for y⁺ 50 to 500 (f_t2 is
    /// below 1e-100 there).
    #[test]
    fn the_log_layer_balances() {
        let (nu, u_tau) = (1.5e-5, 0.8);
        for k in 0..=40 {
            let y_plus = 50.0 * 10f64.powf(k as f64 / 40.0);
            let y = y_plus * nu / u_tau;
            let nt = KAPPA * u_tau * y;
            let omega = u_tau * u_tau / (nu + nu_t(nt, nu));
            assert!((s_tilde(omega, nt, nu, y) * KAPPA * y / u_tau - 1.0).abs() < 1e-13, "y⁺ {y_plus}");
            let s = sources(nt, nu, omega, y);
            let diffusion = (1.0 + C_B2) * KAPPA * KAPPA * u_tau * u_tau / SIGMA;
            let residual = s.production - s.destruction + diffusion;
            assert!(residual.abs() <= 1e-12 * s.production, "y⁺ {y_plus}: {residual} vs P {}", s.production);
        }
    }

    #[test]
    fn f_v1_is_its_closed_form() {
        for chi in [0.0f64, 0.1, 1.0, 3.0, 7.1, 20.0, 1e3] {
            let want = chi.powi(3) / (chi.powi(3) + 357.911);
            assert!((f_v1(chi) - want).abs() <= 1e-15, "{chi}");
        }
        assert!((c_w1() - (0.1355 / 0.1681 + 1.622 * 1.5)).abs() < 1e-15);
    }

    /// The Allmaras limiter: S̃ continuous at S̄ = −c_v2 Ω and positive below it.
    #[test]
    fn the_limiter_is_continuous_and_keeps_s_tilde_positive() {
        let (omega, nu, d) = (2.0, 1.0, 1.0);
        // f_v2 < 0 for χ between about 1 and a few: find ν̃ with S̄ just above and below −c_v2 Ω.
        let s_bar = |nt: f64| nt * f_v2(nt / nu) / (KAPPA * KAPPA * d * d);
        let mut nt_cross = None;
        let mut prev = 0.0;
        for k in 1..100_000 {
            let nt = k as f64 * 1e-4;
            if s_bar(nt) < -C_V2 * omega && s_bar(prev) >= -C_V2 * omega {
                nt_cross = Some((prev, nt));
                break;
            }
            prev = nt;
        }
        let (a, b) = nt_cross.unwrap_or((0.0, 0.0));
        assert!(b > 0.0, "S̄ never went below −c_v2 Ω");
        assert!((s_tilde(omega, a, nu, d) - s_tilde(omega, b, nu, d)).abs() < 1e-3);
        for k in 0..2000 {
            let nt = k as f64 * 5e-3;
            assert!(s_tilde(omega, nt, nu, d) > 0.0, "{nt}");
        }
    }

    /// SA-neg: both sources drive a negative ν̃ to 0 (rate > 0), and the diffusivity stays positive.
    #[test]
    fn sa_neg_pushes_back_to_zero() {
        let nu = 1e-3;
        for k in 1..400 {
            let nt = -(k as f64) * 1e-5;
            let s = sources(nt, nu, 0.05, 0.7);
            assert!(s.production - s.destruction > 0.0, "{nt}");
            assert!(diffusivity(nt, nu) > 0.009 * nu, "{nt}: {}", diffusivity(nt, nu));
            assert_eq!(nu_t(nt, nu), 0.0);
        }
    }

    /// ν_t and the sources are finite and ν_t ≥ 0 over a 10⁴-point hostile grid, d = 0 included.
    #[test]
    fn hostile_inputs_stay_finite() {
        let nu = 1e-2;
        let mut x = 0.371_f64;
        let mut next = || {
            x = (x * 9301.0 + 0.4927).fract();
            x
        };
        let specials = [0.0, f64::NAN, f64::INFINITY, -f64::INFINITY, f64::MIN_POSITIVE, -0.0];
        for k in 0..10_000 {
            let nt = (2.0 * next() - 1.0) * 1e3 * nu;
            let d = if k % 97 == 0 { 0.0 } else { next() * 1e6 };
            let omega = next() * 1e6;
            let v = nu_t(nt, nu);
            assert!(v.is_finite() && v >= 0.0, "{nt}");
            let s = sources(nt, nu, omega, d);
            assert!(s.production.is_finite() && s.destruction.is_finite(), "{nt} {d} {omega}");
            assert!(diffusivity(nt, nu) > 0.0);
        }
        for a in specials {
            for b in specials {
                assert!(nu_t(a, nu) >= 0.0 && nu_t(a, b) >= 0.0);
                let s = sources(a, nu, b, a);
                assert!(s.production.is_finite() && s.destruction.is_finite());
                let s = sources(1e-3, nu, a, b);
                assert!(s.production.is_finite() && s.destruction.is_finite());
            }
        }
    }

    #[test]
    fn vorticity_of_shear_and_rotation() {
        let shear = [[0.0, 0.3, 0.0], [0.0; 3], [0.0; 3]];
        assert!((vorticity(&shear) - 0.3).abs() < 1e-15);
        let rot = [[0.0, 0.5, 0.0], [-0.5, 0.0, 0.0], [0.0; 3]];
        assert!((vorticity(&rot) - 1.0).abs() < 1e-15);
        let strain = [[0.2, 0.0, 0.0], [0.0, -0.2, 0.0], [0.0; 3]];
        assert_eq!(vorticity(&strain), 0.0);
        assert_eq!(vorticity(&[[f64::NAN; 3]; 3]), 0.0);
    }
}
