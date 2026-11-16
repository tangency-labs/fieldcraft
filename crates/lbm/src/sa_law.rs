//! The Spalart–Allmaras constant-stress wall law (`wall_model: "spalart_allmaras_bounce"`): the
//! mean velocity u⁺(y⁺) that the Spalart–Allmaras model itself produces in a layer of constant
//! total stress, as a wall law for the on-node bounce ([`crate::bounce`]) beside Spalding's.
//!
//! # Sources
//! - P. R. Spalart, S. R. Allmaras, *A one-equation turbulence model for aerodynamic flows*,
//!   La Recherche Aérospatiale 1994(1):5–21 (first as AIAA-92-0439, 1992). The model's near-wall
//!   behaviour: its working variable is ν̃ = κ u_τ y in the log layer, where its production and
//!   destruction balance (their §3).
//! - R. H. Nichols, C. C. Nelson, *Wall function boundary conditions including heat transfer and
//!   compressibility*, AIAA J. 42(6):1107–1114, 2004. The idea of a wall function that matches
//!   the first node's ν̃ to the model's own wall behaviour rather than to another model's.
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
//! rest. So the wall law that goes with the model is the integral
//!
//! ```text
//! u⁺(y⁺) = ∫₀^y⁺ dy / (1 + κ y f_v1(κ y)),        ν̃ = κ u_τ y
//! ```
//!
//! with nothing fitted: the same κ as the model (0.41). Its log-law intercept is about 5.0 (u⁺ =
//! 21.89 at y⁺ 1000). It differs from Spalding's law (κ 0.41, B 5.5) in the buffer layer: at y⁺ 15
//! it gives 11.09 where Spalding's gives 10.4 and Lee & Moser's DNS 10.9, and Spalding's law is 0.8
//! low at y⁺ 30. A wall function that imposes ν̃ at the first node and carries a stress at it
//! should use the pair the model is in equilibrium with: with Spalding's law the node's ν_t is 60 %
//! above the model's own at y⁺ 7.5 and the profile above it sits 0.4 to 0.7 low (docs/verification.md).
//!
//! # The table
//! The integral is tabulated once, at first use, on a grid that grows geometrically by 0.2 % from
//! y⁺ = 10⁻³ to 10⁷ (about 11 500 nodes, at most [`MAX_NODES`]), by the trapezoid rule on the
//! integrand, which needs only the four arithmetic operations: the table is the same bit for bit on
//! every platform. Below the first node u⁺ = y⁺ (the integrand differs from 1 by O(y⁴)); above the
//! last, y⁺ saturates (a wall at y⁺ ≥ 10⁷ is outside any lattice Boltzmann grid). Lookup is a binary
//! search and linear interpolation in y⁺.
//!
//! [`friction_velocity`] inverts the law: the friction velocity of a flow of speed U at distance y
//! from the wall with viscosity ν solves u_τ u⁺(y u_τ / ν) = U, i.e. y⁺ u⁺(y⁺) = U y / ν, whose
//! left side rises with y⁺; the root is found in the table cell that brackets it by solving the
//! interpolant's quadratic exactly, so the inverse of the interpolant is exact to round-off.

use crate::sa::{KAPPA, f_v1};
use std::sync::OnceLock;

/// The first tabulated height in wall units; below it u⁺ = y⁺.
pub const Y_PLUS_MIN: f64 = 1e-3;
/// The last tabulated height in wall units; above it the law saturates.
pub const Y_PLUS_MAX: f64 = 1e7;
/// The growth of the grid in y⁺ from one node to the next.
const RATIO: f64 = 1.002;
/// The most nodes the table takes (the grid reaches [`Y_PLUS_MAX`] in about 11 500).
pub const MAX_NODES: usize = 20_000;

/// The eddy viscosity ν_t⁺ = κ y⁺ f_v1(κ y⁺) of the equilibrium layer, for ν̃ = κ y⁺. 0 at
/// y⁺ ≤ 0 or for a number that is not usable.
pub fn eddy_viscosity_plus(y_plus: f64) -> f64 {
    if !(y_plus > 0.0 && y_plus.is_finite()) {
        return 0.0;
    }
    let chi = KAPPA * y_plus;
    let nt = chi * f_v1(chi);
    if nt.is_finite() { nt } else { 0.0 }
}

/// The shear du⁺/dy⁺ = 1 / (1 + ν_t⁺).
fn shear(y_plus: f64) -> f64 {
    1.0 / (1.0 + eddy_viscosity_plus(y_plus))
}

/// The tabulated law: nodes y⁺ᵢ, u⁺ᵢ and the products y⁺ᵢ u⁺ᵢ the inversion searches.
struct Table {
    y: Vec<f64>,
    u: Vec<f64>,
    product: Vec<f64>,
}

fn table() -> &'static Table {
    static TABLE: OnceLock<Table> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut y = vec![0.0, Y_PLUS_MIN];
        let mut u = vec![0.0, Y_PLUS_MIN];
        while y.len() < MAX_NODES {
            let (Some(&y0), Some(&u0)) = (y.last(), u.last()) else { break };
            if y0 >= Y_PLUS_MAX {
                break;
            }
            let y1 = y0 * RATIO;
            u.push(u0 + 0.5 * (shear(y0) + shear(y1)) * (y1 - y0));
            y.push(y1);
        }
        let product = y.iter().zip(&u).map(|(y, u)| y * u).collect();
        Table { y, u, product }
    })
}

/// u⁺(y⁺): the mean velocity in wall units of the equilibrium layer at height `y_plus`. 0 at
/// y⁺ ≤ 0 or for a number that is not usable; the value at [`Y_PLUS_MAX`] above it.
pub fn u_plus(y_plus: f64) -> f64 {
    if !(y_plus > 0.0) {
        return 0.0;
    }
    let t = table();
    let (Some(&y_last), Some(&u_last)) = (t.y.last(), t.u.last()) else { return 0.0 };
    if y_plus <= Y_PLUS_MIN {
        return y_plus;
    }
    if !(y_plus < y_last) {
        return u_last;
    }
    let k = t.y.partition_point(|y| *y < y_plus).clamp(1, t.y.len() - 1);
    let (Some(&ya), Some(&yb), Some(&ua), Some(&ub)) = (t.y.get(k - 1), t.y.get(k), t.u.get(k - 1), t.u.get(k)) else { return 0.0 };
    ua + (ub - ua) * (y_plus - ya) / (yb - ya)
}

/// The y⁺ with y⁺ u⁺(y⁺) = `re`, the inverse of the law's product (a Reynolds number U y / ν).
/// 0 for `re` ≤ 0 or not a number; [`Y_PLUS_MAX`] for a `re` the table does not reach.
pub fn y_plus_for_reynolds(re: f64) -> f64 {
    if !(re > 0.0) {
        return 0.0;
    }
    let t = table();
    let (Some(&y_last), Some(&p_last)) = (t.y.last(), t.product.last()) else { return 0.0 };
    // Below the first node u⁺ = y⁺, so the product is y⁺².
    if re <= Y_PLUS_MIN * Y_PLUS_MIN {
        return re.sqrt();
    }
    if !(re < p_last) {
        return y_last;
    }
    let k = t.product.partition_point(|p| *p <= re).clamp(1, t.y.len() - 1);
    let (Some(&ya), Some(&yb), Some(&ua), Some(&ub)) = (t.y.get(k - 1), t.y.get(k), t.u.get(k - 1), t.u.get(k)) else { return 0.0 };
    // In the cell u⁺ = ua + s (y − ya), so y u⁺ = re is s y² + b y − re = 0 with b = ua − s ya.
    let s = (ub - ua) / (yb - ya);
    let b = ua - s * ya;
    let disc = (b * b + 4.0 * s * re).max(0.0).sqrt();
    // The positive root, in the form that does not cancel.
    let y = if b >= 0.0 { 2.0 * re / (b + disc) } else { (disc - b) / (2.0 * s) };
    if y.is_finite() { y.clamp(ya, yb) } else { ya }
}

/// The growth of the table's grid from one node to the next (the GPU solver reads it to bound a
/// segment from above).
pub const GRID_RATIO: f64 = RATIO;

/// The table as the segments [`y_plus_for_reynolds`] searches, for a device that has no way to
/// share the host's table: segment k (k = 1, …) between nodes k − 1 and k is `[P_k, y_(k−1), s, b]`
/// with P_k = y⁺_k u⁺_k, the largest Reynolds number it holds, the node it starts at, and the slope
/// s and intercept b of u⁺ = b + s y⁺ inside it. The same table, bit for bit, as the host's.
pub fn segments() -> Vec<[f64; 4]> {
    let t = table();
    let mut out = Vec::with_capacity(t.y.len().saturating_sub(1));
    for k in 1..t.y.len() {
        let (Some(&ya), Some(&yb), Some(&ua), Some(&ub), Some(&p)) = (t.y.get(k - 1), t.y.get(k), t.u.get(k - 1), t.u.get(k), t.product.get(k))
        else {
            break;
        };
        let s = (ub - ua) / (yb - ya);
        out.push([p, ya, s, ua - s * ya]);
    }
    out
}

/// The friction velocity of a flow of speed `speed` parallel to a wall at distance `y` with
/// kinematic viscosity `nu`, from the equilibrium layer: u_τ with u_τ u⁺(y u_τ / ν) = U. 0 for a
/// speed, distance or viscosity that is not a positive number.
pub fn friction_velocity(speed: f64, y: f64, nu: f64) -> f64 {
    if !(speed > 0.0 && speed.is_finite() && y > 0.0 && y.is_finite() && nu > 0.0 && nu.is_finite()) {
        return 0.0;
    }
    let y_plus = y_plus_for_reynolds(speed * y / nu);
    let ut = nu * y_plus / y;
    if ut.is_finite() { ut } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa::{C_B1, C_B2, SIGMA, c_w1, f_t2, f_v2, f_w};

    /// The law's limits: u⁺ = y⁺ in the sublayer (to O(y⁴)), the log law of the model's κ above
    /// (u⁺ − ln y⁺ / κ settles to a constant, about 5.0), monotone, and u⁺ ≤ y⁺ everywhere.
    #[test]
    fn the_limits_are_the_sublayer_and_the_log_law() {
        for y in [1e-6, 1e-3, 0.01, 0.1, 0.5, 1.0] {
            // u⁺ = y⁺ − (κ⁴ / c_v1³) y⁺⁵ / 5 + …
            assert!((u_plus(y) - y).abs() <= 2e-5 * y.powi(5) + 1e-12 * y, "{y}: {}", u_plus(y));
        }
        assert_eq!(u_plus(0.0), 0.0);
        let mut last = 0.0;
        for k in 0..=2000 {
            let y = 10f64.powf(-3.0 + 10.0 * k as f64 / 2000.0);
            let u = u_plus(y);
            assert!(u > last && u <= y + 1e-12, "y⁺ {y}: {u}");
            last = u;
        }
        // The intercept B = u⁺ − ln y⁺ / κ: constant to a hundredth from y⁺ 1000 on.
        let b = |y: f64| u_plus(y) - y.ln() / KAPPA;
        for y in [1e3, 1e4, 1e5, 1e6] {
            assert!((b(y) - b(1e3)).abs() < 1e-2, "y⁺ {y}: B {}", b(y));
        }
        assert!((b(1e3) - 5.04).abs() < 0.05, "B {}", b(1e3));
        // The values the wall function is checked at (docs/verification.md).
        for (y, u) in [(7.5, 7.183), (15.0, 11.089), (30.0, 13.383), (60.0, 15.100), (1000.0, 21.890)] {
            assert!((u_plus(y) - u).abs() < 5e-3, "y⁺ {y}: {} vs {u}", u_plus(y));
        }
    }

    /// The segments a device searches give the inversion's answer: the first one whose largest
    /// Reynolds number is above `re`, its quadratic, clamped to the segment.
    #[test]
    fn the_segments_invert_like_the_table() {
        let seg = segments();
        assert!(seg.len() < MAX_NODES && seg.windows(2).all(|w| w[0][0] < w[1][0] && w[0][1] < w[1][1]));
        for e in 0..=60 {
            let re = 2e-6 * 10f64.powf(e as f64 / 5.0);
            let k = seg.partition_point(|s| s[0] <= re);
            let Some(&[_, ya, s, b]) = seg.get(k) else {
                assert!(re >= seg[seg.len() - 1][0] || y_plus_for_reynolds(re) > 0.0);
                continue;
            };
            let disc = (b * b + 4.0 * s * re).max(0.0).sqrt();
            let y = if b >= 0.0 { 2.0 * re / (b + disc) } else { (disc - b) / (2.0 * s) };
            let y = y.clamp(ya, ya * GRID_RATIO);
            let want = y_plus_for_reynolds(re);
            assert!((y - want).abs() <= 1e-12 * want.max(1e-3), "re {re}: {y} vs {want}");
        }
    }

    /// The table against an independent quadrature of the same integral: 64-point composite
    /// Gauss–Legendre panels in ln y (a different rule on a different variable).
    #[test]
    fn the_table_matches_an_independent_quadrature() {
        // 5-point Gauss–Legendre on [-1, 1].
        let nodes = [0.0, 0.538_469_310_105_683_1, -0.538_469_310_105_683_1, 0.906_179_845_938_664, -0.906_179_845_938_664];
        let weights = [0.568_888_888_888_888_9, 0.478_628_670_499_366_5, 0.478_628_670_499_366_5, 0.236_926_885_056_189_1, 0.236_926_885_056_189_1];
        let integral = |top: f64| -> f64 {
            // Panels of ln y from ln 1e-3 to ln top, 32 per decade; below 1e-3 the integrand is 1.
            let (lo, hi) = (Y_PLUS_MIN.ln(), top.ln());
            let panels = ((hi - lo) / std::f64::consts::LN_10 * 32.0).ceil().max(1.0) as usize;
            let h = (hi - lo) / panels as f64;
            let mut sum = Y_PLUS_MIN;
            for p in 0..panels {
                let mid = lo + (p as f64 + 0.5) * h;
                for (x, w) in nodes.iter().zip(weights) {
                    let y = (mid + 0.5 * h * x).exp();
                    sum += 0.5 * h * w * y * shear(y);
                }
            }
            sum
        };
        for y in [0.01, 0.5, 3.0, 7.5, 15.0, 30.0, 60.0, 180.0, 1000.0, 5200.0, 1e5, 1e7 - 1.0] {
            let (a, b) = (u_plus(y), integral(y));
            assert!((a - b).abs() <= 2e-6 * b.max(1.0), "y⁺ {y}: table {a}, quadrature {b}");
        }
    }

    /// χ = κ y⁺ solves the model's one-dimensional steady equation at every height, buffer layer
    /// included, with the crate's own f_v2, f_w, f_t2 and constants: the residual of P − D plus
    /// diffusion is round-off against the production.
    #[test]
    fn the_equilibrium_layer_solves_the_models_equation() {
        for k in 0..=120 {
            let y = 0.05 * 10f64.powf(5.0 * k as f64 / 120.0);
            let chi = KAPPA * y;
            let omega = shear(y);
            let s_bar = chi * f_v2(chi) / (KAPPA * y).powi(2);
            let s_tilde = omega + s_bar;
            assert!((s_tilde * KAPPA * y - 1.0).abs() < 1e-12, "y⁺ {y}: S̃ κ y = {}", s_tilde * KAPPA * y);
            let r = (chi / (s_tilde * KAPPA * KAPPA * y * y)).min(10.0);
            assert!((r - 1.0).abs() < 1e-12, "y⁺ {y}: r {r}");
            let ft2 = f_t2(chi);
            let production = C_B1 * (1.0 - ft2) * s_tilde * chi;
            let destruction = (c_w1() * f_w(r) - C_B1 * ft2 / (KAPPA * KAPPA)) * (chi / y) * (chi / y);
            let diffusion = (1.0 + C_B2) * KAPPA * KAPPA / SIGMA;
            let residual = production - destruction + diffusion;
            assert!(residual.abs() <= 1e-11 * production.abs().max(1e-3), "y⁺ {y}: residual {residual}");
        }
    }

    /// The model's one-dimensional solver relaxed from a different start reaches the law: a column
    /// of SA cells (unit lattice spacing, an on-node wall one cell below the first) with u_τ = 1 and
    /// ν = 1/h, so a cell is h wall units; the shear is recomputed each step from the current ν̃
    /// (constant total stress) and the top face is held at the law's ν̃. The steady ν̃ is κ u_τ y and the
    /// velocity profile integrates to the table.
    #[test]
    fn the_models_one_dimensional_solver_relaxes_to_the_law() {
        use crate::sa::{SaField, SaInput, Side, nu_t};
        let (h, n) = (0.75, 28usize);
        let nu = 1.0 / h;
        let open = Side::Open;
        let sides: Vec<[Side; 6]> = (0..n)
            .map(|j| {
                let below = if j == 0 { Side::Node } else { Side::Cell(j as u32 - 1) };
                let above = if j + 1 == n { Side::Inflow } else { Side::Cell(j as u32 + 1) };
                [open, open, below, above, open, open]
            })
            .collect();
        let mut field = SaField::new(sides, vec![true; n], KAPPA * (n as f64 + 0.5), 0.0).unwrap();
        // Lattice distance from the wall; wall units are h times it.
        let y: Vec<f64> = (0..n).map(|j| (j + 1) as f64).collect();
        // Start from 1.1 times the law's ν̃ in the lower half and 0.9 of it above.
        let start: Vec<f64> = y.iter().map(|y| KAPPA * y * if y * h < 9.0 { 1.1 } else { 0.9 }).collect();
        field.set_nu_tilde(&start).unwrap();
        let u = vec![[0.0; 3]; n];
        // As in the wall function, the first fluid node is imposed.
        let bottom = KAPPA * y[0];
        for _ in 0..4_000 {
            let vorticity: Vec<f64> = field.nu_tilde().iter().map(|v| 1.0 / (nu + nu_t(*v, nu))).collect();
            let input = SaInput { u: &u, vorticity: &vorticity, length: &y, nu, model: true, source: None, substeps: Some(64) };
            field.advance(&input).unwrap();
            field.impose(&[(0, bottom)]).unwrap();
        }
        // ν̃ against κ u_τ y above the imposed first node (the first cells hold the sublayer).
        for (j, yl) in y.iter().enumerate().filter(|(_, yl)| (1.0..=18.0).contains(&(*yl * h))) {
            let v = field.nu_tilde().get(j).copied().unwrap_or(f64::NAN);
            assert!((v / (KAPPA * yl) - 1.0).abs() < 0.03, "y⁺ {}: ν̃ {v} vs κ y {}", yl * h, KAPPA * yl);
        }
        // The velocity the settled ν̃ implies, integrated from the wall, against the table.
        let mut u_col = u_plus(h);
        let mut prev = (h, 1.0 / (1.0 + nu_t(bottom, nu) / nu));
        for (j, yl) in y.iter().enumerate().skip(1) {
            let yp = yl * h;
            if yp > 18.0 {
                break;
            }
            let shear = 1.0 / (1.0 + nu_t(field.nu_tilde().get(j).copied().unwrap_or(0.0), nu) / nu);
            u_col += 0.5 * (prev.1 + shear) * (yp - prev.0);
            prev = (yp, shear);
            if [7.5, 15.0].iter().any(|t| (yp - t).abs() < 1e-9) {
                assert!((u_col - u_plus(yp)).abs() < 0.05, "y⁺ {yp}: u⁺ {u_col} vs {}", u_plus(yp));
            }
        }
    }

    /// The inversion: y⁺ u⁺(y⁺) = Re recovers y⁺ to round-off over the whole table, u_τ solves
    /// u_τ u⁺(y u_τ / ν) = U, and the stress ρ u_τ² rises with the speed.
    #[test]
    fn the_inversion_round_trips() {
        for k in 0..=4000 {
            let y = 10f64.powf(-3.5 + 10.4 * k as f64 / 4000.0);
            let re = y * u_plus(y);
            let back = y_plus_for_reynolds(re);
            assert!((back - y).abs() <= 1e-9 * y, "y⁺ {y}: back {back}");
        }
        let (y, nu) = (1.0, 1.57e-3);
        let mut last = 0.0;
        for k in 1..=200 {
            let speed = 1e-4 * 1.07f64.powi(k);
            let ut = friction_velocity(speed, y, nu);
            assert!(ut > last, "speed {speed}");
            last = ut;
            let check = ut * u_plus(y * ut / nu);
            assert!((check - speed).abs() <= 1e-9 * speed, "speed {speed}: {check}");
        }
    }

    /// Hostile input: zero, negative, infinite and not-a-number arguments give 0 or saturate, never
    /// a panic or a non-finite number; the table is the same on every call.
    #[test]
    fn hostile_input_gives_zero_or_saturates() {
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
        // Above the table the law saturates.
        let last = table().y.last().copied().unwrap_or(Y_PLUS_MAX);
        assert_eq!(u_plus(1e9), u_plus(last));
        assert!(table().y.len() <= MAX_NODES && table().y.last().is_some_and(|y| *y >= Y_PLUS_MAX));
    }
}
