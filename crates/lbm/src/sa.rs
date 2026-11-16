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
/// κ, von Kármán's constant: the one of Han, Ooka & Kikumoto's wall-function bounce, so the
/// model's log layer and the wall's agree.
pub const KAPPA: f64 = crate::wall::BOUNCE_KAPPA;
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
    C_T3 * (-C_T4 * chi * chi).exp()
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
    g * ((1.0 + c6) / (g.powi(6) + c6)).powf(1.0 / 6.0)
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

/// The wall function at the first fluid node (Nichols and Nelson, *Wall function boundary
/// conditions including heat transfer and compressibility*, AIAA J. 42(6):1107–1114, 2004): its
/// eddy viscosity is the one of the law the wall itself applies, at the node's y⁺, and its ν̃
/// the one whose ν̃ f_v1(ν̃/ν) gives that ν_t ([`nu_tilde_for`]). The law's eddy viscosity in
/// wall units is ν_t⁺ = dy⁺/du⁺ − 1 (the total stress is the wall's, (1 + ν_t⁺) du⁺/dy⁺ = 1).
///
/// In the buffer and lower log layer this ν̃ sits 10 to 16 % below SA's own log-layer value
/// κ u_τ y (y⁺ 30 to 100 on Spalding's law): the law's ν_t⁺ lags κ y⁺ by about κ u⁺
/// (Spalding's ν_t⁺ = κ y⁺ − κ u⁺ + κ e^{−κB} (κu⁺)³/6). The node then carries the stress the wall
/// imposes, and SA's transport takes over from the second node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WallLaw {
    /// Spalding's law (κ 0.41, B 5.5), as Han, Ooka and Kikumoto's on-node bounce applies it.
    Spalding,
    /// Musker's profile, as the Musker wall model applies it.
    Musker,
}

/// The law's eddy viscosity ν_t⁺ = dy⁺/du⁺ − 1 at `y_plus`, analytically: Spalding's from its
/// u⁺ (found by bisection, at most [`INVERT_ITERATIONS`] halvings), Musker's the mixing-length ν_t
/// the Musker wall model's own law applies ([`crate::wall::wall_law`]), so that ν̃ and the wall
/// cell's collision agree. 0 at y⁺ ≤ 0 or for a y⁺ that is not a number; Spalding's may be ≤ 0 in
/// the viscous sublayer. An `Err` only if the bisection hits its cap, which takes a y⁺ above about
/// 1e45.
pub fn law_eddy_viscosity(law: WallLaw, y_plus: f64) -> crate::Result<f64> {
    if !(y_plus > 0.0 && y_plus.is_finite()) {
        return Ok(0.0);
    }
    let out = match law {
        WallLaw::Spalding => {
            let (k, e) = (crate::wall::BOUNCE_KAPPA, (-crate::wall::BOUNCE_KAPPA * crate::wall::BOUNCE_B).exp());
            // y⁺(u⁺) = u⁺ + e^{−κB}(e^{κu⁺} − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6) rises from 0 and is
            // at least u⁺, so u⁺ lies in [0, y⁺].
            let y_of = |u: f64| {
                let x = k * u;
                u + e * (x.exp() - 1.0 - x - x * x / 2.0 - x * x * x / 6.0)
            };
            let u = bisect(0.0, y_plus, |u| y_of(u) < y_plus, "Spalding's u⁺")?;
            let x = k * u;
            k * e * (x.exp() - 1.0 - x - x * x / 2.0)
        }
        WallLaw::Musker => crate::wall::musker_eddy_viscosity(y_plus),
    };
    Ok(if out.is_finite() { out } else { 0.0 })
}

/// The ν̃ ≥ 0 whose eddy viscosity ν̃ f_v1(ν̃/ν) is `nu_t`: [`WallLaw`]'s first-node ν̃. ν̃ f_v1
/// rises strictly from 0, so the root is unique, and with χ = ν̃/ν and t = ν_t/ν it lies in
/// [t, t + c_v1] (f_v1 ≤ 1 below; above, χ f_v1 ≥ χ − c_v1³/χ² ≥ t at χ = t + c_v1). Found by
/// bisection on χ. 0 where `nu_t` ≤ 0 (the viscous sublayer) or `nu_t`, `nu` is not a usable
/// number; an `Err` only if the bisection hits its cap.
pub fn nu_tilde_for(nu_t: f64, nu: f64) -> crate::Result<f64> {
    let t = nu_t / nu;
    if !(nu > 0.0 && nu.is_finite() && t > 0.0 && t.is_finite()) {
        return Ok(0.0);
    }
    Ok(bisect(t, t + C_V1, |chi| chi * f_v1(chi) < t, "the wall function's ν̃")? * nu)
}

/// The most halvings a wall-function bisection takes; one stops once its bracket is below
/// [`INVERT_TOLERANCE`] relative (about 50 halvings for these brackets).
pub const INVERT_ITERATIONS: u32 = 200;
/// The relative width at which a wall-function bisection stops.
pub const INVERT_TOLERANCE: f64 = 1e-15;

/// The root in [lo, hi] of a rising function, `below(x)` true left of it: the bracket's upper end
/// once it is narrower than [`INVERT_TOLERANCE`] relative or at f64 resolution; an `Err` naming
/// `what` after [`INVERT_ITERATIONS`] halvings.
fn bisect(mut lo: f64, mut hi: f64, below: impl Fn(f64) -> bool, what: &str) -> crate::Result<f64> {
    for _ in 0..INVERT_ITERATIONS {
        let mid = 0.5 * (lo + hi);
        if hi - lo <= INVERT_TOLERANCE * hi.abs() || !(mid > lo && mid < hi) {
            return Ok(hi);
        }
        if below(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(crate::Error::Unstable(format!("Spalart–Allmaras: {what} did not converge in {INVERT_ITERATIONS} halvings")))
}

/// The source terms of a cell, split for a point-implicit update: the rate of ν̃ is
/// `production − destruction`. Both 0 where an input is not usable or `d` is not above 0.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Sources {
    pub production: f64,
    pub destruction: f64,
}

/// P and D of the cell (see the module docs) for ν̃ `nu_tilde`, molecular viscosity `nu`,
/// vorticity magnitude `omega` and length `d` (the wall distance).
pub fn sources(nu_tilde: f64, nu: f64, omega: f64, d: f64) -> Sources {
    coefficients(nu_tilde, nu, omega, d)
        .map(|(p, k)| Sources { production: p * nu_tilde, destruction: k * nu_tilde })
        .filter(|s| s.production.is_finite() && s.destruction.is_finite())
        .unwrap_or_default()
}

/// P and D per unit ν̃, `(p, k)` with P = p ν̃ and D = k ν̃, formed without dividing by ν̃ (the
/// point-implicit update needs D/ν̃ and, for SA-neg, (P − D)/ν̃; a quotient would lose them near
/// ν̃ = 0, and the GPU's f32 more so). `None` where an input is not usable or `d` is not above 0.
fn coefficients(nu_tilde: f64, nu: f64, omega: f64, d: f64) -> Option<(f64, f64)> {
    if !(nu > 0.0 && d > 0.0 && omega >= 0.0 && nu_tilde.is_finite() && omega.is_finite() && !d.is_nan()) {
        return None;
    }
    let (p, k) = if nu_tilde >= 0.0 {
        let ft2 = f_t2(nu_tilde / nu);
        let st = s_tilde(omega, nu_tilde, nu, d);
        let denominator = st * KAPPA * KAPPA * d * d;
        let r = if denominator > 0.0 { (nu_tilde / denominator).min(R_MAX) } else { R_MAX };
        (C_B1 * (1.0 - ft2) * st, (c_w1() * f_w(r) - C_B1 / (KAPPA * KAPPA) * ft2) * nu_tilde / (d * d))
    } else {
        (C_B1 * (1.0 - C_T3) * omega, -c_w1() * nu_tilde / (d * d))
    };
    (p.is_finite() && k.is_finite()).then_some((p, k))
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

/// What lies across one face of a cell, for ν̃'s transport.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// A fluid cell (its index).
    Cell(u32),
    /// A no-slip wall half a cell from the cell centre (halfway bounce-back): ν̃ = 0 there, no
    /// flow through it.
    Wall,
    /// A no-slip wall on the boundary node one full cell away (Han, Ooka & Kikumoto's on-node
    /// wall): ν̃ = 0 on the node, no flow through it.
    Node,
    /// An inflow face half a cell away: ν̃ is the field's inflow value there.
    Inflow,
    /// An outflow, slip or far-field face: zero gradient of ν̃, no diffusion through it.
    Open,
}

/// What a face of the box is, for [`SaField::boxed`].
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoxFace {
    Periodic,
    Wall,
    Inflow,
    Open,
}

/// The six faces of a cell in the order −x, +x, −y, +y, −z, +z.
const AXIS_OF_SIDE: [usize; 6] = [0, 0, 1, 1, 2, 2];
const SIGN_OF_SIDE: [f64; 6] = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0];

/// The most transport sub-steps one lattice step may take.
pub const MAX_SUBSTEPS: u32 = 64;

/// Cells per parallel block (fixed, so the result does not depend on the thread count).
const BLOCK: usize = 4096;

/// The inputs of one lattice step of ν̃ (every slice one entry per cell, lattice units).
#[derive(Debug, Clone, Copy)]
pub struct SaInput<'a> {
    /// The cell velocities of the pass before the collision.
    pub u: &'a [[f64; 3]],
    /// The vorticity magnitude Ω of each cell.
    pub vorticity: &'a [f64],
    /// The length of the destruction term, lattice units: the wall distance d (RANS, from
    /// [`WallDistance::as_slice`](crate::distance::WallDistance::as_slice)).
    pub length: &'a [f64],
    /// The molecular kinematic viscosity.
    pub nu: f64,
    /// Whether the model's production and destruction act (off: pure transport, for verification).
    pub model: bool,
    /// An extra explicit source per cell (rate per lattice step), for manufactured solutions.
    pub source: Option<&'a [f64]>,
    /// `Some(n)`: exactly n equal sub-steps (1 to [`MAX_SUBSTEPS`]), the reference a device
    /// backend with a fixed count is checked against, refused as unstable if n breaks the
    /// monotone bound; `None`: as many as the bound needs.
    pub substeps: Option<u32>,
}

/// What one lattice step of ν̃ did.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SaStats {
    pub substeps: u32,
    /// The largest ν_t of a fluid cell (lattice units).
    pub max_nu_t: f64,
    /// The smallest ν̃ of a fluid cell.
    pub min_nu_tilde: f64,
    /// Fluid cells with ν̃ below 0 (SA-neg's branch).
    pub negative_cells: usize,
}

/// The Spalart–Allmaras working variable ν̃ on the lattice's cells, and its transport.
///
/// **Scheme.** Per sub-step, explicit on the cell grid (Δx = 1): first-order upwind advection in
/// flux form; conservative central diffusion `(1/σ) Σ_f η_f (ν̃_nb − ν̃_c)` with `η = ν + ν̃`
/// (`ν + ν̃ f_n` below 0), the face value half a cell away at a wall or inflow face; the
/// `c_b2 |∇ν̃|² / σ` term as a non-negative source; production and destruction point-implicitly
/// ([`sources`]), so a positive ν̃ stays positive and a negative one decays towards 0.
///
/// Every neighbour's coefficient is non-negative, so with `Δt Σ coefficients ≤ 1` (the sub-step) a
/// field that starts non-negative stays non-negative, and ν_t = max(0, ν̃ f_v1) always is. The
/// sub-steps are at most [`MAX_SUBSTEPS`].
#[derive(Debug, Clone, PartialEq)]
pub struct SaField {
    nu_tilde: Vec<f64>,
    sides: Vec<[Side; 6]>,
    fluid: Vec<bool>,
    inflow: f64,
    /// Scratch of [`SaField::advance`], kept for its capacity: each cell's transport rate and
    /// coefficient sum, and the next sub-step's ν̃.
    rates: Vec<(f64, f64)>,
    next: Vec<f64>,
}

/// The bytes a field takes per cell: ν̃, its six faces, its fluid flag and the scratch.
pub const BYTES_PER_CELL: u64 = (size_of::<f64>() + size_of::<[Side; 6]>() + size_of::<bool>() + size_of::<(f64, f64)>() + size_of::<f64>()) as u64;

impl SaField {
    /// A field over cells with the given faces, every fluid cell at `initial`, `inflow` at inflow
    /// faces. A side's cell must be a fluid cell of the field.
    pub fn new(sides: Vec<[Side; 6]>, fluid: Vec<bool>, inflow: f64, initial: f64) -> crate::Result<Self> {
        if sides.len() != fluid.len() {
            return Err(crate::Error::Invalid(format!("Spalart–Allmaras: {} cells of faces for {} cells", sides.len(), fluid.len())));
        }
        if !(inflow.is_finite() && initial.is_finite()) {
            return Err(crate::Error::Invalid(format!("Spalart–Allmaras: the inflow ν̃ {inflow} and initial ν̃ {initial} must be finite")));
        }
        for (c, s) in sides.iter().enumerate() {
            for side in s {
                if let Side::Cell(n) = *side
                    && !fluid.get(n as usize).copied().unwrap_or(false)
                {
                    return Err(crate::Error::Invalid(format!("Spalart–Allmaras: cell {c} faces {n}, which is not a fluid cell")));
                }
            }
        }
        let n = fluid.len();
        let mut nu_tilde = crate::flow::filled(n, 0.0, "Spalart–Allmaras ν̃")?;
        for (v, &f) in nu_tilde.iter_mut().zip(&fluid) {
            if f {
                *v = initial;
            }
        }
        let rates = crate::flow::filled(n, (0.0, 0.0), "Spalart–Allmaras rates")?;
        let next = crate::flow::filled(n, 0.0, "Spalart–Allmaras ν̃")?;
        Ok(SaField { nu_tilde, sides, fluid, inflow, rates, next })
    }

    /// A field on a box of `dims` cells (x fastest), `fluid` per cell, each box face of `faces`
    /// (−x, +x, −y, +y, −z, +z; periodic in pairs). A face to a non-fluid cell is a wall.
    #[doc(hidden)]
    pub fn boxed(dims: [usize; 3], fluid: &[bool], faces: [BoxFace; 6], inflow: f64, initial: f64) -> crate::Result<Self> {
        let n = dims[0].checked_mul(dims[1]).and_then(|v| v.checked_mul(dims[2]));
        if n != Some(fluid.len()) || fluid.len() > u32::MAX as usize - 1 {
            return Err(crate::Error::Invalid(format!("Spalart–Allmaras: a {dims:?} box for {} cells", fluid.len())));
        }
        for a in 0..3 {
            if (faces[2 * a] == BoxFace::Periodic) != (faces[2 * a + 1] == BoxFace::Periodic) {
                return Err(crate::Error::Invalid(format!("Spalart–Allmaras: box axis {a} periodic on one face only")));
            }
        }
        let stride = [1, dims[0], dims[0] * dims[1]];
        let mut sides = vec![[Side::Open; 6]; fluid.len()];
        for (c, s) in sides.iter_mut().enumerate() {
            let at = [c % dims[0], (c / dims[0]) % dims[1], c / (dims[0] * dims[1])];
            for k in 0..6 {
                let a = AXIS_OF_SIDE[k];
                let up = k % 2 == 1;
                let edge = if up { at[a] + 1 == dims[a] } else { at[a] == 0 };
                let nb = if !edge {
                    Some(if up { c + stride[a] } else { c - stride[a] })
                } else if faces[k] == BoxFace::Periodic {
                    Some(if up { c + stride[a] - dims[a] * stride[a] } else { c + (dims[a] - 1) * stride[a] })
                } else {
                    None
                };
                s[k] = match nb {
                    Some(nb) if fluid.get(nb).copied().unwrap_or(false) => Side::Cell(nb as u32),
                    Some(_) => Side::Wall,
                    None => match faces[k] {
                        BoxFace::Wall | BoxFace::Periodic => Side::Wall,
                        BoxFace::Inflow => Side::Inflow,
                        BoxFace::Open => Side::Open,
                    },
                };
            }
        }
        SaField::new(sides, fluid.to_vec(), inflow, initial)
    }

    /// ν̃ per cell (lattice units; 0 in non-fluid cells).
    pub fn nu_tilde(&self) -> &[f64] {
        &self.nu_tilde
    }

    /// Each cell's six faces (−x, +x, −y, +y, −z, +z): what a device backend checks its own
    /// neighbour rule against.
    pub fn sides(&self) -> &[[Side; 6]] {
        &self.sides
    }

    /// ν̃ at inflow faces.
    pub fn inflow(&self) -> f64 {
        self.inflow
    }

    /// Set ν̃ per cell (a checkpoint, a test); refused unless one finite value per cell.
    pub fn set_nu_tilde(&mut self, values: &[f64]) -> crate::Result<()> {
        if values.len() != self.nu_tilde.len() || values.iter().any(|v| !v.is_finite()) {
            return Err(crate::Error::Invalid("Spalart–Allmaras: ν̃ needs one finite value per cell".into()));
        }
        for ((out, &v), &f) in self.nu_tilde.iter_mut().zip(values).zip(&self.fluid) {
            *out = if f { v } else { 0.0 };
        }
        Ok(())
    }

    /// Fix ν̃ at the given fluid cells, `(cell, ν̃)`: the wall function's first-node values. A
    /// cell outside the field or not fluid, or a value that is not finite, is refused and nothing
    /// changes.
    pub fn impose(&mut self, values: &[(usize, f64)]) -> crate::Result<()> {
        if let Some((c, v)) = values.iter().find(|(c, v)| !(v.is_finite() && self.fluid.get(*c).copied().unwrap_or(false))) {
            return Err(crate::Error::Invalid(format!("Spalart–Allmaras: ν̃ {v} at cell {c}, which is not a fluid cell of the field or not finite")));
        }
        for &(c, v) in values {
            if let Some(out) = self.nu_tilde.get_mut(c) {
                *out = v;
            }
        }
        Ok(())
    }

    /// The eddy viscosity ν_t = max(0, ν̃ f_v1) of `cell` (0 outside the field).
    pub fn nu_t(&self, cell: usize, nu: f64) -> f64 {
        nu_t(self.nu_tilde.get(cell).copied().unwrap_or(0.0), nu)
    }

    /// The bytes its arrays hold (their capacities).
    pub fn heap_bytes(&self) -> u64 {
        (self.nu_tilde.capacity() * size_of::<f64>()
            + self.sides.capacity() * size_of::<[Side; 6]>()
            + self.fluid.capacity()
            + self.rates.capacity() * size_of::<(f64, f64)>()
            + self.next.capacity() * size_of::<f64>()) as u64
    }

    /// Save or restore ν̃ (the faces, fluid cells and inflow value are fixed by the setup and
    /// must match). A restored ν̃ must be finite, and 0 outside the fluid cells.
    pub(crate) fn archive(&mut self, a: &mut dyn crate::checkpoint::Archive) -> crate::Result<()> {
        crate::checkpoint::fixed(a, self.inflow.to_bits())?;
        let n = self.fluid.len();
        crate::checkpoint::doubles(a, &mut self.nu_tilde, n)?;
        if self.nu_tilde.iter().zip(&self.fluid).any(|(v, &f)| !v.is_finite() || (!f && *v != 0.0)) {
            return Err(crate::Error::Invalid("checkpoint Spalart–Allmaras ν̃ is not finite, or not 0 outside the fluid".into()));
        }
        Ok(())
    }

    /// The transport rate of cell `c` (advection, diffusion, c_b2 term, extra source) and the sum
    /// of its neighbours' coefficients (the bound on a monotone sub-step is its inverse).
    fn rate(&self, c: usize, input: &SaInput) -> (f64, f64) {
        let (Some(sides), Some(&v)) = (self.sides.get(c), self.nu_tilde.get(c)) else {
            return (0.0, 0.0);
        };
        let nu = input.nu;
        let uc = input.u.get(c).copied().unwrap_or([0.0; 3]);
        let eta_c = diffusivity(v, nu);
        let (mut rate, mut coef) = (0.0, 0.0);
        let mut slopes = [[None::<f64>; 2]; 3];
        for (k, side) in sides.iter().enumerate() {
            let a = AXIS_OF_SIDE[k];
            let sign = SIGN_OF_SIDE[k];
            // Face value, outward normal velocity, diffusion weight and slope per half spacing.
            let (face, un, d) = match *side {
                Side::Cell(n) => {
                    let n = n as usize;
                    let vn = self.nu_tilde.get(n).copied().unwrap_or(0.0);
                    let un = input.u.get(n).map_or(0.0, |un| 0.5 * (uc[a] + un[a]) * sign);
                    let d = 0.5 * (eta_c + diffusivity(vn, nu)) / SIGMA;
                    slopes[a][k % 2] = Some(vn - v);
                    (Some(vn), un, d)
                }
                Side::Wall => {
                    slopes[a][k % 2] = Some(-2.0 * v);
                    (Some(0.0), 0.0, 2.0 * nu / SIGMA)
                }
                Side::Node => {
                    slopes[a][k % 2] = Some(-v);
                    (Some(0.0), 0.0, 0.5 * (eta_c + nu) / SIGMA)
                }
                Side::Inflow => {
                    slopes[a][k % 2] = Some(2.0 * (self.inflow - v));
                    (Some(self.inflow), uc[a] * sign, 2.0 * diffusivity(self.inflow, nu) / SIGMA)
                }
                Side::Open => (None, uc[a] * sign, 0.0),
            };
            let upwind = (-un).max(0.0);
            if let Some(face) = face {
                rate += (upwind + d) * (face - v);
                coef += upwind + d;
            }
        }
        let mut grad2 = 0.0;
        for s in slopes {
            grad2 += match s {
                [Some(a), Some(b)] => 0.5 * (a * a + b * b),
                [Some(a), None] | [None, Some(a)] => a * a,
                [None, None] => 0.0,
            };
        }
        rate += C_B2 / SIGMA * grad2;
        if let Some(s) = input.source {
            rate += s.get(c).copied().unwrap_or(0.0);
        }
        (rate, coef)
    }

    /// One lattice step (Δt = 1) of ν̃, in as many sub-steps as monotone transport needs.
    ///
    /// An `Err` (bad input lengths, or [`Error::Unstable`](crate::Error::Unstable) after
    /// [`MAX_SUBSTEPS`]) may leave ν̃ partly advanced: the caller discards the field, as
    /// `Session::execute` restores the model on any error.
    pub fn advance(&mut self, input: &SaInput) -> crate::Result<SaStats> {
        let cells = self.nu_tilde.len();
        let bad = |what: &str, len: usize| crate::Error::Invalid(format!("Spalart–Allmaras: {len} {what} for {cells} cells"));
        if input.u.len() != cells {
            return Err(bad("velocities", input.u.len()));
        }
        if input.vorticity.len() != cells {
            return Err(bad("vorticities", input.vorticity.len()));
        }
        if input.length.len() != cells {
            return Err(bad("lengths", input.length.len()));
        }
        if let Some(s) = input.source
            && s.len() != cells
        {
            return Err(bad("sources", s.len()));
        }
        if !(input.nu > 0.0 && input.nu.is_finite()) {
            return Err(crate::Error::Invalid(format!("Spalart–Allmaras: the viscosity {} must be above 0", input.nu)));
        }
        if let Some(n) = input.substeps.filter(|n| !(1..=MAX_SUBSTEPS).contains(n)) {
            return Err(crate::Error::Invalid(format!("Spalart–Allmaras: {n} sub-steps; give 1 to {MAX_SUBSTEPS}")));
        }
        let mut rates = std::mem::take(&mut self.rates);
        let mut next = std::mem::take(&mut self.next);
        rates.resize(cells, (0.0, 0.0));
        next.resize(cells, 0.0);
        let result = self.substeps(input, &mut rates, &mut next);
        self.rates = rates;
        self.next = next;
        result
    }

    fn substeps(&mut self, input: &SaInput, rates: &mut [(f64, f64)], next: &mut Vec<f64>) -> crate::Result<SaStats> {
        let too_many = |field: &[f64]| {
            crate::Error::Unstable(format!(
                "Spalart–Allmaras: ν̃'s transport needs more than {MAX_SUBSTEPS} sub-steps per lattice step (the largest ν̃/ν is {:.3e}): refine the grid or lower the inflow ν̃",
                field.iter().fold(0.0f64, |m, &v| m.max(v)) / input.nu
            ))
        };
        let mut remaining = 1.0f64;
        let mut substeps = 0;
        while remaining > 0.0 {
            if substeps == MAX_SUBSTEPS {
                return Err(too_many(&self.nu_tilde));
            }
            let this = &*self;
            fieldcraft_parallel::chunks_mut(rates, BLOCK, |b, out| {
                for (i, r) in out.iter_mut().enumerate() {
                    let c = b * BLOCK + i;
                    *r = if this.fluid.get(c).copied().unwrap_or(false) { this.rate(c, input) } else { (0.0, 0.0) };
                }
            });
            let most = rates.iter().fold(0.0f64, |m, r| if r.1 > m { r.1 } else { m });
            if !most.is_finite() {
                return Err(crate::Error::Unstable("Spalart–Allmaras: ν̃'s transport coefficients are not finite".into()));
            }
            if most > f64::from(MAX_SUBSTEPS) {
                return Err(too_many(&self.nu_tilde));
            }
            // A fixed count takes 1/n each time, n times.
            let parts = match input.substeps {
                None => (most * remaining).ceil().max(1.0),
                Some(n) if most <= f64::from(n) => f64::from(n - substeps),
                Some(n) => {
                    return Err(crate::Error::Unstable(format!(
                        "Spalart–Allmaras: {n} fixed sub-steps break ν̃'s monotone bound, which needs {}",
                        most.ceil()
                    )));
                }
            };
            let dt = if parts <= 1.0 { remaining } else { remaining / parts };
            let this = &*self;
            let rates = &*rates;
            fieldcraft_parallel::chunks_mut(next, BLOCK, |b, out| {
                for (i, o) in out.iter_mut().enumerate() {
                    let c = b * BLOCK + i;
                    *o = this.update(c, rates.get(c).map_or(0.0, |r| r.0), dt, input);
                }
            });
            if let Some(c) = next.iter().position(|v| !v.is_finite()) {
                return Err(crate::Error::Unstable(format!("Spalart–Allmaras: ν̃ is not finite at cell {c}")));
            }
            std::mem::swap(&mut self.nu_tilde, next);
            remaining = if parts <= 1.0 { 0.0 } else { remaining - dt };
            substeps += 1;
        }
        let mut stats = SaStats { substeps, min_nu_tilde: f64::INFINITY, ..SaStats::default() };
        for (&v, &f) in self.nu_tilde.iter().zip(&self.fluid) {
            if f {
                stats.max_nu_t = stats.max_nu_t.max(nu_t(v, input.nu));
                stats.min_nu_tilde = stats.min_nu_tilde.min(v);
                stats.negative_cells += usize::from(v < 0.0);
            }
        }
        if stats.min_nu_tilde == f64::INFINITY {
            stats.min_nu_tilde = 0.0;
        }
        Ok(stats)
    }

    /// Cell `c` after a sub-step `dt` with transport rate `rate`: explicit transport, then the
    /// model's sources point-implicitly.
    fn update(&self, c: usize, rate: f64, dt: f64, input: &SaInput) -> f64 {
        if !self.fluid.get(c).copied().unwrap_or(false) {
            return 0.0;
        }
        let v = self.nu_tilde.get(c).copied().unwrap_or(0.0) + dt * rate;
        if !input.model {
            return v;
        }
        let omega = input.vorticity.get(c).copied().unwrap_or(0.0);
        let d = input.length.get(c).copied().unwrap_or(0.0);
        let Some((p, k)) = coefficients(v, input.nu, omega, d) else { return v };
        if v > 0.0 {
            // D = k ν̃ is implicit, below 0 only through f_t2, then explicit; P is below 0 where
            // f_t2 > 1 (χ below 0.604), then a sink, implicit too.
            (v + dt * v * (p.max(0.0) + (-k).max(0.0))) / (1.0 + dt * (k.max(0.0) + (-p).max(0.0)))
        } else if v < 0.0 {
            // Both SA-neg sources are a decay towards 0: rate = −a ν̃ with a = k − p ≥ 0
            // (k = −c_w1 ν̃/d² ≥ 0 and p = c_b1 (1 − c_t3) Ω ≤ 0).
            v / (1.0 + dt * (k - p).max(0.0))
        } else {
            v
        }
    }
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

#[cfg(test)]
mod transport_tests {
    use super::*;
    use std::f64::consts::PI;

    /// The manufactured field (physical units on the unit square, periodic):
    /// ν̃ = ν (A + B cos(ωt) sin(kx) cos(ky)), k = ω = 2π, frozen velocity
    /// u = U (0.3 + 0.2 sin(ky), −0.2 + 0.1 cos(kx)). The forcing, derived by hand:
    /// S = ∂ν̃/∂t + u·∇ν̃ − (1/σ)[η′(ν̃)|∇ν̃|² + η(ν̃)∇²ν̃ + c_b2 |∇ν̃|²] (∇·(η∇ν̃) expanded),
    /// with ∇²ν̃ = −2k² (ν̃ − νA), η′ = 1 above 0 and f_n + χ f_n′(χ), f_n′ = 96 χ² / (16 − χ³)²,
    /// below.
    struct Mms {
        nu: f64,
        a: f64,
        b: f64,
        u: f64,
    }

    const K: f64 = 2.0 * PI;
    const W: f64 = 2.0 * PI;

    impl Mms {
        fn value(&self, x: f64, y: f64, t: f64) -> f64 {
            self.nu * (self.a + self.b * (W * t).cos() * (K * x).sin() * (K * y).cos())
        }

        fn velocity(&self, x: f64, y: f64) -> [f64; 2] {
            [self.u * (0.3 + 0.2 * (K * y).sin()), self.u * (-0.2 + 0.1 * (K * x).cos())]
        }

        fn eta_prime(&self, v: f64) -> f64 {
            if v >= 0.0 {
                return 1.0;
            }
            let chi = v / self.nu;
            let c3 = chi * chi * chi;
            f_n(chi) + chi * 96.0 * chi * chi / ((16.0 - c3) * (16.0 - c3))
        }

        fn forcing(&self, x: f64, y: f64, t: f64) -> f64 {
            let (nu, b) = (self.nu, self.b);
            let v = self.value(x, y, t);
            let dt = -nu * b * W * (W * t).sin() * (K * x).sin() * (K * y).cos();
            let gx = nu * b * (W * t).cos() * K * (K * x).cos() * (K * y).cos();
            let gy = -nu * b * (W * t).cos() * K * (K * x).sin() * (K * y).sin();
            let lap = -2.0 * K * K * (v - nu * self.a);
            let [ux, uy] = self.velocity(x, y);
            let g2 = gx * gx + gy * gy;
            dt + ux * gx + uy * gy - (self.eta_prime(v) * g2 + diffusivity(v, nu) * lap + C_B2 * g2) / SIGMA
        }

        /// The L2 error (over ν) at `t_end` on an n × n grid, in lattice units with Δt = h²
        /// (diffusive scaling: ν and ν̃ keep their values, velocities scale by h, rates by h²).
        fn error(&self, n: usize, t_end: f64) -> f64 {
            let h = 1.0 / n as f64;
            let cells = n * n;
            let at = |c: usize| (((c % n) as f64 + 0.5) * h, ((c / n) as f64 + 0.5) * h);
            let mut field = SaField::boxed([n, n, 1], &vec![true; cells], [BoxFace::Periodic; 6], 0.0, 0.0).unwrap();
            let init: Vec<f64> = (0..cells).map(|c| self.value(at(c).0, at(c).1, 0.0)).collect();
            field.set_nu_tilde(&init).unwrap();
            let u: Vec<[f64; 3]> = (0..cells)
                .map(|c| {
                    let [ux, uy] = self.velocity(at(c).0, at(c).1);
                    [ux * h, uy * h, 0.0]
                })
                .collect();
            let zeros = vec![0.0; cells];
            let steps = (t_end / (h * h)).round() as usize;
            for s in 0..steps {
                let t = s as f64 * h * h;
                let source: Vec<f64> = (0..cells).map(|c| self.forcing(at(c).0, at(c).1, t) * h * h).collect();
                let input = SaInput { u: &u, vorticity: &zeros, length: &zeros, nu: self.nu, model: false, source: Some(&source), substeps: None };
                field.advance(&input).unwrap();
            }
            let t = steps as f64 * h * h;
            let sum: f64 = (0..cells).map(|c| (field.nu_tilde()[c] - self.value(at(c).0, at(c).1, t)).powi(2)).sum();
            (sum / cells as f64).sqrt() / self.nu
        }

        /// The observed orders 16→32 and 32→64, and the three errors.
        fn orders(&self, t_end: f64) -> (f64, f64, [f64; 3]) {
            let e = [16, 32, 64].map(|n| self.error(n, t_end));
            let (p1, p2) = ((e[0] / e[1]).log2(), (e[1] / e[2]).log2());
            println!("MMS A {} U {} T {t_end}: L2 errors {e:?}, orders {p1:.4} {p2:.4}", self.a, self.u);
            (p1, p2, e)
        }
    }

    /// Advection, diffusion and source: first order (upwind), in [0.9, 1.2] on 16/32/64. The
    /// advection is strong (U = 4) so its first-order error dominates the second-order ones, over
    /// a short time (T = 0.025) so it stays in its linear range (k² U h T / 2 small): with weaker
    /// advection the two cancel on these grids, with longer times the damping saturates.
    #[test]
    fn mms_advection_diffusion_source_is_first_order() {
        let (p1, p2, e) = Mms { nu: 0.02, a: 1.5, b: 1.0, u: 4.0 }.orders(0.025);
        assert!((0.9..=1.2).contains(&p2), "orders {p1:.3} {p2:.3}, errors {e:?}");
    }

    /// Pure diffusion (the c_b2 term and the ν̃-dependent diffusivity): second order, in [1.8, 2.2].
    #[test]
    fn mms_pure_diffusion_is_second_order() {
        let (p1, p2, e) = Mms { nu: 0.02, a: 1.5, b: 1.0, u: 0.0 }.orders(0.25);
        assert!((1.8..=2.2).contains(&p2), "orders {p1:.3} {p2:.3}, errors {e:?}");
    }

    /// A field that crosses zero (SA-neg's diffusivity below 0) converges at the same orders.
    #[test]
    fn mms_through_zero_keeps_its_orders() {
        let (p1, p2, e) = Mms { nu: 0.02, a: 0.2, b: 1.0, u: 4.0 }.orders(0.025);
        assert!((0.9..=1.2).contains(&p2), "advected: orders {p1:.3} {p2:.3}, errors {e:?}");
        let (p1, p2, e) = Mms { nu: 0.02, a: 0.2, b: 1.0, u: 0.0 }.orders(0.25);
        assert!((1.8..=2.2).contains(&p2), "diffused: orders {p1:.3} {p2:.3}, errors {e:?}");
    }

    const NU: f64 = 1e-3;

    /// A channel-like box: walls at y = −½ and y = 15½, inflow at −x, open at +x, periodic in z,
    /// a block obstacle, a sheared velocity and a random initial ν̃.
    fn channel(seed: f64) -> (SaField, Vec<[f64; 3]>, Vec<f64>, Vec<f64>) {
        let dims = [24, 16, 12];
        let cells = 24 * 16 * 12;
        let fluid: Vec<bool> = (0..cells)
            .map(|c| {
                let (i, j) = (c % 24, (c / 24) % 16);
                !((10..13).contains(&i) && (6..9).contains(&j))
            })
            .collect();
        let faces = [BoxFace::Inflow, BoxFace::Open, BoxFace::Wall, BoxFace::Wall, BoxFace::Periodic, BoxFace::Periodic];
        let mut field = SaField::boxed(dims, &fluid, faces, 3.0 * NU, 3.0 * NU).unwrap();
        let mut x = seed;
        let init: Vec<f64> = (0..cells)
            .map(|_| {
                x = (x * 9301.0 + 0.4927).fract();
                x * 50.0 * NU
            })
            .collect();
        field.set_nu_tilde(&init).unwrap();
        let u: Vec<[f64; 3]> = (0..cells)
            .map(|c| {
                let y = (((c / 24) % 16) as f64 + 0.5) / 16.0;
                [0.15 * 4.0 * y * (1.0 - y), 0.01 * ((c % 7) as f64 - 3.0), 0.0]
            })
            .collect();
        let vort: Vec<f64> = (0..cells).map(|c| 0.6 * (1.0 - 2.0 * (((c / 24) % 16) as f64 + 0.5) / 16.0).abs() / 16.0).collect();
        // The exact distance to every wall face (the channel's and the obstacle's), as the
        // solver's own field gives it: the destruction needs it, the plane walls alone would
        // leave the cells behind the obstacle with d ≈ 7.5 and ν̃ running away there.
        let at = |c: usize| [c % 24, (c / 24) % 16, c / (24 * 16)].map(|v| v as f64);
        let mut points = Vec::new();
        for (c, sides) in field.sides.iter().enumerate() {
            for (k, side) in sides.iter().enumerate() {
                if *side == Side::Wall && fluid[c] {
                    let mut p = at(c);
                    p[AXIS_OF_SIDE[k]] += 0.5 * SIGN_OF_SIDE[k];
                    points.push(p);
                }
            }
        }
        let d: Vec<f64> = (0..cells)
            .map(|c| {
                let p = at(c);
                points
                    .iter()
                    .map(|q: &[f64; 3]| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt())
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        assert_eq!(d[24 * 7 + 13], 0.5);
        (field, u, vort, d)
    }

    fn input<'a>(u: &'a [[f64; 3]], vort: &'a [f64], d: &'a [f64]) -> SaInput<'a> {
        SaInput { u, vorticity: vort, length: d, nu: NU, model: true, source: None, substeps: None }
    }

    /// From a non-negative start, with the model on, ν̃ never goes below 0 and ν_t ≥ 0, every step.
    #[test]
    fn nu_tilde_stays_non_negative() {
        let (mut field, u, vort, d) = channel(0.37);
        for step in 0..300 {
            let stats = field.advance(&input(&u, &vort, &d)).unwrap();
            assert_eq!(stats.negative_cells, 0, "step {step}");
            assert!(stats.min_nu_tilde >= 0.0 && stats.substeps >= 1);
            assert!(field.nu_tilde().iter().all(|&v| v >= 0.0 && v.is_finite()), "step {step}");
            assert!((0..field.nu_tilde().len()).all(|c| field.nu_t(c, NU) >= 0.0));
        }
    }

    /// From a field with negative cells (SA-neg's branch), the most negative ν̃ only rises towards
    /// 0 and, once there, stays at or above it; ν_t is 0 wherever ν̃ ≤ 0.
    #[test]
    fn negative_nu_tilde_recovers_and_has_no_eddy_viscosity() {
        let (mut field, u, vort, d) = channel(0.61);
        let shifted: Vec<f64> = field.nu_tilde().iter().map(|v| v - 25.0 * NU).collect();
        field.set_nu_tilde(&shifted).unwrap();
        let mut lowest = field.nu_tilde().iter().fold(0.0f64, |m, &v| m.min(v));
        assert!(lowest < 0.0);
        for step in 0..300 {
            let stats = field.advance(&input(&u, &vort, &d)).unwrap();
            assert!(stats.min_nu_tilde >= lowest.min(0.0), "step {step}: {} below {lowest}", stats.min_nu_tilde);
            lowest = stats.min_nu_tilde;
            for (c, &v) in field.nu_tilde().iter().enumerate() {
                assert!(field.nu_t(c, NU) >= 0.0);
                if v <= 0.0 {
                    assert_eq!(field.nu_t(c, NU), 0.0);
                }
            }
        }
    }

    /// The same result bit for bit on 1, 2 and 7 threads (the box has 4608 cells, two blocks).
    #[test]
    fn threads_do_not_change_the_result() {
        let run = |threads: usize| {
            fieldcraft_parallel::with_threads(threads, || {
                let (mut field, u, vort, d) = channel(0.37);
                for _ in 0..40 {
                    field.advance(&input(&u, &vort, &d)).unwrap();
                }
                field.nu_tilde().iter().map(|v| v.to_bits()).collect::<Vec<u64>>()
            })
        };
        let one = run(1);
        assert_eq!(one, run(2));
        assert_eq!(one, run(7));
    }

    /// The inverse of ν̃ f_v1 returns ν̃ whose ν_t is the target to round-off over eighteen decades
    /// of ν_t/ν, and 0 for targets that are not positive finite numbers.
    #[test]
    fn nu_tilde_for_inverts_the_eddy_viscosity() {
        for k in -30..=60 {
            let t = 10f64.powf(k as f64 / 5.0);
            let v = nu_tilde_for(t * NU, NU).unwrap();
            let back = nu_t(v, NU);
            assert!((back / (t * NU) - 1.0).abs() < 1e-13, "ν_t/ν {t}: {back} vs {}", t * NU);
        }
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(nu_tilde_for(bad, NU).unwrap(), 0.0);
            assert_eq!(nu_tilde_for(NU, bad).unwrap(), 0.0);
        }
        let (mut field, _, _, _) = channel(0.37);
        assert!(field.impose(&[(24 * 7 + 11, f64::NAN)]).is_err());
        assert!(field.impose(&[(24 * 7 + 11, 1.0)]).is_err(), "a solid cell of the obstacle");
        assert!(field.impose(&[(usize::MAX, 1.0)]).is_err());
        field.impose(&[(0, 0.25)]).unwrap();
        assert_eq!(field.nu_tilde()[0], 0.25);
    }

    /// The wall laws' eddy viscosity: Spalding's ν_t⁺ = dy⁺/du⁺ − 1 matches a centred difference
    /// of its y⁺(u⁺); Musker's is the wall model's wall_law; Spalding's is ≈ 0 deep in the
    /// sublayer and both rise with y⁺; hostile y⁺ give 0.
    #[test]
    fn the_wall_laws_eddy_viscosity() {
        let (k, b) = (crate::wall::BOUNCE_KAPPA, crate::wall::BOUNCE_B);
        let y_of = |u: f64| {
            let x = k * u;
            u + (-k * b).exp() * (x.exp() - 1.0 - x - x * x / 2.0 - x * x * x / 6.0)
        };
        for y_plus in [1.0, 5.0, 11.0, 30.0, 100.0, 1000.0] {
            let s = law_eddy_viscosity(WallLaw::Spalding, y_plus).unwrap();
            let u = bisect(0.0, y_plus, |u| y_of(u) < y_plus, "test").unwrap();
            assert!((y_of(u) / y_plus - 1.0).abs() < 1e-12);
            let h = 1e-5 * u;
            let fd = (y_of(u + h) - y_of(u - h)) / (2.0 * h) - 1.0;
            assert!((s - fd).abs() < 1e-6 * (1.0 + fd.abs()), "Spalding y⁺ {y_plus}: {s} vs {fd}");
            // Musker's is the wall model's own: wall_law at u_τ = 1, ν = 1 (so y = y⁺).
            let m = law_eddy_viscosity(WallLaw::Musker, y_plus).unwrap();
            let own = crate::wall::wall_law(crate::wall::musker(3.0 * y_plus).0, y_plus, 3.0 * y_plus, 1.0).1;
            assert!((m / own - 1.0).abs() < 1e-9, "Musker y⁺ {y_plus}: {m} vs {own}");
        }
        assert!(law_eddy_viscosity(WallLaw::Spalding, 1.0).unwrap().abs() < 1e-3);
        for bad in [0.0, -3.0, f64::NAN, f64::INFINITY] {
            for law in [WallLaw::Spalding, WallLaw::Musker] {
                assert_eq!(law_eddy_viscosity(law, bad).unwrap(), 0.0);
            }
        }
    }

    /// The first node's ν̃ against SA's log-layer value κ u_τ y: frozen from an independent
    /// evaluation (Spalding κ 0.41, B 5.5; Musker's mixing length as wall.rs writes it, κ 0.384, A⁺ 26), ν̃/(κ u_τ y) − 1
    /// within ±0.005 of the table at each y⁺.
    #[test]
    fn the_first_node_nu_tilde_against_kappa_u_tau_y_is_frozen() {
        let y = [30.0, 40.0, 50.0, 70.0, 100.0, 200.0, 500.0];
        let spalding = [-0.1443, -0.1648, -0.1618, -0.1376, -0.1045, -0.0538, -0.0201];
        let musker = [-0.2458, -0.2580, -0.2402, -0.1892, -0.1458, -0.1284, -0.1288];
        for (law, table) in [(WallLaw::Spalding, spalding), (WallLaw::Musker, musker)] {
            for (&y_plus, &want) in y.iter().zip(&table) {
                // Wall units: ν = 1, u_τ = 1, so y = y⁺ and ν_t = ν_t⁺.
                let nu_t_law = law_eddy_viscosity(law, y_plus).unwrap();
                let v = nu_tilde_for(nu_t_law, 1.0).unwrap();
                assert!((nu_t(v, 1.0) / nu_t_law - 1.0).abs() < 1e-12);
                let got = v / (KAPPA * y_plus) - 1.0;
                println!("{law:?} y⁺ {y_plus}: ν_t⁺ {nu_t_law:.4}, ν̃/ν {v:.4}, ν̃/(κ u_τ y) − 1 {got:+.4} (frozen {want:+.4})");
                assert!((got - want).abs() <= 0.005, "{law:?} y⁺ {y_plus}: {got} vs {want}");
            }
        }
    }

    /// The on-node wall (Han, Ooka & Kikumoto's boundary node, a full cell from the first fluid
    /// node) is a neighbour held at ν̃ = 0 and at rest: each cell's rate in a column between two
    /// such walls is, bit for bit, its rate in the same column with the wall nodes as cells at 0.
    /// The halfway wall, half a cell away, gives another rate (twice the slope).
    #[test]
    fn the_on_node_wall_is_a_zero_neighbour_a_full_cell_away() {
        let open = Side::Open;
        let col = |below: Side, above: Side| [open, open, below, above, open, open];
        let walled = |wall: Side| {
            SaField::new(vec![col(wall, Side::Cell(1)), col(Side::Cell(0), Side::Cell(2)), col(Side::Cell(1), wall)], vec![true; 3], 0.0, 0.0)
                .unwrap()
        };
        let mut on_node = walled(Side::Node);
        let mut halfway = walled(Side::Wall);
        let mut padded = SaField::new(
            vec![
                col(open, Side::Cell(1)),
                col(Side::Cell(0), Side::Cell(2)),
                col(Side::Cell(1), Side::Cell(3)),
                col(Side::Cell(2), Side::Cell(4)),
                col(Side::Cell(3), open),
            ],
            vec![true; 5],
            0.0,
            0.0,
        )
        .unwrap();
        let inner = [2e-3, 5e-3, 1e-3];
        on_node.set_nu_tilde(&inner).unwrap();
        halfway.set_nu_tilde(&inner).unwrap();
        padded.set_nu_tilde(&[0.0, inner[0], inner[1], inner[2], 0.0]).unwrap();
        let (u3, u5) = (vec![[0.1, 0.0, 0.0]; 3], vec![[0.1, 0.0, 0.0]; 5]);
        let (z3, z5) = (vec![0.0; 3], vec![0.0; 5]);
        let at = |u, z| SaInput { u, vorticity: z, length: z, nu: NU, model: true, source: None, substeps: None };
        for c in 0..3 {
            let (node, cell) = (on_node.rate(c, &at(&u3, &z3)), padded.rate(c + 1, &at(&u5, &z5)));
            assert_eq!((node.0.to_bits(), node.1.to_bits()), (cell.0.to_bits(), cell.1.to_bits()), "cell {c}");
        }
        assert_ne!(halfway.rate(0, &at(&u3, &z3)), on_node.rate(0, &at(&u3, &z3)));
    }

    /// Where f_t2 > 1 (χ below 0.604) the production is a sink, and a long step takes it
    /// implicitly: ν̃ decays towards 0 and stays positive.
    #[test]
    fn a_production_below_zero_is_taken_implicitly() {
        let mut field = SaField::new(vec![[Side::Open; 6]], vec![true], 0.0, 0.0).unwrap();
        let v = 0.1 * NU;
        field.set_nu_tilde(&[v]).unwrap();
        let (u, vort, d) = (vec![[0.0; 3]], vec![100.0], vec![1.0]);
        let input = SaInput { u: &u, vorticity: &vort, length: &d, nu: NU, model: true, source: None, substeps: None };
        let (p, _) = coefficients(v, NU, 100.0, 1.0).unwrap();
        assert!(p < -1.0, "{p}");
        let next = field.update(0, 0.0, 10.0, &input);
        assert!(next > 0.0 && next < v, "{next} from {v}");
    }

    /// A fixed count of n sub-steps takes exactly n (with n the adaptive count, the same field bit
    /// for bit: both take 1/n each time); a count below the monotone bound, 0 or above the cap is
    /// refused.
    #[test]
    fn a_fixed_sub_step_count_is_taken_or_refused() {
        let (field, u, vort, d) = channel(0.37);
        let bits = |f: &SaField| f.nu_tilde().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        let mut adaptive = field.clone();
        let taken = adaptive.advance(&input(&u, &vort, &d)).unwrap().substeps;
        let mut fixed = field.clone();
        assert_eq!(fixed.advance(&SaInput { substeps: Some(taken), ..input(&u, &vort, &d) }).unwrap().substeps, taken);
        assert_eq!(bits(&fixed), bits(&adaptive));
        let mut more = field.clone();
        assert_eq!(more.advance(&SaInput { substeps: Some(taken + 3), ..input(&u, &vort, &d) }).unwrap().substeps, taken + 3);
        let mut hot = field.clone();
        hot.set_nu_tilde(&vec![0.2; u.len()]).unwrap();
        let needed = hot.clone().advance(&input(&u, &vort, &d)).unwrap().substeps;
        println!("fixed sub-steps: the channel takes {taken}, the hot field {needed}");
        assert!(needed > 1);
        assert!(matches!(hot.advance(&SaInput { substeps: Some(1), ..input(&u, &vort, &d) }), Err(crate::Error::Unstable(_))));
        for n in [0, MAX_SUBSTEPS + 1] {
            let mut f = field.clone();
            assert!(matches!(f.advance(&SaInput { substeps: Some(n), ..input(&u, &vort, &d) }), Err(crate::Error::Invalid(_))), "{n}");
        }
    }

    /// Hostile inputs are refused, not panicked on; a runaway diffusivity stops at the sub-step cap.
    #[test]
    fn hostile_inputs_are_refused() {
        let (mut field, u, vort, d) = channel(0.37);
        let short = vec![0.0; 3];
        assert!(matches!(field.advance(&input(&u, &short, &d)), Err(crate::Error::Invalid(_))));
        for nu in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let bad = SaInput { nu, ..input(&u, &vort, &d) };
            assert!(field.advance(&bad).is_err(), "{nu}");
        }
        assert!(field.set_nu_tilde(&vec![f64::NAN; u.len()]).is_err());
        field.set_nu_tilde(&vec![1e3; u.len()]).unwrap();
        assert!(matches!(field.advance(&input(&u, &vort, &d)), Err(crate::Error::Unstable(_))));
        assert!(SaField::boxed([2, 2, 2], &[true; 7], [BoxFace::Open; 6], 0.0, 0.0).is_err());
        let one_sided = [BoxFace::Periodic, BoxFace::Open, BoxFace::Open, BoxFace::Open, BoxFace::Open, BoxFace::Open];
        assert!(SaField::boxed([2, 2, 2], &[true; 8], one_sided, 0.0, 0.0).is_err());
        assert!(SaField::new(vec![[Side::Cell(5); 6]], vec![true], 0.0, 0.0).is_err());
        assert!(SaField::new(vec![[Side::Open; 6]], vec![true], f64::NAN, 0.0).is_err());
    }
}
