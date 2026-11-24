//! Helpers for the Spalart–Allmaras channel cases: the integral of the wall law over an on-node
//! wall's interval, the settling rule, and a gate table that prints each result against its
//! reference and tolerance.

/// The unresolved interval between an on-node wall (Han et al.'s `spalding_bounce`, the wall at
/// the boundary node's centre) and the first fluid node at distance `y1`, whose speed parallel to
/// the wall is `u1`: its integrals `∫₀^{y1} u dy` and `∫₀^{y1} u² dy`, from the wall law the wall
/// applies (Spalding's, κ and B of [`fieldcraft_lbm::wall::BOUNCE_KAPPA`] and `BOUNCE_B`), and
/// from the trapezoid between u = 0 at the wall and u1 (what a sum over the rows with the wall
/// row at u = 0 gives). Every reading that integrates a profile from the wall reports both, and
/// gates on the law's: the trapezoid misses the curvature of the law over the wall interval.
///
/// Spalding's law gives y⁺(u⁺) explicitly, so the integrals are closed forms by parts:
/// `∫₀^{y⁺₁} u⁺ dy⁺ = u⁺₁ y⁺₁ − ∫₀^{u⁺₁} y⁺ du⁺` and `∫₀^{y⁺₁} u⁺² dy⁺ = u⁺₁² y⁺₁ − 2 ∫₀^{u⁺₁} u⁺ y⁺ du⁺`;
/// in physical units `∫ u dy = ν ∫ u⁺ dy⁺` and `∫ u² dy = ν u_τ ∫ u⁺² dy⁺`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallInterval {
    /// The friction velocity the law gives at the first node (what the wall applies).
    pub u_tau: f64,
    pub y1: f64,
    /// `[∫ u dy, ∫ u² dy]` from the wall law.
    pub law: [f64; 2],
    /// The same from the trapezoid, `[u1 y1 / 2, u1² y1 / 2]`.
    pub trapezoid: [f64; 2],
}

impl WallInterval {
    /// `None` for a speed, distance or viscosity that is not a positive finite number.
    pub fn new(u1: f64, y1: f64, nu: f64) -> Option<WallInterval> {
        use fieldcraft_lbm::wall::{BOUNCE_B, BOUNCE_KAPPA, spalding_friction};
        if !(u1.is_finite() && y1.is_finite() && nu.is_finite() && u1 > 0.0 && y1 > 0.0 && nu > 0.0) {
            return None;
        }
        let u_tau = spalding_friction(u1, y1, nu, BOUNCE_KAPPA, BOUNCE_B);
        if !(u_tau > 0.0) {
            return None;
        }
        let (k, e) = (BOUNCE_KAPPA, (-BOUNCE_KAPPA * BOUNCE_B).exp());
        let u = u1 / u_tau;
        let ku = k * u;
        let y_plus = u + e * (ku.exp() - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
        // ∫₀^u y⁺ du⁺ and ∫₀^u u⁺ y⁺ du⁺.
        let y1_int = u * u / 2.0 + e * ((ku.exp() - 1.0) / k - u - k * u * u / 2.0 - k * k * u.powi(3) / 6.0 - k.powi(3) * u.powi(4) / 24.0);
        let y2_int = u.powi(3) / 3.0
            + e * (ku.exp() * (u / k - 1.0 / (k * k)) + 1.0 / (k * k)
                - u * u / 2.0
                - k * u.powi(3) / 3.0
                - k * k * u.powi(4) / 8.0
                - k.powi(3) * u.powi(5) / 30.0);
        let first = u * y_plus - y1_int;
        let second = u * u * y_plus - 2.0 * y2_int;
        Some(WallInterval { u_tau, y1, law: [nu * first, nu * u_tau * second], trapezoid: [0.5 * u1 * y1, 0.5 * u1 * u1 * y1] })
    }

    /// [`WallInterval::new`] with the Spalart–Allmaras constant-stress law: the friction velocity
    /// from [`fieldcraft_lbm::sa_law::friction_velocity`] and the integrals `∫ u⁺ dy⁺`, `∫ u⁺² dy⁺`
    /// of the tabulated u⁺ by Simpson's rule on 2000 uniform panels (u⁺ is smooth on [0, y⁺₁]).
    pub fn new_sa_law(u1: f64, y1: f64, nu: f64) -> Option<WallInterval> {
        use fieldcraft_lbm::sa_law::{friction_velocity, u_plus};
        if !(u1.is_finite() && y1.is_finite() && nu.is_finite() && u1 > 0.0 && y1 > 0.0 && nu > 0.0) {
            return None;
        }
        let u_tau = friction_velocity(u1, y1, nu);
        if !(u_tau > 0.0) {
            return None;
        }
        let y_plus = y1 * u_tau / nu;
        let panels = 2000;
        let h = y_plus / panels as f64;
        let (mut first, mut second) = (0.0, 0.0);
        for k in 0..=panels {
            let w = if k == 0 || k == panels {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let u = u_plus(k as f64 * h);
            first += w * u;
            second += w * u * u;
        }
        let (first, second) = (first * h / 3.0, second * h / 3.0);
        Some(WallInterval { u_tau, y1, law: [nu * first, nu * u_tau * second], trapezoid: [0.5 * u1 * y1, 0.5 * u1 * u1 * y1] })
    }
}

/// When a gated case has settled: after at least `min_flow_throughs` flow-throughs, the time
/// mean of every gated quantity over the last averaging window differs from its mean over the
/// window before by less than `drift` of its magnitude. Each case sets its averaging window to
/// whole periods of its slowest oscillation, so a stationary oscillation does not read as drift.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settling {
    pub min_flow_throughs: f64,
    pub drift: f64,
}

/// The settling rule for the gated benchmark cases.
pub const GATED: Settling = Settling { min_flow_throughs: 2.0, drift: 1e-3 };

impl Settling {
    /// Whether `history` (one row of quantities per sample) has settled after `flow_throughs`,
    /// with averaging windows of `window` samples.
    pub fn settled(&self, history: &[Vec<f64>], window: usize, flow_throughs: f64) -> bool {
        flow_throughs >= self.min_flow_throughs && means_settled(history, window, self.drift)
    }
}

/// The mean of each column of `rows` (`None` for no rows or rows of different lengths).
pub fn column_means(rows: &[Vec<f64>]) -> Option<Vec<f64>> {
    let first = rows.first()?;
    if rows.iter().any(|r| r.len() != first.len()) {
        return None;
    }
    Some((0..first.len()).map(|c| rows.iter().map(|r| r[c]).sum::<f64>() / rows.len() as f64).collect())
}

/// Whether the means of the last `window` rows of `history` differ from the means of the `window`
/// rows before them by less than `relative` of their magnitude in every column. Too short a
/// history, or a non-finite mean, is not settled.
pub fn means_settled(history: &[Vec<f64>], window: usize, relative: f64) -> bool {
    let Some(start) = history.len().checked_sub(2 * window) else {
        return false;
    };
    let (Some(then), Some(now)) = (column_means(&history[start..start + window]), column_means(&history[start + window..])) else {
        return false;
    };
    window > 0 && now.len() == then.len() && now.iter().zip(&then).all(|(a, b)| a.is_finite() && b.is_finite() && (a - b).abs() <= relative * a.abs())
}

/// One quantity of a case against its reference.
#[derive(Debug, Clone)]
pub struct Gate {
    pub case: String,
    pub quantity: String,
    pub value: f64,
    pub reference: f64,
    /// The allowed relative error.
    pub relative: f64,
}

impl Gate {
    pub fn new(case: &str, quantity: &str, value: f64, reference: f64, relative: f64) -> Gate {
        Gate { case: case.into(), quantity: quantity.into(), value, reference, relative }
    }

    /// The relative error (value − reference) / |reference|.
    pub fn error(&self) -> f64 {
        (self.value - self.reference) / self.reference.abs()
    }

    pub fn passes(&self) -> bool {
        self.value.is_finite() && self.error().abs() <= self.relative
    }
}

/// The gates of one or more cases, printed as a table.
#[derive(Debug, Clone, Default)]
pub struct Score {
    pub gates: Vec<Gate>,
}

impl Score {
    pub fn push(&mut self, gate: Gate) {
        self.gates.push(gate);
    }

    pub fn passes(&self) -> bool {
        !self.gates.is_empty() && self.gates.iter().all(Gate::passes)
    }

    /// The table: one row per gate, then the verdict.
    pub fn table(&self) -> String {
        let mut s = String::from("| case | quantity | value | reference | error | tolerance | |\n|---|---|---|---|---|---|---|\n");
        for g in &self.gates {
            s += &format!(
                "| {} | {} | {:.5e} | {:.5e} | {:+.2} % | ±{:.1} % | {} |\n",
                g.case,
                g.quantity,
                g.value,
                g.reference,
                100.0 * g.error(),
                100.0 * g.relative,
                if g.passes() { "pass" } else { "FAIL" }
            );
        }
        let passing = self.gates.iter().filter(|g| g.passes()).count();
        s += &format!("\n{passing} of {} gates pass\n", self.gates.len());
        s
    }
}
