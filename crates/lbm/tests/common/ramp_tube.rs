//! The ideal tube behind `inlet_ramp`: one-dimensional isothermal flow (p = c_s² ρ, the lattice's equation of state)
//! between a velocity face, u(0, t) = U r(t) with the raised cosine r of [`fieldcraft_lbm`]'s inlet ramp, and a
//! pressure face, ρ(L, t) = 1, started at rest with ρ = 1. Solved by the method of characteristics: the Riemann
//! invariants J± = u ± c_s ln ρ are constant along dx/dt = u ± c_s (B. Riemann, Abh. Kön. Ges. Wiss. Göttingen 8,
//! 1860; G. B. Whitham, *Linear and Nonlinear Waves*, Wiley, 1974, §6.8). It keeps what linear acoustics drops: the
//! waves travel at c_s ± u, so the round trip between the faces changes while the inlet rises, and the period at the
//! full speed is T₀ / (1 − Ma²).
//!
//! Discretisation: `intervals + 1` nodes at spacing h = L / intervals, a time step Δt = h / (c_s (1 + 2 Ma)) so that
//! the foot of every characteristic stays inside the neighbouring interval while |u| ≤ 2 U (twice the inlet's speed,
//! what one reflection can add to it). The foot of a characteristic is found from its trip over one step: the mean of
//! the speed at the foot (old time) and at the arrival node (new time). The arrival speed at the new time is predicted
//! by one pass of the same scheme and corrected by the second pass (predictor–corrector, so the trip is second order in
//! Δt). J± at the foot is a cubic Lagrange interpolation of the last step's values. At the velocity face J+ =
//! 2 U r(t) − J−; at the pressure face J− = J+ (ln ρ = 0). The solution must stay smooth: a start without a ramp
//! steepens into shocks, which this scheme does not capture.
//!
//! Lengths are in cells and times in lattice steps (c_s = 1/√3 cells per step).

/// The bound the reference's own error may not exceed, relative. The heavy convergence test asserts the Richardson
/// estimate of the reference's error is within it; the fast tests take it as the reference's error term.
pub const REFERENCE_ERROR: f64 = 1e-4;

/// Passes of the fixed-point iteration for a characteristic's foot. The map contracts by Δt |∂u/∂x| per pass, below
/// Ma² here, so two passes leave an error below Ma⁴ of a cell.
const FOOT_PASSES: usize = 2;

/// The tube and its ramp.
#[derive(Debug, Clone, Copy)]
pub struct RampTube {
    /// The distance between the velocity face and the pressure face, in cells.
    pub length: f64,
    /// The inlet's full speed U, in cells per step.
    pub speed: f64,
    /// The ramp's length in steps; zero or less is refused by [`RampTube::density`] (no ramp makes shocks).
    pub ramp_steps: f64,
    /// The probe's distance from the velocity face, in cells.
    pub probe: f64,
    /// Intervals of the grid.
    pub intervals: usize,
}

/// The cubic Lagrange interpolation of `a` at the fractional index `s`, its four-node stencil kept inside the array.
fn cubic(a: &[f64], s: f64) -> f64 {
    let last = a.len() as isize - 3;
    let i = (s.floor() as isize).clamp(1, last) as usize;
    let d = s - i as f64;
    a[i - 1] * (-d * (d - 1.0) * (d - 2.0) / 6.0)
        + a[i] * ((d + 1.0) * (d - 1.0) * (d - 2.0) / 2.0)
        + a[i + 1] * (-(d + 1.0) * d * (d - 2.0) / 2.0)
        + a[i + 2] * ((d + 1.0) * d * (d - 1.0) / 6.0)
}

/// One step of the characteristics from the last step's `plus` and `minus` into `next_plus` and `next_minus`. `arrival[j]`
/// is the speed the trip of node j's characteristics is averaged with at the arrival (the predictor's, or the last
/// step's for the predictor itself).
fn advance(plus: &[f64], minus: &[f64], arrival: &[f64], inlet: f64, courant: f64, per_h: f64, next_plus: &mut [f64], next_minus: &mut [f64]) {
    let n = plus.len() - 1;
    for j in 0..=n {
        let u_arrival = arrival[j];
        // dir = −1: the + characteristic, which comes from the left; dir = +1: the − characteristic.
        let foot = |dir: f64| {
            let mut s = j as f64 + dir * courant - u_arrival * per_h;
            for _ in 0..FOOT_PASSES {
                let u_foot = 0.5 * (cubic(plus, s) + cubic(minus, s));
                s = j as f64 + dir * courant - 0.5 * (u_foot + u_arrival) * per_h;
            }
            s
        };
        let from_left = if j > 0 { cubic(plus, foot(-1.0)) } else { 0.0 };
        let from_right = if j < n { cubic(minus, foot(1.0)) } else { 0.0 };
        next_plus[j] = if j == 0 { 2.0 * inlet - from_right } else { from_left };
        next_minus[j] = if j == n { from_left } else { from_right };
    }
}

impl RampTube {
    /// The inlet's share of its full speed at time `t` (steps): the raised cosine of the solver's ramp.
    fn share(&self, t: f64) -> f64 {
        if t >= self.ramp_steps { 1.0 } else { 0.5 * (1.0 - (std::f64::consts::PI * t / self.ramp_steps).cos()) }
    }

    /// The density at the probe after each of `steps` lattice steps (the solver's own steps are longer; the value at
    /// a lattice step is linear between the two solver steps around it).
    pub fn density(&self, steps: u64) -> Vec<f64> {
        assert!(self.ramp_steps > 0.0 && self.intervals >= 4, "the ideal tube needs a ramp and at least four intervals");
        let cs = 1.0 / 3f64.sqrt();
        let n = self.intervals;
        let h = self.length / n as f64;
        let dt = h / (cs * (1.0 + 2.0 * self.speed / cs));
        let (courant, per_h) = (cs * dt / h, dt / h);
        let at_probe = |plus: &[f64], minus: &[f64]| ((cubic(plus, self.probe / h) - cubic(minus, self.probe / h)) / (2.0 * cs)).exp();
        let (mut plus, mut minus) = (vec![0.0; n + 1], vec![0.0; n + 1]);
        let (mut next_plus, mut next_minus) = (vec![0.0; n + 1], vec![0.0; n + 1]);
        let mut arrival = vec![0.0; n + 1];
        let mut out = Vec::with_capacity(steps as usize);
        let (mut t, mut last) = (0.0, 1.0);
        while (out.len() as u64) < steps {
            let inlet = self.speed * self.share(t + dt);
            for j in 0..=n {
                arrival[j] = 0.5 * (plus[j] + minus[j]);
            }
            // The predictor: the arrival speeds of the new step, from the old ones.
            advance(&plus, &minus, &arrival, inlet, courant, per_h, &mut next_plus, &mut next_minus);
            for j in 0..=n {
                arrival[j] = 0.5 * (next_plus[j] + next_minus[j]);
            }
            // The corrector: the trips averaged with the predicted arrival speeds.
            advance(&plus, &minus, &arrival, inlet, courant, per_h, &mut next_plus, &mut next_minus);
            std::mem::swap(&mut plus, &mut next_plus);
            std::mem::swap(&mut minus, &mut next_minus);
            let now = at_probe(&plus, &minus);
            while (out.len() as u64) < steps && (out.len() + 1) as f64 <= t + dt {
                let w = ((out.len() + 1) as f64 - t) / dt;
                out.push(last + w * (now - last));
            }
            t += dt;
            last = now;
        }
        out
    }
}

/// The rms of `v` about its mean.
pub fn rms(v: &[f64]) -> f64 {
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
}

/// The observed order p of three results on grids of spacing h, h/2 and h/4 (c_h − c_{h/2} against c_{h/2} − c_{h/4}),
/// and the Richardson bound on the error of the coarsest result: c_h − c* ≤ |c_h − c_{h/2}| · 2^p / (2^p − 1).
pub fn richardson(coarse: f64, middle: f64, fine: f64) -> (f64, f64) {
    let (d_coarse, d_fine) = ((coarse - middle).abs(), (middle - fine).abs());
    let order = (d_coarse / d_fine).log2();
    let two_p = 2f64.powf(order);
    (order, d_coarse * two_p / (two_p - 1.0))
}

/// The period of the closed–open tube at rest, T₀ = 4 L / c_s, in steps.
pub fn rest_period(length: f64) -> f64 {
    4.0 * length * 3f64.sqrt()
}

/// Linear acoustics (Ma → 0) for the same tube: the rms of the pressure coefficient at `probe` after a ramp of
/// `ramp_steps`, over whole periods. The step's square wave has the odd harmonics 4 / (π n) of 2 / Ma; the ramp of
/// m = T / T₀ periods keeps |cos(π n m)| / |4 n² m² − 1| of harmonic n (the Fourier transform of the rise's slope),
/// and the mode's shape at the probe is cos(n k x), k = π / (2 L). `harmonics` = 1 is the manual's formula.
pub fn linear_rms(length: f64, speed: f64, ramp_steps: f64, probe: f64, harmonics: usize) -> f64 {
    let pi = std::f64::consts::PI;
    let mach = 3f64.sqrt() * speed;
    let m = ramp_steps / rest_period(length);
    let k = pi / (2.0 * length);
    let mut sum = 0.0;
    for i in 0..harmonics {
        let n = (2 * i + 1) as f64;
        let amplitude = 4.0 / (pi * n) * 2.0 / mach * (pi * n * m).cos().abs() / (4.0 * n * n * m * m - 1.0).abs() * (n * k * probe).cos();
        sum += 0.5 * amplitude * amplitude;
    }
    sum.sqrt()
}
