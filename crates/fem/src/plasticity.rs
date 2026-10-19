//! Hypoelastic–plastic stress update of an isotropic material: von Mises (J2) plasticity with
//! isotropic, linear kinematic or Johnson–Cook hardening, integrated by radial return, and
//! Johnson–Cook damage.
//!
//! **Radial return** (J. C. Simo, T. J. R. Hughes, *Computational Inelasticity*, Springer 1998,
//! §3.3, Box 3.1 and 3.2; R. D. Krieg, D. B. Krieg, J. Press. Vessel Technol. 99(4):510–515,
//! 1977, for its accuracy): with shear modulus G and bulk modulus K, the pressure takes the
//! volumetric increment elastically, p ← p + K Δθ, and the deviator the elastic trial
//! s* = s + 2G Δe. With the back stress α (kinematic hardening), ξ* = s* − α and the trial
//! von Mises stress q* = √(3/2)‖ξ*‖. If q* exceeds the current yield stress, the equivalent
//! plastic strain increment Δε_p solves the consistency condition
//! q* − 3G Δε_p − H_k Δε_p = σ_y(ε_p + Δε_p, Δε_p / Δt, T), and with the flow direction
//! n = ξ*/‖ξ*‖: s = s* − 2G √(3/2) Δε_p n, α ← α + √(2/3) H_k Δε_p n (Prager's linear kinematic
//! rule; H_k the kinematic modulus, so the uniaxial stress rises by H_k Δε_p either way).
//! - linear hardening: Δε_p = (q* − σ_y)/(3G + H) in closed form;
//! - a hardening table: piecewise linear, solved exactly segment by segment;
//! - Johnson–Cook (rate and temperature dependent): a bracketed scalar Newton iteration with
//!   bisection safeguard on [0, (q* − σ_y)/(3G)] (the root is in it because σ_y only grows with
//!   Δε_p), capped at [`MAX_ITERATIONS`].
//!
//! **Johnson–Cook** flow stress σ_y = (A + B ε_pⁿ)(1 + C ln ε̇*)(1 − T*ᵐ) (G. R. Johnson,
//! W. H. Cook, Proc. 7th Int. Symp. Ballistics, 1983) with ε̇* = ε̇_p / ε̇₀ taken as at least 1
//! (the logarithm stays ≥ 0, which also removes its singularity at rest) and
//! T* = (T − T_room)/(T_melt − T_room) clamped to [0, 1]. The temperature is that at the start
//! of the increment; plastic work heats the material adiabatically, ΔT = χ q Δε_p / (ρ c_p).
//!
//! **Johnson–Cook damage** (G. R. Johnson, W. H. Cook, Eng. Fract. Mech. 21(1):31–48, 1985):
//! D ← D + Δε_p / ε_f with ε_f = (D1 + D2 e^{D3 σ*})(1 + D4 ln ε̇*)(1 + D5 T*), σ* the mean
//! stress over the von Mises stress after the increment. A fracture strain below
//! [`MIN_FRACTURE_STRAIN`] counts as that (some published sets give negative values at high
//! triaxiality, where the point fails at once). The point has failed at D ≥ 1.
//!
//! **One return map for both solvers.** [`J2::radial_return`] is the return above; explicit
//! dynamics calls it through [`J2::update`] on its hypoelastic stress, the implicit material
//! point ([`crate::elastoplastic`]) on its trial elastic strain, small or logarithmic, and adds
//! the consistent tangent [`J2::tangent`].
//!
//! The strain increment the update takes is the **mechanical** one: a caller with initial or
//! thermal strains subtracts their increment first.
//!
//! Stress 6-vectors hold tensor components (xx, yy, zz, xy, xz, yz), strain increments
//! engineering shears.

use fieldcraft_model::{Damage, Elastic, Hardening, Plastic};

use crate::material::von_mises;
use crate::{Error, Result};

/// Most iterations of the rate-dependent return mapping.
pub const MAX_ITERATIONS: u32 = 100;

/// The smallest Johnson–Cook fracture strain.
pub const MIN_FRACTURE_STRAIN: f64 = 1e-6; // unit-factors: a strain

/// The flow stress of a plastic material.
#[derive(Debug, Clone, PartialEq)]
pub enum Flow {
    /// σ_y = y0 + h ε_p (isotropic) or a surface of radius y0 moving with modulus h (kinematic).
    Linear {
        y0: f64,
        h: f64,
        kinematic: bool,
    },
    /// Isotropic hardening by (plastic strain, yield stress) points, constant after the last.
    Table(Vec<(f64, f64)>),
    JohnsonCook {
        a: f64,
        b: f64,
        n: f64,
        c: f64,
        m: f64,
        rate0: f64,
        room: f64,
        melt: f64,
    },
}

/// A hypoelastic, possibly plastic, isotropic material.
#[derive(Debug, Clone, PartialEq)]
pub struct J2 {
    pub shear: f64,
    pub bulk: f64,
    pub flow: Option<Flow>,
    /// Johnson–Cook damage constants D1–D5.
    pub damage: Option<[f64; 5]>,
    /// χ / c_p: the temperature rise per plastic work per mass (0 without heating).
    pub heating: f64,
}

/// The state of a material point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointState {
    pub stress: [f64; 6],
    pub back: [f64; 6],
    pub plastic_strain: f64,
    pub temperature: f64,
    pub damage: f64,
}

impl PointState {
    pub fn new(temperature: f64) -> PointState {
        PointState { stress: [0.0; 6], back: [0.0; 6], plastic_strain: 0.0, temperature, damage: 0.0 }
    }
}

/// What a radial return did: what the consistent tangent needs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Return {
    /// The equivalent plastic strain increment Δε_p.
    pub plastic_strain: f64,
    /// The trial von Mises stress of ξ* = s* − α, q* = √(3/2)‖ξ*‖.
    pub trial: f64,
    /// The flow direction ξ*/‖ξ*‖ (tensor components).
    pub normal: [f64; 6],
    /// The kinematic hardening modulus H_k.
    pub kinematic: f64,
    /// The slope of the flow stress against the plastic strain at the end of the increment.
    pub slope: f64,
}

/// What one update did.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Increment {
    /// The equivalent plastic strain increment.
    pub plastic_strain: f64,
    /// Plastic work per volume, q Δε_p.
    pub plastic_work: f64,
}

/// ‖s‖² of a symmetric tensor's 6-vector (tensor components).
fn norm2(s: &[f64; 6]) -> f64 {
    s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])
}

impl Flow {
    /// The flow stress at plastic strain `e`, plastic strain rate `rate` and temperature `t`, and
    /// its derivative with respect to the plastic strain increment over an increment `dt` (for
    /// the Newton iteration).
    fn stress(&self, e: f64, rate: f64, t: f64, dt: f64) -> (f64, f64) {
        match *self {
            Flow::Linear { y0, h, kinematic } => {
                if kinematic {
                    (y0, 0.0)
                } else {
                    (y0 + h * e, h)
                }
            }
            Flow::Table(ref points) => table(points, e),
            Flow::JohnsonCook { a, b, n, c, m, rate0, room, melt } => {
                let hardening = a + b * e.max(0.0).powf(n);
                let d_hardening = if e > 0.0 {
                    b * n * e.powf(n - 1.0)
                } else if n < 1.0 && b > 0.0 {
                    f64::INFINITY
                } else {
                    b * n
                };
                let ratio = rate / rate0;
                let (rate_term, d_rate) = if ratio > 1.0 && dt > 0.0 { (1.0 + c * ratio.ln(), c / (rate * dt)) } else { (1.0, 0.0) };
                let homologous = ((t - room) / (melt - room)).clamp(0.0, 1.0);
                let thermal = if homologous > 0.0 { 1.0 - homologous.powf(m) } else { 1.0 };
                (hardening * rate_term * thermal, (d_hardening * rate_term + hardening * d_rate) * thermal)
            }
        }
    }

    /// The kinematic hardening modulus.
    fn kinematic(&self) -> f64 {
        match *self {
            Flow::Linear { h, kinematic: true, .. } => h,
            _ => 0.0,
        }
    }

    /// The rate sensitivity's reference rate and the homologous temperature (for damage).
    fn rate_and_temperature(&self, t: f64) -> (f64, f64) {
        match *self {
            Flow::JohnsonCook { rate0, room, melt, .. } => (rate0, ((t - room) / (melt - room)).clamp(0.0, 1.0)),
            _ => (1.0, 0.0),
        }
    }
}

/// The yield stress of a hardening table at plastic strain `e`, and its slope.
fn table(points: &[(f64, f64)], e: f64) -> (f64, f64) {
    for w in points.windows(2) {
        if let [(e0, s0), (e1, s1)] = *w
            && e < e1
        {
            let slope = (s1 - s0) / (e1 - e0);
            return (s0 + slope * (e - e0).max(0.0), slope);
        }
    }
    (points.last().map_or(0.0, |p| p.1), 0.0)
}

impl J2 {
    /// The material of a model's isotropic elastic constants with its plasticity and damage.
    pub fn from_model(elastic: &Elastic, plastic: Option<&Plastic>, damage: Option<&Damage>) -> Result<J2> {
        let Elastic::Isotropic { young, poisson } = *elastic else {
            return Err(Error::NotSupported("plasticity of an orthotropic material".into()));
        };
        let shear = young / (2.0 * (1.0 + poisson));
        let bulk = young / (3.0 * (1.0 - 2.0 * poisson));
        let mut heating = 0.0;
        let flow = plastic.map(|p| match p {
            Plastic::Bilinear { yield_stress, hardening_modulus, hardening } => {
                Flow::Linear { y0: *yield_stress, h: *hardening_modulus, kinematic: *hardening == Hardening::Kinematic }
            }
            Plastic::Table { points } => Flow::Table(points.iter().map(|p| (p.strain, p.stress)).collect()),
            Plastic::JohnsonCook { a, b, n, c, m, reference_strain_rate, room_temperature, melt_temperature, specific_heat, taylor_quinney } => {
                heating = specific_heat.map_or(0.0, |cp| taylor_quinney / cp);
                Flow::JohnsonCook {
                    a: *a,
                    b: *b,
                    n: *n,
                    c: *c,
                    m: *m,
                    rate0: *reference_strain_rate,
                    room: *room_temperature,
                    melt: *melt_temperature,
                }
            }
        });
        let damage = damage.map(|Damage::JohnsonCook { d1, d2, d3, d4, d5 }| [*d1, *d2, *d3, *d4, *d5]);
        if !(shear > 0.0 && bulk > 0.0 && shear.is_finite() && bulk.is_finite()) {
            return Err(Error::Invalid(format!("E = {young}, ν = {poisson} are not a stable material")));
        }
        Ok(J2 { shear, bulk, flow, damage, heating })
    }

    /// Whether the flow stress depends on the temperature (Johnson–Cook softening).
    pub fn is_temperature_dependent(&self) -> bool {
        matches!(self.flow, Some(Flow::JohnsonCook { .. }))
    }

    /// The temperature a point starts at: Johnson–Cook's room temperature; 0 for a law that
    /// does not read it (a placeholder, never reported: see `is_temperature_dependent`).
    pub fn initial_temperature(&self) -> f64 {
        match self.flow {
            Some(Flow::JohnsonCook { room, .. }) => room,
            _ => 0.0,
        }
    }

    /// Update `state` over an increment of length `dt` with the deviatoric strain increment
    /// `dev` (engineering shears; its trace is ignored) and the volumetric strain increment
    /// `vol`, at the current density `density` (for heating).
    pub fn update(&self, state: &mut PointState, dev: &[f64; 6], vol: f64, dt: f64, density: f64) -> Increment {
        let s = &mut state.stress;
        let mean = (s[0] + s[1] + s[2]) / 3.0;
        let tr = (dev[0] + dev[1] + dev[2]) / 3.0;
        let g2 = 2.0 * self.shear;
        // Trial deviator.
        let mut trial = [s[0] - mean + g2 * (dev[0] - tr), s[1] - mean + g2 * (dev[1] - tr), s[2] - mean + g2 * (dev[2] - tr), 0.0, 0.0, 0.0];
        for k in 3..6 {
            trial[k] = s[k] + self.shear * dev[k];
        }
        let pressure = mean + self.bulk * vol;
        let mut out = Increment::default();
        if let Some(r) = self.radial_return(state, &mut trial, dt) {
            let q_new = r.trial - (3.0 * self.shear + r.kinematic) * r.plastic_strain;
            out = Increment { plastic_strain: r.plastic_strain, plastic_work: q_new * r.plastic_strain };
            if self.heating > 0.0 && density > 0.0 {
                state.temperature += self.heating * out.plastic_work / density;
            }
        }
        for k in 0..3 {
            trial[k] += pressure;
        }
        state.stress = trial;
        if let (Some(d), Some(flow)) = (self.damage, &self.flow)
            && out.plastic_strain > 0.0
        {
            let q = von_mises(&state.stress);
            let triaxiality = if q > 0.0 { pressure / q } else { 0.0 };
            let (rate0, homologous) = flow.rate_and_temperature(state.temperature);
            let rate = if dt > 0.0 { out.plastic_strain / dt / rate0 } else { 0.0 };
            let log_rate = if rate > 1.0 { rate.ln() } else { 0.0 };
            let fracture = (d[0] + d[1] * (d[2] * triaxiality).exp()) * (1.0 + d[3] * log_rate) * (1.0 + d[4] * homologous);
            state.damage += out.plastic_strain / fracture.max(MIN_FRACTURE_STRAIN);
        }
        out
    }

    /// The radial return of a trial deviator `trial` (tensor components) from `state`: when the
    /// trial von Mises stress exceeds the yield stress, `trial` is moved back to the yield
    /// surface, the back stress and the equivalent plastic strain of `state` are advanced, and
    /// what the consistent tangent needs is returned; `None` when the point stays elastic.
    pub fn radial_return(&self, state: &mut PointState, trial: &mut [f64; 6], dt: f64) -> Option<Return> {
        let flow = self.flow.as_ref()?;
        let xi: [f64; 6] = std::array::from_fn(|k| trial[k] - state.back[k]);
        let norm = norm2(&xi).sqrt();
        let q = (1.5f64).sqrt() * norm;
        let (yield_now, _) = flow.stress(state.plastic_strain, 0.0, state.temperature, dt);
        if !(q > yield_now && norm > 0.0) {
            return None;
        }
        let hk = flow.kinematic();
        let de = self.plastic_increment(flow, state, q, yield_now, hk, dt);
        let n: [f64; 6] = xi.map(|x| x / norm);
        let radial = 2.0 * self.shear * (1.5f64).sqrt() * de;
        let shift = (2.0f64 / 3.0).sqrt() * hk * de;
        for k in 0..6 {
            trial[k] -= radial * n[k];
            state.back[k] += shift * n[k];
        }
        let rate = if dt > 0.0 { de / dt } else { 0.0 };
        let (_, slope) = flow.stress(state.plastic_strain + de, rate, state.temperature, dt);
        state.plastic_strain += de;
        Some(Return { plastic_strain: de, trial: q, normal: n, kinematic: hk, slope })
    }

    /// The consistent (algorithmic) tangent of the return map, dσ/dε of the stress after the
    /// update against the total strain of the increment (engineering-shear columns): the
    /// elasticity K I ⊗ I + 2G ℙ while elastic, and after a return
    /// K I ⊗ I + 2G θ ℙ − 2G θ̄ n ⊗ n with θ = 1 − 3G Δε_p / q*, θ̄ = 1/(1 + H/3G) − (1 − θ) and
    /// H the slope of the flow stress plus the kinematic modulus at the end of the increment
    /// (Simo & Hughes 1998, Box 3.2, eq. 3.3.10, in the equivalent plastic strain). `bulk`
    /// leaves out the volumetric part K I ⊗ I when false (the isochoric part a mixed element
    /// takes). `None` when 3G + H is not positive: the point softens faster than it unloads
    /// elastically and has no tangent.
    pub fn tangent(&self, r: Option<&Return>, bulk: bool) -> Option<[[f64; 6]; 6]> {
        let g = self.shear;
        let (theta, theta_bar, n) = match r {
            None => (1.0, 0.0, [0.0; 6]),
            Some(r) => {
                let stiff = 3.0 * g + r.slope + r.kinematic;
                if !(stiff > 0.0 && stiff.is_finite() && r.trial > 0.0) {
                    return None;
                }
                let theta = 1.0 - 3.0 * g * r.plastic_strain / r.trial;
                (theta, 3.0 * g / stiff - (1.0 - theta), r.normal)
            }
        };
        let k = if bulk { self.bulk } else { 0.0 };
        let identity = crate::constitutive::identity4();
        let mut d = [[0.0; 6]; 6];
        for a in 0..6 {
            for b in 0..6 {
                let (ia, ib) = (if a < 3 { 1.0 } else { 0.0 }, if b < 3 { 1.0 } else { 0.0 });
                d[a][b] = k * ia * ib + 2.0 * g * theta * (identity[a][b] - ia * ib / 3.0) - 2.0 * g * theta_bar * n[a] * n[b];
            }
        }
        Some(d)
    }

    /// The plastic strain increment of a trial state that yields.
    fn plastic_increment(&self, flow: &Flow, state: &PointState, q: f64, yield_now: f64, hk: f64, dt: f64) -> f64 {
        let g3 = 3.0 * self.shear;
        let (e0, t) = (state.plastic_strain, state.temperature);
        match flow {
            Flow::Linear { h, kinematic, .. } => (q - yield_now) / (g3 + if *kinematic { hk } else { *h }),
            Flow::Table(points) => {
                // Walk the segments from the current one: on a segment of slope h starting at
                // (e_a, s_a), q − 3G Δ = s_a + h (e0 + Δ − e_a).
                for w in points.windows(2) {
                    let [(ea, sa), (eb, sb)] = *w else { continue };
                    if eb <= e0 {
                        continue;
                    }
                    let h = (sb - sa) / (eb - ea);
                    let at = sa + h * (e0 - ea);
                    let de = (q - at) / (g3 + h);
                    if e0 + de <= eb {
                        return de.max(0.0);
                    }
                }
                // Beyond the table: perfectly plastic at the last stress.
                let last = points.last().map_or(yield_now, |p| p.1);
                ((q - last) / g3).max(0.0)
            }
            Flow::JohnsonCook { .. } => {
                let residual = |de: f64| {
                    let rate = if dt > 0.0 { de / dt } else { 0.0 };
                    let (y, dy) = flow.stress(e0 + de, rate, t, dt);
                    (q - g3 * de - y, -g3 - dy)
                };
                let (mut lo, mut hi) = (0.0, ((q - yield_now) / g3).max(0.0));
                let mut de = 0.5 * hi;
                for _ in 0..MAX_ITERATIONS {
                    let (r, dr) = residual(de);
                    if r > 0.0 {
                        lo = de;
                    } else {
                        hi = de;
                    }
                    if r.abs() <= 1e-13 * q || hi - lo <= 1e-15 * hi {
                        break;
                    }
                    let newton = de - r / dr;
                    de = if newton.is_finite() && newton > lo && newton < hi { newton } else { 0.5 * (lo + hi) };
                }
                de
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::HardeningPoint;

    const E: f64 = 200_000.0;
    const NU: f64 = 0.3;

    fn material(plastic: Plastic) -> J2 {
        J2::from_model(&Elastic::Isotropic { young: E, poisson: NU }, Some(&plastic), None).unwrap()
    }

    /// Uniaxial stress at a constant strain rate: the axial strain increments by `de` while the
    /// lateral strain increment is found by a secant iteration on the point update itself so
    /// the lateral stresses vanish. Returns (axial strain, axial stress, plastic strain) per step.
    fn uniaxial(m: &J2, steps: &[f64], dt: f64) -> Vec<(f64, f64, f64)> {
        let mut s = PointState::new(m.initial_temperature());
        let mut out = Vec::new();
        let mut strain = 0.0;
        for de in steps {
            let run = |lateral: f64| {
                let mut t = s;
                let d = [*de, lateral * de, lateral * de, 0.0, 0.0, 0.0];
                m.update(&mut t, &d, d[0] + d[1] + d[2], dt, 7.8e-9);
                t
            };
            let (mut a, mut b) = (-NU, -0.5);
            let (mut fa, mut fb) = (run(a).stress[1], run(b).stress[1]);
            for _ in 0..100 {
                if fb == fa || fb.abs() < 1e-12 {
                    break;
                }
                let c = b - fb * (b - a) / (fb - fa);
                (a, fa) = (b, fb);
                b = c;
                fb = run(b).stress[1];
            }
            s = run(b);
            strain += de;
            out.push((strain, s.stress[0], s.plastic_strain));
        }
        out
    }

    #[test]
    fn bilinear_uniaxial_curve_and_unloading() {
        let (y0, h) = (250.0, 2000.0);
        for kinematic in [false, true] {
            let hardening = if kinematic { Hardening::Kinematic } else { Hardening::Isotropic };
            let m = material(Plastic::Bilinear { yield_stress: y0, hardening_modulus: h, hardening });
            // Load to 1% strain, unload to −1%.
            let mut steps = vec![1e-4; 100];
            steps.extend(vec![-1e-4; 200]);
            let curve = uniaxial(&m, &steps, 1e-3);
            // Loading: σ = E ε below yield, then σ = y0 + H ε_p with ε = σ/E + ε_p, that is
            // σ = y0 + E_t (ε − y0/E), E_t = E H/(E + H).
            let et = E * h / (E + h);
            let (e1, s1, p1) = curve[99];
            let exact = y0 + et * (e1 - y0 / E);
            assert!((s1 - exact).abs() < 1e-9 * exact, "{kinematic}: {s1} vs {exact}");
            assert!((p1 - (s1 - y0) / h).abs() < 1e-12, "{p1}");
            // Unloading is elastic until the reversed yield: isotropic at −σ1, kinematic at
            // σ1 − 2 y0 (the surface keeps its size and moves).
            let reverse = if kinematic { s1 - 2.0 * y0 } else { -s1 };
            for (e, s, _) in &curve[100..] {
                let elastic = s1 - E * (e1 - e);
                if elastic > reverse + 1.0 {
                    assert!((s - elastic).abs() < 1e-8 * s1, "{kinematic}: {s} vs elastic {elastic}");
                }
            }
            let (_, s_end, _) = curve[299];
            assert!(s_end < reverse, "{kinematic}: {s_end} vs {reverse}");
            eprintln!("bilinear kinematic = {kinematic}: σ(1%) = {s1:.6} (exact {exact:.6}); reversed yield at {reverse:.3}; σ(−1%) = {s_end:.4}");
        }
    }

    #[test]
    fn table_and_johnson_cook_follow_their_flow_curves() {
        let points = vec![
            HardeningPoint { strain: 0.0, stress: 200.0 },
            HardeningPoint { strain: 0.01, stress: 300.0 },
            HardeningPoint { strain: 0.05, stress: 340.0 },
        ];
        let m = material(Plastic::Table { points });
        let curve = uniaxial(&m, &vec![2e-4; 400], 1e-3);
        for (_, s, p) in &curve {
            if *p > 0.0 {
                let want = table(&[(0.0, 200.0), (0.01, 300.0), (0.05, 340.0)], *p).0;
                assert!((s - want).abs() < 1e-8 * want, "{s} vs {want} at {p}");
            }
        }
        assert!(curve.last().unwrap().2 > 0.05, "past the table's end");
        // Johnson–Cook at a constant plastic rate: the stress is the flow stress at that rate.
        let (a, b, n, c) = (350.0, 275.0, 0.36, 0.022);
        let jc = Plastic::JohnsonCook {
            a,
            b,
            n,
            c,
            m: 1.0,
            reference_strain_rate: 1.0,
            room_temperature: 293.0,
            melt_temperature: 1356.0,
            specific_heat: None,
            taylor_quinney: 0.9,
        };
        let m = material(jc);
        let dt = 1e-6;
        let curve = uniaxial(&m, &vec![1e-3; 200], dt);
        let mut worst: f64 = 0.0;
        for w in curve.windows(2).skip(20) {
            let (_, s, p) = w[1];
            let rate = (p - w[0].2) / dt;
            let want = (a + b * p.powf(n)) * (1.0 + c * rate.ln());
            worst = worst.max((s - want).abs() / want);
        }
        eprintln!("Johnson–Cook uniaxial: largest relative error of the flow stress {worst:.2e}");
        assert!(worst < 1e-9, "{worst}");
    }

    #[test]
    fn heating_softens_and_damage_accumulates() {
        let jc = |cp: Option<f64>| Plastic::JohnsonCook {
            a: 350.0,
            b: 275.0,
            n: 0.36,
            c: 0.0,
            m: 1.0,
            reference_strain_rate: 1.0,
            room_temperature: 293.0,
            melt_temperature: 1356.0,
            specific_heat: cp,
            taylor_quinney: 0.9,
        };
        let el = Elastic::Isotropic { young: E, poisson: NU };
        let damage = Damage::JohnsonCook { d1: 0.5, d2: 0.0, d3: 0.0, d4: 0.0, d5: 0.0 };
        let cold = J2::from_model(&el, Some(&jc(None)), Some(&damage)).unwrap();
        let hot = J2::from_model(&el, Some(&jc(Some(3.83e8))), Some(&damage)).unwrap();
        let a = uniaxial(&cold, &vec![1e-3; 300], 1e-3);
        let b = uniaxial(&hot, &vec![1e-3; 300], 1e-3);
        assert!(b.last().unwrap().1 < a.last().unwrap().1);
        // Constant fracture strain 0.5: D = ε_p / 0.5.
        let mut s = PointState::new(293.0);
        for _ in 0..300 {
            cold.update(&mut s, &[1e-3, -5e-4, -5e-4, 0.0, 0.0, 0.0], 0.0, 1e-3, 7.8e-9);
        }
        assert!((s.damage - s.plastic_strain / 0.5).abs() < 1e-12, "{} vs {}", s.damage, s.plastic_strain);
    }

    #[test]
    fn hostile_increments_stay_finite() {
        let m = material(Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 0.0, hardening: Hardening::Isotropic });
        let mut s = PointState::new(0.0);
        m.update(&mut s, &[1e3, -1e3, 0.0, 5e2, 0.0, 0.0], 0.0, 0.0, 0.0);
        assert!(s.stress.iter().all(|x| x.is_finite()) && (von_mises(&s.stress) - 250.0).abs() < 1e-9 * 250.0);
    }
}
