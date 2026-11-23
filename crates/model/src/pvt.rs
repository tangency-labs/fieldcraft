//! Pressure–volume–temperature (PVT) behaviour of a material: its specific volume, and so its
//! density, as a function of temperature and pressure.
//!
//! The model is the modified two-domain Tait equation for polymers, which splits the
//! temperature range at the transition temperature T_t(p) (the glass transition of an amorphous
//! polymer, the crystallisation of a semi-crystalline one) into a solid and a melt domain, each
//! with its own coefficients b1..b4:
//!
//! ```text
//! v(T, p) = v0(T) · [1 − C · ln(1 + p / B(T))] + vt(T, p),   C = 0.0894
//! v0(T)   = b1 + b2 · (T − b5)
//! B(T)    = b3 · exp(−b4 · (T − b5))
//! T_t(p)  = b5 + b6 · p
//! vt      = b7 · exp(b8 · (T − b5) − b9 · p)   in the solid domain (T ≤ T_t(p)), 0 in the melt
//! ```
//!
//! C is the universal constant of the Tait equation for polymers. The two-domain form goes back to
//! R. Y. Chang, C. H. Chen, K. S. Su, *Modifying the Tait equation with cooling-rate effects to
//! predict the pressure-volume-temperature behaviors of amorphous polymers: modeling and
//! experiments*, Polym. Eng. Sci. 36(13):1789–1795, 1996, doi:10.1002/pen.10574; it is written
//! out as above, with the role of each coefficient, in T. Willerer, T. Brinkmann, K. Drechsler,
//! *Development and Application of a Cooling Rate Dependent PVT Model for Injection Molding
//! Simulation of Semi Crystalline Thermoplastics*, Polymers 16(22):3194, 2024,
//! doi:10.3390/polym16223194 (eqs. 1–9, open access). Temperatures are absolute (kelvin) and
//! pressures in pascal, as everything the model stores; the coefficients' units follow (b1 and b7
//! m³/kg, b2 m³/(kg·K), b3 Pa, b4 and b8 1/K, b5 K, b6 K/Pa, b9 1/Pa).

use serde::{Deserialize, Serialize};

use crate::{Error, Result, finite, nonnegative, positive};

/// Tait's universal constant C (dimensionless).
pub const TAIT_C: f64 = 0.0894;

/// One domain's coefficients of the Tait equation: the zero-pressure specific volume
/// v0 = b1 + b2 (T − b5) and the pressure sensitivity B = b3 exp(−b4 (T − b5)).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TaitDomain {
    /// Specific volume at the transition temperature b5 and zero pressure (m³/kg).
    pub b1: f64,
    /// Rise of the zero-pressure specific volume with temperature (m³/(kg·K)).
    pub b2: f64,
    /// Pressure sensitivity B at the transition temperature b5 (Pa).
    pub b3: f64,
    /// Fall of B with temperature (1/K).
    pub b4: f64,
}

impl TaitDomain {
    fn validate(&self, domain: &str) -> Result<()> {
        positive(&format!("{domain} b1"), self.b1)?;
        nonnegative(&format!("{domain} b2"), self.b2)?;
        positive(&format!("{domain} b3"), self.b3)?;
        nonnegative(&format!("{domain} b4"), self.b4)?;
        Ok(())
    }
}

/// A material's PVT behaviour.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Pvt {
    /// The modified two-domain Tait equation (module docs). b7, b8 and b9 describe the extra
    /// volume of a semi-crystalline polymer's solid domain and default to 0 (an amorphous one).
    Tait {
        melt: TaitDomain,
        solid: TaitDomain,
        /// Transition temperature at zero pressure (K).
        b5: f64,
        /// Rise of the transition temperature with pressure (K/Pa).
        b6: f64,
        /// Transition volume term at b5 and zero pressure (m³/kg).
        #[serde(default)]
        b7: f64,
        /// Its rate with temperature (1/K).
        #[serde(default)]
        b8: f64,
        /// Its rate with pressure (1/Pa).
        #[serde(default)]
        b9: f64,
    },
}

/// The domain of the Tait equation a state is in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PvtDomain {
    /// Above the transition temperature.
    Melt,
    /// At or below it.
    Solid,
}

/// A material's specific volume and density at one temperature and pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PvtState {
    /// Specific volume v (m³/kg).
    pub specific_volume: f64,
    /// Density ρ = 1/v (kg/m³).
    pub density: f64,
    /// The transition temperature T_t(p) at this pressure (K).
    pub transition_temperature: f64,
    pub domain: PvtDomain,
}

impl Pvt {
    pub fn validate(&self) -> Result<()> {
        match self {
            Pvt::Tait { melt, solid, b5, b6, b7, b8, b9 } => {
                melt.validate("melt")?;
                solid.validate("solid")?;
                positive("b5 (the transition temperature, K)", *b5)?;
                for (what, v) in [("b6", b6), ("b7", b7), ("b8", b8), ("b9", b9)] {
                    nonnegative(what, *v)?;
                }
            }
        }
        Ok(())
    }

    /// The specific volume and density at absolute temperature `temperature` (K) and pressure
    /// `pressure` (Pa, ≥ 0). Refuses a state outside the equation's range: a non-finite or
    /// non-positive temperature, a negative pressure, or one where the equation gives no
    /// positive, finite volume (a pressure far beyond the fitted range).
    pub fn evaluate(&self, temperature: f64, pressure: f64) -> Result<PvtState> {
        positive("temperature (K)", temperature)?;
        nonnegative("pressure (Pa)", pressure)?;
        match *self {
            Pvt::Tait { melt, solid, b5, b6, b7, b8, b9 } => {
                let transition = finite("transition temperature", b5 + b6 * pressure)?;
                let dt = temperature - b5;
                let (domain, c, vt) = if temperature > transition {
                    (PvtDomain::Melt, melt, 0.0)
                } else {
                    (PvtDomain::Solid, solid, b7 * libm::exp(b8 * dt - b9 * pressure))
                };
                let v0 = c.b1 + c.b2 * dt;
                let b = c.b3 * libm::exp(-c.b4 * dt);
                let v = v0 * (1.0 - TAIT_C * libm::log1p(pressure / b)) + vt;
                if !(v0 > 0.0 && b > 0.0 && b.is_finite() && v > 0.0 && v.is_finite()) {
                    return Err(Error::Invalid(format!(
                        "the Tait equation has no positive, finite specific volume at {temperature} K and {pressure} Pa (outside its range)"
                    )));
                }
                Ok(PvtState { specific_volume: v, density: 1.0 / v, transition_temperature: transition, domain })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A constructed polypropylene-like parameter set (not a published fit), chosen so that
    /// realistic states give the densities of unfilled polypropylene in handbooks: about
    /// 900–910 kg/m³ solid at room temperature, 750–775 kg/m³ melt at 200–230 °C, and a melt
    /// about 10 % denser at 100 MPa. Transition at 400 K (127 °C) rising 25 K per 100 MPa;
    /// b9 > b8 · b6 keeps the solid's transition volume from growing along T_t(p).
    const B5: f64 = 400.0;
    const B6: f64 = 2.5e-7;
    const B7: f64 = 5.0e-5;
    const B8: f64 = 0.1;
    const B9: f64 = 4.0e-8;
    const MELT: [f64; 4] = [1.23e-3, 9.0e-7, 8.0e7, 4.5e-3];
    const SOLID: [f64; 4] = [1.14e-3, 4.0e-7, 1.5e8, 3.5e-3];

    fn pp() -> Pvt {
        let domain = |[b1, b2, b3, b4]: [f64; 4]| TaitDomain { b1, b2, b3, b4 };
        Pvt::Tait { melt: domain(MELT), solid: domain(SOLID), b5: B5, b6: B6, b7: B7, b8: B8, b9: B9 }
    }

    /// The equation written out again, term by term.
    fn tait_by_hand(t: f64, p: f64) -> f64 {
        let melt = t > B5 + B6 * p;
        let [b1, b2, b3, b4] = if melt { MELT } else { SOLID };
        let v0 = b1 + b2 * (t - B5);
        let big_b = b3 * f64::exp(-b4 * (t - B5));
        let vt = if melt { 0.0 } else { B7 * f64::exp(B8 * (t - B5) - B9 * p) };
        v0 * (1.0 - 0.0894 * f64::ln(1.0 + p / big_b)) + vt
    }

    fn v(pvt: &Pvt, t: f64, p: f64) -> f64 {
        pvt.evaluate(t, p).unwrap().specific_volume
    }

    #[test]
    fn polypropylene_like_set_gives_realistic_densities() {
        let pp = pp();
        pp.validate().unwrap();
        // The same equation evaluated separately in Python (math.log1p): a check of this
        // implementation, not of the model against measured data.
        for (t, p, want) in [
            (296.15, 0.0, 1.0984615446205544e-3),
            (296.15, 1e8, 1.0610612827608727e-3),
            (503.15, 0.0, 1.322835e-3),
            (503.15, 1e8, 1.1933706069147847e-3),
        ] {
            let got = v(&pp, t, p);
            assert!((got - want).abs() <= 1e-13 * want, "T {t} K, p {p} Pa: {got} vs {want}");
        }
        // Realistic densities: unfilled PP solid at 23 °C, melt at 200 and 230 °C.
        let rho = |t: f64, p: f64| pp.evaluate(t, p).unwrap().density;
        assert!((900.0..920.0).contains(&rho(296.15, 0.0)), "{}", rho(296.15, 0.0));
        assert!((765.0..780.0).contains(&rho(473.15, 0.0)), "{}", rho(473.15, 0.0));
        assert!((750.0..765.0).contains(&rho(503.15, 0.0)), "{}", rho(503.15, 0.0));
        let hot = pp.evaluate(473.15, 0.0).unwrap();
        assert!((hot.specific_volume - (1.23e-3 + 9.0e-7 * 73.15)).abs() < 1e-18, "{hot:?}");
        assert!((hot.density - 1.0 / hot.specific_volume).abs() < 1e-12);
        assert_eq!(hot.domain, PvtDomain::Melt);
        // At b5 and zero pressure the solid branch is b1s + b7.
        assert!((v(&pp, B5, 0.0) - (1.14e-3 + 5.0e-5)).abs() < 1e-18);
        // Across 0–200 MPa and 20–250 °C the model is the equation.
        for t in [293.15, 323.15, 373.15, 400.0, 420.0, 450.0, 473.15, 523.15] {
            for p in [0.0, 5e6, 5e7, 1e8, 1.5e8, 2e8] {
                let (got, want) = (v(&pp, t, p), tait_by_hand(t, p));
                assert!((got - want).abs() <= 1e-14 * want, "T {t} K, p {p} Pa: {got} vs {want}");
            }
        }
        // The transition temperature rises with pressure: 400 K + 2.5e-7 K/Pa · 100 MPa.
        let s = pp.evaluate(420.0, 1e8).unwrap();
        assert!((s.transition_temperature - 425.0).abs() < 1e-9);
        assert_eq!(s.domain, PvtDomain::Solid, "420 K is below T_t at 100 MPa");
        assert_eq!(pp.evaluate(420.0, 0.0).unwrap().domain, PvtDomain::Melt);
        // The solid is denser than the melt, also just below the transition.
        for p in [0.0, 5e7, 1.5e8] {
            let tt = B5 + B6 * p;
            assert!(pp.evaluate(tt, p).unwrap().density > pp.evaluate(tt + 1e-6, p).unwrap().density, "{p}");
        }
    }

    #[test]
    fn zero_pressure_gives_the_zero_pressure_volume() {
        let pp = pp();
        for t in [300.0, 380.0, 450.0, 520.0] {
            let melt = t > B5;
            let [b1, b2, ..] = if melt { MELT } else { SOLID };
            let vt = if melt { 0.0 } else { B7 * f64::exp(B8 * (t - B5)) };
            assert!((v(&pp, t, 0.0) - (b1 + b2 * (t - B5) + vt)).abs() < 1e-18, "{t}");
        }
    }

    #[test]
    fn volume_falls_with_pressure_and_rises_with_temperature() {
        let pp = pp();
        let pressures: Vec<f64> = (0..=40).map(|i| f64::from(i) * 5e6).collect();
        for t in [300.0, 350.0, 400.0, 430.0, 480.0, 530.0] {
            for w in pressures.windows(2) {
                if let [a, b] = *w {
                    assert!(v(&pp, t, b) < v(&pp, t, a), "T {t} K: {a} → {b} Pa");
                }
            }
        }
        // Within each domain (the volume steps down at the transition on cooling).
        for p in [0.0, 5e7, 1.5e8] {
            let tt = B5 + B6 * p;
            for (lo, hi) in [(280.0, tt), (tt + 1e-6, 560.0)] {
                let steps: Vec<f64> = (0..=50).map(|i| lo + (hi - lo) * f64::from(i) / 50.0).collect();
                for w in steps.windows(2) {
                    if let [a, b] = *w {
                        assert!(v(&pp, b, p) > v(&pp, a, p), "p {p} Pa: {a} → {b} K");
                    }
                }
            }
        }
    }

    #[test]
    fn branches_meet_at_the_transition_when_the_volumes_match() {
        // Zero pressure: the solid branch at b5 is b1s + b7, the melt's limit b1m; with
        // b1m = b1s + b7 they meet whatever b2, b3, b4 each domain has.
        let solid = TaitDomain { b1: 8.4e-4, b2: 3.0e-7, b3: 1.5e8, b4: 3.0e-3 };
        let melt = TaitDomain { b1: 8.4e-4 + 7.0e-5, b2: 7.0e-7, b3: 7.0e7, b4: 4.0e-3 };
        let pvt = Pvt::Tait { melt, solid, b5: 400.0, b6: 2.0e-7, b7: 7.0e-5, b8: 0.2, b9: 5.0e-8 };
        let at = pvt.evaluate(400.0, 0.0).unwrap();
        let above = pvt.evaluate(400.0 + 1e-9, 0.0).unwrap();
        assert_eq!((at.domain, above.domain), (PvtDomain::Solid, PvtDomain::Melt));
        assert!((at.specific_volume - above.specific_volume).abs() < 1e-15, "{at:?} {above:?}");
        // Under pressure, with no transition volume (vt = 0, an amorphous polymer) and a shared
        // B: the branches meet at T_t(p) when v0 does, b1s = b1m + (b2m − b2s)(T_t − b5).
        let p = 1e8;
        let tt = 400.0 + 2.0e-7 * p;
        let melt = TaitDomain { b1: 9.0e-4, b2: 6.0e-7, b3: 1.0e8, b4: 4.0e-3 };
        let solid = TaitDomain { b1: 9.0e-4 + (6.0e-7 - 2.5e-7) * (tt - 400.0), b2: 2.5e-7, ..melt };
        let pvt = Pvt::Tait { melt, solid, b5: 400.0, b6: 2.0e-7, b7: 0.0, b8: 0.0, b9: 0.0 };
        pvt.validate().unwrap();
        let below = pvt.evaluate(tt, p).unwrap();
        let above = pvt.evaluate(tt + 1e-9, p).unwrap();
        assert_eq!((below.domain, above.domain), (PvtDomain::Solid, PvtDomain::Melt));
        assert!((below.specific_volume - above.specific_volume).abs() < 1e-15, "{below:?} {above:?}");
        assert!((below.transition_temperature - tt).abs() < 1e-12);
    }

    #[test]
    fn hostile_coefficients_and_states_are_refused() {
        let ok = pp();
        // Every coefficient of a copy of the set, in a fixed order, to spoil one.
        let with = |i: usize, bad: f64| {
            let mut p = ok;
            let Pvt::Tait { melt, solid, b5, b6, b7, b8, b9 } = &mut p;
            let fields = [
                &mut melt.b1,
                &mut melt.b2,
                &mut melt.b3,
                &mut melt.b4,
                &mut solid.b1,
                &mut solid.b2,
                &mut solid.b3,
                &mut solid.b4,
                b5,
                b6,
                b7,
                b8,
                b9,
            ];
            if let Some(f) = fields.into_iter().nth(i) {
                *f = bad;
            }
            p
        };
        // Non-finite: every coefficient. Negative: every one (b1, b3 and b5 also at zero).
        for i in 0..13 {
            for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e-9] {
                assert!(with(i, bad).validate().is_err(), "coefficient {i} = {bad}");
            }
        }
        for i in [0, 2, 4, 6, 8] {
            assert!(with(i, 0.0).validate().is_err(), "coefficient {i} = 0");
        }
        // States outside the equation's range are errors, never NaN or a panic.
        for (t, p) in [
            (f64::NAN, 0.0),
            (0.0, 0.0),
            (-10.0, 0.0),
            (f64::INFINITY, 0.0),
            (300.0, -1.0),
            (300.0, f64::NAN),
            (300.0, f64::INFINITY),
            (300.0, 1e300),
            (1e300, 0.0),
        ] {
            assert!(ok.evaluate(t, p).is_err(), "T {t}, p {p}");
        }
        // A pressure far beyond the fitted range makes the bracket negative: refused.
        assert!(ok.evaluate(480.0, 1e16).is_err());
        // JSON: b7..b9 default to 0, unknown keys are refused.
        let amorphous: Pvt = serde_json::from_value(serde_json::json!({
            "type": "tait", "b5": 373.0, "b6": 2e-7,
            "melt": {"b1": 9.7e-4, "b2": 5.7e-7, "b3": 1.8e8, "b4": 4.3e-3},
            "solid": {"b1": 9.7e-4, "b2": 2.3e-7, "b3": 2.6e8, "b4": 3.0e-3}
        }))
        .unwrap();
        assert!(matches!(amorphous, Pvt::Tait { b7, b8, b9, .. } if b7 == 0.0 && b8 == 0.0 && b9 == 0.0));
        assert!(serde_json::from_value::<Pvt>(serde_json::json!({"type": "tait", "b5": 1, "b6": 0, "b10": 1})).is_err());
        assert!(serde_json::from_value::<Pvt>(serde_json::json!({"type": "spencer_gilmore"})).is_err());
    }
}
