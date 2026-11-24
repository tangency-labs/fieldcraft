//! Fatigue assessments: how the load of a solved linear static step repeats, for stress-life
//! (S-N) post-processing of its stress field (`fieldcraft_fem::fatigue`, the engine's
//! `results.fatigue`).
//!
//! An assessment belongs to one step, by the step's name: the step's stress field is the response
//! to one reference load, and every stress of the repeated load is that field times a load
//! factor (proportional loading). The loading is a constant amplitude cycle between two factors,
//! a stress ratio, or a history of factors counted by rainflow. Nothing here is read by a solver:
//! an assessment is left out of what a step's result depends on ([`Model::step_inputs`]), so
//! defining or changing one keeps the step's result current.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, Model, Procedure, Result, context, finite, positive};

/// The most load factors a history holds.
pub const MAX_HISTORY: usize = 100_000;

/// The largest load factor (in magnitude) an assessment takes: a cap on hostile input.
pub const MAX_FACTOR: f64 = 1e9; // unit-factors: a cap on a load factor, not a conversion

/// How a history of load factors is counted.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Counting {
    /// The history is one block of a loading that repeats: counted from its largest peak or
    /// deepest valley round to it again, so every cycle closes (Downing and Socie's method for a repeating history).
    #[default]
    Repeated,
    /// The history happens once: rainflow counting by the three-point procedure of ASTM E1049-85, the ranges
    /// left over at the end counted as half cycles.
    Once,
}

/// The load the step's stress field is scaled by, cycle after cycle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Loading {
    /// One cycle repeated: the load goes from `min` to `max` times the step's load and back.
    ConstantAmplitude { max: f64, min: f64 },
    /// One cycle repeated, by its stress ratio R = σ_min / σ_max: from `scale` times the step's
    /// load to R times that (R = −1 fully reversed, R = 0 from zero, 0 < R < 1 fluctuating).
    StressRatio {
        ratio: f64,
        #[serde(default = "one", skip_serializing_if = "is_one")]
        scale: f64,
    },
    /// A sequence of load factors, one block of the loading, cycle-counted by rainflow.
    History {
        factors: Vec<f64>,
        #[serde(default, skip_serializing_if = "is_repeated")]
        counting: Counting,
    },
}

fn one() -> f64 {
    1.0
}

fn is_one(v: &f64) -> bool {
    *v == 1.0
}

fn is_repeated(c: &Counting) -> bool {
    *c == Counting::Repeated
}

impl Loading {
    /// The load factors of a constant amplitude cycle, (max, min); `None` for a history.
    pub fn constant(&self) -> Option<(f64, f64)> {
        match *self {
            Loading::ConstantAmplitude { max, min } => Some((max, min)),
            Loading::StressRatio { ratio, scale } => {
                let (a, b) = (scale, ratio * scale);
                Some((a.max(b), a.min(b)))
            }
            Loading::History { .. } => None,
        }
    }

    /// The smallest and largest load factors (min, max) the loading reaches.
    pub fn extent(&self) -> (f64, f64) {
        match self {
            Loading::History { factors, .. } => factors.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), f| (lo.min(*f), hi.max(*f))),
            _ => self.constant().map_or((0.0, 0.0), |(max, min)| (min, max)),
        }
    }

    /// The largest load factor in magnitude.
    pub fn peak(&self) -> f64 {
        match self {
            Loading::History { factors, .. } => factors.iter().fold(0.0, |m: f64, f| m.max(f.abs())),
            _ => self.constant().map_or(0.0, |(a, b)| a.abs().max(b.abs())),
        }
    }

    pub fn validate(&self) -> Result<()> {
        let factor = |what: &str, v: f64| -> Result<f64> {
            finite(what, v)?;
            if v.abs() > MAX_FACTOR {
                return Err(Error::Invalid(format!("{what} {v} is beyond ±{MAX_FACTOR:e}")));
            }
            Ok(v)
        };
        match self {
            Loading::ConstantAmplitude { max, min } => {
                factor("the largest load factor `max`", *max)?;
                factor("the smallest load factor `min`", *min)?;
                if max < min {
                    return Err(Error::Invalid(format!("`max` ({max}) is below `min` ({min})")));
                }
            }
            Loading::StressRatio { ratio, scale } => {
                factor("the stress ratio R", *ratio)?;
                if ratio.abs() > 1.0 {
                    return Err(Error::Invalid(format!(
                        "the stress ratio R = {ratio} is beyond -1 to 1: R = σ_min/σ_max with |σ_min| ≤ |σ_max|; give a compressive cycle as `constant_amplitude` instead"
                    )));
                }
                positive("the load scale", *scale)?;
                factor("the load scale", *scale)?;
            }
            Loading::History { factors, .. } => {
                if factors.len() < 2 {
                    return Err(Error::Invalid("a load history needs at least two load factors".into()));
                }
                if factors.len() > MAX_HISTORY {
                    return Err(Error::Invalid(format!("{} load factors; a history holds at most {MAX_HISTORY}", factors.len())));
                }
                for (i, f) in factors.iter().enumerate() {
                    factor(&format!("load factor {}", i + 1), *f)?;
                }
            }
        }
        Ok(())
    }
}

/// The scalar stress the S-N curve is read with, from the stress tensor at a node. Both are signed
/// and scale with the load factor, so a proportional load's cycles are those of its factors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StressMeasure {
    /// The von Mises stress with the sign of the principal stress of largest magnitude (ductile
    /// metals).
    SignedVonMises,
    /// The principal stress of largest magnitude, with its sign (brittle materials, and cracks
    /// that grow normal to the largest principal stress).
    MaxAbsPrincipal,
}

impl StressMeasure {
    pub fn id(self) -> &'static str {
        match self {
            StressMeasure::SignedVonMises => "signed_von_mises",
            StressMeasure::MaxAbsPrincipal => "max_abs_principal",
        }
    }
}

/// How the mean stress of a cycle changes the amplitude the S-N curve is read at.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MeanStress {
    /// No correction: the curve is read at the cycle's own amplitude.
    #[serde(rename = "none")]
    Uncorrected,
    /// Goodman's line through the ultimate strength.
    Goodman,
    /// Gerber's parabola through the ultimate strength.
    Gerber,
    /// Soderberg's line through the yield strength.
    Soderberg,
    /// Smith, Watson and Topper's parameter: the curve is read at √(σ_max σ_a), σ_max = σ_a + σ_m
    /// (no strength needed; no damage where σ_max is not tensile).
    Swt,
}

impl MeanStress {
    pub fn id(self) -> &'static str {
        match self {
            MeanStress::Uncorrected => "none",
            MeanStress::Goodman => "goodman",
            MeanStress::Gerber => "gerber",
            MeanStress::Soderberg => "soderberg",
            MeanStress::Swt => "swt",
        }
    }
}

/// A fatigue assessment of one linear static step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FatigueAssessment {
    pub loading: Loading,
    pub measure: StressMeasure,
    pub mean_stress: MeanStress,
    /// Multiplies the whole stress history (every amplitude and every mean alike: the result is
    /// linear in the load) before the S-N curve and the strength checks are read: a stress
    /// concentration the mesh does not resolve, or 1/k for a surface or size factor k of the part
    /// against the specimens of the curve. It is one number for both uses; choose it with that in
    /// mind (a fatigue notch factor is usually applied to the amplitude and, for a notch that
    /// does not yield, to the mean too). 1 (the default) leaves the stresses as solved.
    #[serde(default = "one", skip_serializing_if = "is_one")]
    pub stress_factor: f64,
}

impl FatigueAssessment {
    /// The largest stress factor taken.
    pub const MAX_STRESS_FACTOR: f64 = 100.0;

    pub fn validate(&self) -> Result<()> {
        self.loading.validate()?;
        positive("the stress factor", self.stress_factor)?;
        if self.stress_factor > Self::MAX_STRESS_FACTOR {
            return Err(Error::Invalid(format!("the stress factor {} is above {}", self.stress_factor, Self::MAX_STRESS_FACTOR)));
        }
        Ok(())
    }
}

impl Model {
    /// The fatigue assessments, by the name of the step they assess.
    pub fn fatigue_assessments(&self) -> &BTreeMap<String, FatigueAssessment> {
        &self.fatigue
    }

    /// The fatigue assessment of step `step`.
    pub fn fatigue(&self, step: &str) -> Option<&FatigueAssessment> {
        self.fatigue.get(step)
    }

    /// Give step `step` its fatigue assessment, or replace the one it has. The step must be a
    /// linear static step.
    pub fn set_fatigue(&mut self, step: &str, assessment: FatigueAssessment) -> Result<()> {
        let s = self.step(step).ok_or_else(|| Error::NotFound(format!("step `{step}`")))?;
        if !matches!(s.procedure, Procedure::Static {}) {
            return Err(Error::Invalid(format!(
                "step `{step}` is a {} step: fatigue assesses linear static steps (their stress is proportional to the load)",
                s.procedure.name()
            )));
        }
        assessment.validate().map_err(|e| context("fatigue assessment of step", step, e))?;
        self.fatigue.insert(step.to_string(), assessment);
        Ok(())
    }

    /// Remove step `step`'s fatigue assessment.
    pub fn remove_fatigue(&mut self, step: &str) -> Result<FatigueAssessment> {
        self.fatigue.remove(step).ok_or_else(|| Error::NotFound(format!("fatigue assessment of step `{step}`")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Step;
    use serde_json::json;

    fn assessment(loading: serde_json::Value) -> serde_json::Result<FatigueAssessment> {
        serde_json::from_value(json!({"loading": loading, "measure": "signed_von_mises", "mean_stress": "goodman"}))
    }

    fn model() -> Model {
        let mut m = Model::default();
        m.add_step(Step::new("Load", Procedure::Static {}), None).unwrap();
        m.add_step(Step::new("Modes", Procedure::Frequency { modes: 3, min_hz: None, max_hz: None }), None).unwrap();
        m
    }

    #[test]
    fn an_assessment_round_trips_and_belongs_to_a_static_step() {
        let mut m = model();
        let a = assessment(json!({"type": "stress_ratio", "ratio": -1})).unwrap();
        assert_eq!(a.loading.constant(), Some((1.0, -1.0)));
        assert_eq!(a.stress_factor, 1.0);
        m.set_fatigue("Load", a.clone()).unwrap();
        assert_eq!(m.fatigue("Load"), Some(&a));
        let saved = serde_json::to_value(&m).unwrap();
        assert_eq!(
            saved["fatigue"]["Load"],
            json!({"loading": {"type": "stress_ratio", "ratio": -1.0}, "measure": "signed_von_mises", "mean_stress": "goodman"})
        );
        let back: Model = serde_json::from_value(saved).unwrap();
        assert_eq!(back, m);
        // A model without one writes nothing new.
        assert!(serde_json::to_value(model()).unwrap().get("fatigue").is_none());
        let e = m.set_fatigue("Modes", a.clone()).unwrap_err().to_string();
        assert!(e.contains("frequency step") && e.contains("linear static"), "{e}");
        assert!(m.set_fatigue("Nope", a).is_err());
        // Renaming the step carries its assessment; deleting it takes the assessment along.
        // A procedure that is not linear static is refused while the step has an assessment, and
        // nothing changes; deleting the assessment first allows it.
        let before = m.clone();
        let mut nonlinear = Step::new("Load", Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
        let e = m.replace_step(nonlinear.clone()).unwrap_err().to_string();
        assert!(e.contains("fatigue assessment") && e.contains("model.fatigue.delete") && e.contains("`Load`"), "{e}");
        assert_eq!(m, before);
        nonlinear.procedure = Procedure::Static {};
        m.replace_step(nonlinear).unwrap();
        m.rename(crate::rename::Kind::Step, "Load", "Service").unwrap();
        assert!(m.fatigue("Load").is_none() && m.fatigue("Service").is_some());
        m.remove_step("Service").unwrap();
        assert!(m.fatigue_assessments().is_empty());
    }

    #[test]
    fn a_ratio_gives_its_cycle_and_a_compressive_scale_is_refused() {
        let a = assessment(json!({"type": "stress_ratio", "ratio": 0.1, "scale": 2})).unwrap();
        assert_eq!(a.loading.constant(), Some((2.0, 0.2)));
        a.validate().unwrap();
        assert_eq!(a.loading.peak(), 2.0);
        let h = assessment(json!({"type": "history", "factors": [0, -3, 2]})).unwrap();
        assert_eq!(h.loading.peak(), 3.0);
        assert_eq!(h.loading.constant(), None);
        // Every mean stress choice is read by its id and written back as it.
        for c in [MeanStress::Uncorrected, MeanStress::Goodman, MeanStress::Gerber, MeanStress::Soderberg, MeanStress::Swt] {
            assert_eq!(serde_json::to_value(c).unwrap(), json!(c.id()));
            assert_eq!(serde_json::from_value::<MeanStress>(json!(c.id())).unwrap(), c);
        }
    }

    #[test]
    fn hostile_assessments_are_refused() {
        let mut m = model();
        let before = m.clone();
        for bad in [
            json!({"type": "constant_amplitude", "max": -1, "min": 1}),
            json!({"type": "constant_amplitude", "max": 1e300, "min": 0}),
            json!({"type": "stress_ratio", "ratio": 2}),
            json!({"type": "stress_ratio", "ratio": -1, "scale": 0}),
            json!({"type": "stress_ratio", "ratio": -1, "scale": -2}),
            json!({"type": "history", "factors": [1]}),
            json!({"type": "history", "factors": []}),
            json!({"type": "history", "factors": vec![0.5; MAX_HISTORY + 1]}),
        ] {
            let a = assessment(bad.clone()).unwrap();
            assert!(m.set_fatigue("Load", a).is_err(), "{bad}");
        }
        // Not numbers, unknown keys and unknown choices are refused when read.
        for bad in [
            json!({"type": "constant_amplitude", "max": "one", "min": 0}),
            json!({"type": "constant_amplitude", "max": 1, "mni": 0}),
            json!({"type": "history", "factors": [0, 1], "counting": "twice"}),
            json!({"type": "spectrum"}),
        ] {
            assert!(assessment(bad.clone()).is_err(), "{bad}");
        }
        let mut a = assessment(json!({"type": "stress_ratio", "ratio": 0})).unwrap();
        for f in [0.0, -1.0, f64::NAN, f64::INFINITY, 101.0] {
            a.stress_factor = f;
            assert!(m.set_fatigue("Load", a.clone()).is_err(), "{f}");
        }
        a.loading = Loading::ConstantAmplitude { max: f64::NAN, min: 0.0 };
        a.stress_factor = 1.0;
        assert!(m.set_fatigue("Load", a).is_err());
        assert_eq!(m, before, "refusals change nothing");
        assert!(
            serde_json::from_value::<FatigueAssessment>(
                json!({"loading": {"type": "stress_ratio", "ratio": 0}, "measure": "tresca", "mean_stress": "none"})
            )
            .is_err()
        );
        assert!(
            serde_json::from_value::<FatigueAssessment>(json!({"loading": {"type": "stress_ratio", "ratio": 0}, "measure": "signed_von_mises"}))
                .is_err(),
            "the mean stress choice is required"
        );
    }
}
