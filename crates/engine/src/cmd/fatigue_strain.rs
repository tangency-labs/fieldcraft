//! Counted total-strain cycles to CMB life and Miner damage; no plastic correction is inferred.
use super::fatigue::field;
use crate::schema::{Field, json as field_json};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_solver::fatigue::{
    MAX_HISTORY_POINTS,
    life::{StrainLife, StrainMean},
};
use serde::Deserialize;
use serde_json::{Value, json};

const FORM: &[Field] = &[
    field_json("curve", "Strain-life curve (SI)")
        .req()
        .help("{young_pa,strength_pa,ductility,b,c}: CMB coefficients; total strain is dimensionless, stress in Pa. One of strength_pa/ductility may be zero for an exact single-term limit; not both. Mean corrections require positive strength_pa."),
    field_json("cycles", "Counted total-strain cycles")
        .req()
        .help("[{total_strain_amplitude,count,mean_stress_pa?,maximum_stress_pa?}]; stress and strain must describe the same local cycle."),
    field_json("correction", "Strain mean correction").req().help("none, morrow (needs mean_stress_pa), or swt (needs positive maximum_stress_pa)."),
];
pub const STRAIN_CYCLES: CommandSpec = CommandSpec {
    id: "solve.fatigue.strain_cycles",
    label: "Strain-cycle Fatigue",
    params: "{curve,cycles,correction} — Coffin–Manson–Basquin life and Miner damage from explicitly counted local total-strain cycles; Morrow or SWT; exact elastic/plastic single-term limits; no inferred plastic strain or automatic stress/strain pairing",
    form: Some(FORM),
    kind: CommandKind::Model,
    requires: &[],
    run,
};
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Curve {
    young_pa: f64,
    strength_pa: f64,
    ductility: f64,
    b: f64,
    c: f64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Cycle {
    total_strain_amplitude: f64,
    count: f64,
    mean_stress_pa: Option<f64>,
    maximum_stress_pa: Option<f64>,
}
#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Correction {
    None,
    Morrow,
    Swt,
}
fn run(_: &mut Session, p: &Value) -> Result<Value> {
    let curve: Curve = field(p, "curve")?;
    let cycles: Vec<Cycle> = field(p, "cycles")?;
    let correction: Correction = field(p, "correction")?;
    if cycles.is_empty() || cycles.len() > MAX_HISTORY_POINTS {
        return Err(EngineError::invalid("cycles", "need 1..1000000 counted cycles"));
    }
    let law = StrainLife { young: curve.young_pa, strength: curve.strength_pa, ductility: curve.ductility, b: curve.b, c: curve.c };
    let mut damage = 0.0;
    let mut total_count = 0.0;
    let mut rows = Vec::with_capacity(cycles.len());
    for cycle in cycles {
        if !cycle.count.is_finite()
            || cycle.count <= 0.0
            || cycle.mean_stress_pa.is_some_and(|x| !x.is_finite())
            || cycle.maximum_stress_pa.is_some_and(|x| !x.is_finite())
        {
            return Err(EngineError::invalid("cycles", "counts must be finite positive; supplied stresses must be finite"));
        }
        if let (Some(mean), Some(maximum)) = (cycle.mean_stress_pa, cycle.maximum_stress_pa)
            && maximum < mean
        {
            return Err(EngineError::invalid("cycles", "maximum stress cannot be below mean stress for the same cycle"));
        }
        let mean = match correction {
            Correction::None => StrainMean::None,
            Correction::Morrow => StrainMean::Morrow {
                mean_stress: cycle.mean_stress_pa.ok_or_else(|| EngineError::invalid("cycles", "Morrow requires the cycle mean_stress_pa"))?,
            },
            Correction::Swt => StrainMean::Swt {
                maximum_stress: cycle.maximum_stress_pa.ok_or_else(|| EngineError::invalid("cycles", "SWT requires the cycle maximum_stress_pa"))?,
            },
        };
        let life = law.cycles(cycle.total_strain_amplitude, mean)?;
        let increment = life.map_or(0.0, |n| cycle.count / n);
        damage += increment;
        total_count += cycle.count;
        if !damage.is_finite() || !total_count.is_finite() {
            return Err(EngineError::invalid("cycles", "count or damage overflow"));
        }
        rows.push(json!({"total_strain_amplitude":cycle.total_strain_amplitude,"count":cycle.count,"life_cycles":life,"damage":increment,"cumulative_damage":damage}));
    }
    Ok(
        json!({"method":"counted_total_strain_cmb_miner","cycles":rows,"damage":damage,"cycle_count":total_count,"units":{"stress":"Pa","strain":"dimensionless","life":"cycles"},"limitations":"CMB uses reversals 2N internally and reports cycles N; input cycles must pair local total strain with stress; no elastic-to-plastic conversion, no automatic counting, no time-life or field inference; null life means zero amplitude, not a validated endurance limit"}),
    )
}
