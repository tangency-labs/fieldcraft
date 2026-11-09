//! Read-only fatigue postprocessing of current stored linear reference results.
use crate::schema::{Field, json as field_json};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};
use fieldcraft_solver::{
    StepResult as ImplicitResult,
    fatigue::{curve::SnCurve, life::StressMean, stress::Criterion},
};
use serde::Deserialize;
use serde_json::{Value, json};

const FORM: &[Field] = &[
    field_json("channels", "Reference channels")
        .req()
        .help("[{step,coefficients:[dimensionless load factors]}]; one value per time, same homogeneous supports, current model."),
    field_json("times_seconds", "Times (s)").req().help("Strictly increasing finite SI seconds; no interpolation."),
    field_json("nodes", "Node labels").req().help("Node labels shared by every reference result, or \"all\" for every reference node."),
    field_json("weights", "Stress projection").help("Fixed [xx,yy,zz,xy,xz,yz] weights; applied after tensor summation."),
    field_json("criterion", "Tensor stress criterion").help("Use instead of weights: type=projection with weights, signed_von_mises (sign of absolute-largest principal; ties refused), or max_principal (algebraic maximum). Tensor criteria require proportional complete global Pa histories."),
    field_json("curve", "S-N curve (SI)").req().help("Basquin: amplitude_ref Pa, cycles_ref, exponent, endurance_limit Pa; or points [{cycles,amplitude Pa}], endurance_limit, surface_factor, size_factor, extrapolation end_slopes/reject."),
    field_json("correction", "Mean stress correction").req().help("Tagged type none/goodman/gerber/soderberg and strength in Pa."),
];
pub const HISTORY: CommandSpec = CommandSpec {
    id: "solve.fatigue.history",
    label: "Reference-load Fatigue",
    params: r#"{channels:[{step,coefficients}],times_seconds,nodes,weights?|criterion?,curve,correction} — sums stored linear static reference stress tensors on the current model, then scalarizes (projection or proportional signed von Mises/algebraic max principal), rainflow counts finite records and sums Miner damage. Coefficients multiply full reference load patterns, not assumed 1 N. Requires identical homogeneous supports; excludes thermal, contact and nonlinear base states. Stores fatigue_damage/fatigue_life contours on the first reference step, with explicit node coverage and ranked nodes; no transient dynamics."#,
    form: Some(FORM),
    kind: CommandKind::Model,
    requires: &[],
    run,
};
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Channel {
    step: String,
    coefficients: Vec<f64>,
}
struct Input {
    channels: Vec<Channel>,
    times_seconds: Vec<f64>,
    nodes: Value,
    criterion: Criterion,
    curve: SnCurve,
    correction: StressMean,
}
pub(super) fn field<T: serde::de::DeserializeOwned>(p: &Value, key: &str) -> Result<T> {
    let value = p.get(key).filter(|v| !v.is_null()).ok_or_else(|| EngineError::param(key, crate::schema::ProblemKind::Required))?;
    serde_json::from_value(value.clone()).map_err(|e| EngineError::invalid(key, e.to_string()))
}
pub(super) fn criterion(p: &Value) -> Result<Criterion> {
    if p.get("criterion").is_some() {
        if p.get("weights").is_some() {
            return Err(EngineError::invalid("criterion", "use criterion or legacy weights, not both"));
        }
        let c: Criterion = field(p, "criterion")?;
        c.validate()?;
        Ok(c)
    } else {
        Ok(Criterion::Projection { weights: field(p, "weights")? })
    }
}
fn run(s: &mut Session, p: &Value) -> Result<Value> {
    let input = Input {
        channels: field(p, "channels")?,
        times_seconds: field(p, "times_seconds")?,
        nodes: field(p, "nodes")?,
        criterion: criterion(p)?,
        curve: field(p, "curve")?,
        correction: field(p, "correction")?,
    };
    if input.channels.is_empty() || input.channels.len() > 64 {
        return Err(EngineError::invalid("channels", "need 1..64 reference channels"));
    }
    let mut boundary = None;
    let mut held = Vec::new();
    let mut coefficients = Vec::new();
    for channel in &input.channels {
        let step = s.model.step(&channel.step).ok_or_else(|| EngineError::invalid("channels", "reference step absent"))?;
        if !matches!(step.procedure, fieldcraft_model::Procedure::Static {}) || !step.output.contains(&fieldcraft_model::OutputVariable::Stress) {
            return Err(EngineError::invalid("channels", "requires linear static steps with stress output"));
        }
        if boundary.as_ref().is_some_and(|b| b != &step.boundary_conditions) {
            return Err(EngineError::invalid("channels", "reference support sets differ"));
        }
        boundary = Some(step.boundary_conditions.clone());
        for name in &step.boundary_conditions {
            let bc = s.model.boundary_condition(name).ok_or_else(|| EngineError::invalid("channels", "missing support"))?;
            if bc.drives() || bc.is_thermal() {
                return Err(EngineError::invalid("channels", "nonhomogeneous or thermal supports would repeat an affine base state"));
            }
        }
        held.push(s.result(&channel.step).ok_or_else(|| {
            EngineError::invalid("channels", "no current stored reference result; solve all reference steps after defining the model")
        })?);
        coefficients.push(channel.coefficients.clone());
    }
    let results = held
        .iter()
        .map(|r| match r.as_ref() {
            StepResult::Implicit(ImplicitResult::Static(r)) => Ok(r),
            _ => Err(EngineError::invalid("channels", "stored reference is not static")),
        })
        .collect::<Result<Vec<_>>>()?;
    let labels: Vec<u32> = if input.nodes == "all" {
        results.first().ok_or_else(|| EngineError::invalid("channels", "no reference"))?.nodes.clone()
    } else {
        serde_json::from_value(input.nodes.clone()).map_err(|e| EngineError::invalid("nodes", e.to_string()))?
    };
    let nodes = fieldcraft_solver::fatigue::history::evaluate(
        &results,
        &coefficients,
        &input.times_seconds,
        &labels,
        input.criterion,
        &input.curve,
        input.correction,
    )?;
    let stored = results.first().ok_or_else(|| EngineError::invalid("channels", "no reference"))?;
    let mut stored = (**stored).clone();
    let map = nodes.iter().map(|n| (n.node, n)).collect::<std::collections::BTreeMap<_, _>>();
    let dependencies = input
        .channels
        .iter()
        .map(|c| {
            fieldcraft_io::results::definition_hash(&s.model, &c.step, &s.sources)
                .map(|h| (c.step.clone(), h))
                .ok_or_else(|| EngineError::invalid("channels", "missing source identity"))
        })
        .collect::<Result<Vec<_>>>()?;
    stored.fatigue = Some(fieldcraft_solver::fatigue::field::FatigueField {
        damage: stored.nodes.iter().map(|n| map.get(n).map(|v| v.damage)).collect(),
        life_seconds: stored.nodes.iter().map(|n| map.get(n).and_then(|v| v.repeated_record_life_seconds)).collect(),
        inputs_json: p.to_string(),
        dependencies,
    });
    if let Some(f) = &stored.fatigue {
        f.validate(stored.nodes.len())?;
    }
    let coverage_all = nodes.len() == stored.nodes.len();
    let stored_step = stored.step.clone();
    s.store_result(StepResult::Implicit(ImplicitResult::Static(stored)));
    let mut ranked = nodes.iter().collect::<Vec<_>>();
    ranked.sort_by(|a, b| b.damage.total_cmp(&a.damage).then(a.node.cmp(&b.node)));
    Ok(json!({"method":"linear_reference_superposition_rainflow_miner","nodes":nodes,"dangerous_selected_nodes":ranked,
        "stored_step":stored_step,"coverage_all":coverage_all,"source_steps":input.channels.iter().map(|c|&c.step).collect::<Vec<_>>(),"criterion":input.criterion,"weights":match input.criterion { Criterion::Projection {weights} => Some(weights), _ => None },
        "units":{"stress":"Pa","time":"s","life":"s","coefficients":"dimensionless"},
        "limitations":"finite record with half-cycle residues; repeated-record life assumes repetition of counted spectrum; no transient inertia or general multiaxial criterion"}))
}
