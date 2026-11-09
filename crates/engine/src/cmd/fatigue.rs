//! Read-only fatigue postprocessing of current stored linear reference results.
use crate::schema::{Field, json as field_json};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};
use fieldcraft_solver::{
    StepResult as ImplicitResult,
    fatigue::{Basquin, life::StressMean},
};
use serde::Deserialize;
use serde_json::{Value, json};

const FORM: &[Field] = &[
    field_json("channels", "Reference channels")
        .req()
        .help("[{step,coefficients:[dimensionless load factors]}]; one value per time, same homogeneous supports, current model."),
    field_json("times_seconds", "Times (s)").req().help("Strictly increasing finite SI seconds; no interpolation."),
    field_json("nodes", "Node labels").req().help("Selected node labels shared by every reference result."),
    field_json("weights", "Stress projection").req().help("Fixed [xx,yy,zz,xy,xz,yz] weights; applied after tensor summation."),
    field_json("curve", "Basquin curve (SI)").req().help("amplitude_ref Pa, cycles_ref cycles, exponent, endurance_limit Pa."),
    field_json("correction", "Mean stress correction").req().help("Tagged type none/goodman/gerber/soderberg and strength in Pa."),
];
pub const HISTORY: CommandSpec = CommandSpec {
    id: "solve.fatigue.history",
    label: "Reference-load Fatigue",
    params: r#"{channels:[{step,coefficients}],times_seconds,nodes,weights,curve,correction} — sums stored linear static reference stress tensors on the current model, then projects, rainflow counts finite records and sums Miner damage. Coefficients multiply full reference load patterns, not assumed 1 N. Requires identical homogeneous supports; excludes thermal, contact and nonlinear base states. Outputs SI damage/life and ranked selected nodes; no transient dynamics."#,
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
    nodes: Vec<u32>,
    weights: [f64; 6],
    curve: Basquin,
    correction: StressMean,
}
fn field<T: serde::de::DeserializeOwned>(p: &Value, key: &str) -> Result<T> {
    let value = p.get(key).filter(|v| !v.is_null()).ok_or_else(|| EngineError::param(key, crate::schema::ProblemKind::Required))?;
    serde_json::from_value(value.clone()).map_err(|e| EngineError::invalid(key, e.to_string()))
}
fn run(s: &mut Session, p: &Value) -> Result<Value> {
    let input = Input {
        channels: field(p, "channels")?,
        times_seconds: field(p, "times_seconds")?,
        nodes: field(p, "nodes")?,
        weights: field(p, "weights")?,
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
    let nodes = fieldcraft_solver::fatigue::history::evaluate(
        &results,
        &coefficients,
        &input.times_seconds,
        &input.nodes,
        input.weights,
        &input.curve,
        input.correction,
    )?;
    let mut ranked = nodes.iter().collect::<Vec<_>>();
    ranked.sort_by(|a, b| b.damage.total_cmp(&a.damage).then(a.node.cmp(&b.node)));
    Ok(json!({"method":"linear_reference_superposition_rainflow_miner","nodes":nodes,"dangerous_selected_nodes":ranked,
        "source_steps":input.channels.iter().map(|c|&c.step).collect::<Vec<_>>(),"weights":input.weights,
        "units":{"stress":"Pa","time":"s","life":"s","coefficients":"dimensionless"},
        "limitations":"finite record with half-cycle residues; repeated-record life assumes repetition of counted spectrum; no transient inertia or general multiaxial criterion"}))
}
