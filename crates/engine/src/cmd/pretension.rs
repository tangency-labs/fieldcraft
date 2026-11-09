//! Pretension cut creation, loading and locking a current equilibrium adjustment.
use crate::requirement::Requirement;
use crate::schema::{Collection as C, Dim, Field, num, reference, text, vector};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};
use fieldcraft_model::{PretensionControl, PretensionSection, units::Quantity};
use serde_json::{Value, json};

const NAME: Field = text("name", "Section name").req().help("The pretension section name.");
const STEP: Field = reference("step", "Step", C::Step).req().help("Static or nonlinear static target step.");
const CREATE_FORM: &[Field] = &[
    NAME,
    reference("first", "First cut", C::Set).req().help("First conforming cut-node set, with duplicate nodes at the second cut."),
    reference("second", "Second cut", C::Set).req().help("Second cut-node set, geometrically coincident but not sharing nodes."),
    reference("control", "Scalar control point", C::ReferencePoint).req().help("Dedicated reference point; Ux is cut adjustment, not global X."),
    vector("normal", "Cut normal", Dim::None, &["x", "y", "z"]).req().help("Unit direction from first half toward second half."),
];
pub const CREATE: CommandSpec = CommandSpec {
    id: "model.pretension.create",
    label: "New Pretension Section",
    params: r#"{"name":"Bolt","first":"cut1","second":"cut2","normal":[0,0,1],"control":"BoltControl"} — coincident conforming cut nodes; dedicated reference Ux represents uniform axial adjustment"#,
    form: Some(CREATE_FORM),
    kind: CommandKind::Model,
    requires: &[Requirement::Set],
    run: |s, p| {
        let name = super::name(p)?;
        s.model.add_pretension(
            name,
            PretensionSection {
                first: super::req(p, "first")?,
                second: super::req(p, "second")?,
                normal: super::req(p, "normal")?,
                control: super::req(p, "control")?,
            },
        )?;
        Ok(json!({"name":name}))
    },
};
fn set_control(s: &mut Session, p: &Value, force: bool) -> Result<Value> {
    let name = super::name(p)?;
    let step: &str = super::req(p, "step")?;
    let value = super::req(p, "value")?;
    let mut target = s.model.step(step).cloned().ok_or_else(|| EngineError::Params(format!("no step `{step}`")))?;
    target.pretension.insert(name.into(), if force { PretensionControl::Force { value } } else { PretensionControl::Adjustment { value } });
    s.model.replace_step(target)?;
    Ok(json!({"name":name,"step":step}))
}
pub const FORCE: CommandSpec = CommandSpec {
    id: "model.pretension.force",
    label: "Pretension Force",
    params: r#"{"name":"Bolt","step":"Preload","value":1000} — force conjugate to cut adjustment, positive tightens the bolt; display force units"#,
    form: Some(&[NAME, STEP, num("value", "Force", Dim::Of(Quantity::Force)).req().help("Positive tightening force.")]),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: |s, p| set_control(s, p, true),
};
pub const ADJUST: CommandSpec = CommandSpec {
    id: "model.pretension.adjust",
    label: "Pretension Adjustment",
    params: r#"{"name":"Bolt","step":"Service","value":0.001} — prescribed cut adjustment; display length units; not a total bolt length"#,
    form: Some(&[NAME, STEP, num("value", "Adjustment", Dim::Length).req().help("Prescribed relative cut adjustment.")]),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: |s, p| set_control(s, p, false),
};
fn equilibrium<'a>(s: &'a Session, step: &str) -> Result<&'a fieldcraft_solver::StaticResult> {
    match s.result(step).map(|r| r.as_ref()) {
        Some(StepResult::Implicit(fieldcraft_solver::StepResult::Static(r))) => Ok(r),
        _ => Err(EngineError::Params(format!("step `{step}` needs a current static equilibrium result"))),
    }
}
pub const LOCK: CommandSpec = CommandSpec {
    id: "model.pretension.lock",
    label: "Lock Pretension Adjustment",
    params: r#"{"name":"Bolt","step":"Service","from":"Preload"} — snapshot the current solved earlier equilibrium cut adjustment into this step; solve the base again and repeat this command after changing preload"#,
    form: Some(&[NAME, STEP, reference("from", "Solved base", C::Step).req().help("Earlier equilibrium step with a current result.")]),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: |s, p| {
        let name = super::name(p)?;
        let step: &str = super::req(p, "step")?;
        let from: &str = super::req(p, "from")?;
        let source_index = s.model.steps().iter().position(|v| v.name == from);
        let target_index = s.model.steps().iter().position(|v| v.name == step);
        if !matches!((source_index,target_index),(Some(a),Some(b)) if a<b) {
            return Err(EngineError::Params("lock source must precede the target".into()));
        }
        if s.model.materials().values().any(|m| m.plastic.is_some()) {
            return Err(EngineError::Params("locking with plastic history requires a restart state".into()));
        }
        let section = s.model.pretensions().get(name).ok_or_else(|| EngineError::Params(format!("no pretension `{name}`")))?;
        let value = equilibrium(s, from)?
            .reference_points
            .iter()
            .find(|v| v.name == section.control)
            .map(|v| v.displacement[0])
            .ok_or_else(|| EngineError::Params("base result has no pretension control".into()))?;
        let mut target = s.model.step(step).cloned().ok_or_else(|| EngineError::Params(format!("no step `{step}`")))?;
        target.pretension.insert(name.into(), PretensionControl::Adjustment { value });
        s.model.replace_step(target)?;
        Ok(json!({"name":name,"step":step,"from":from,"adjustment":value}))
    },
};
pub const RESULT: CommandSpec = CommandSpec {
    id: "solve.pretension",
    label: "Pretension Result",
    params: r#"{"name":"Bolt","step":"Preload"} — cut adjustment, commanded force and generalized section force; display units"#,
    form: Some(&[NAME, STEP]),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: |s, p| {
        let name = super::name(p)?;
        let step: &str = super::req(p, "step")?;
        let result = equilibrium(s, step)?;
        let cut = result.pretension.iter().find(|v| v.name == name).ok_or_else(|| EngineError::Params("no pretension section result".into()))?;
        Ok(json!({"name":name,"step":step,"adjustment":cut.adjustment,"section_force":cut.section_force,"applied_force":cut.applied_force}))
    },
};

/// Store the resolved adjustment, so replay does not depend on a transient solved result.
pub(crate) fn journal_entry(s: &Session, p: &Value) -> Option<(String, Value)> {
    let name = p.get("name")?.as_str()?;
    let step = p.get("step")?.as_str()?;
    let PretensionControl::Adjustment { value } = s.model.step(step)?.pretension.get(name)? else { return None };
    Some(("model.pretension.adjust".into(), json!({"name":name,"step":step,"value":format!("{value:.17e} m")})))
}
