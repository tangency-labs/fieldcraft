//! Surface foundations: all numeric parameters pass through the unit boundary.
use crate::requirement::Requirement;
use crate::schema::{Collection, Dim, Field, num, reference, text};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_model::{Foundation, units::Quantity};
use serde_json::{Value, json};
const FORM: &[Field] = &[
    text("name", "Name").req().help("The unique foundation definition name."),
    reference("set", "Face set", Collection::Set).req().help("An existing set of element faces."),
    num("normal", "Normal stiffness per area", Dim::Of(Quantity::ForcePerVolume)).req().help("Normal pressure per unit normal displacement."),
    num("tangential", "Tangential stiffness per area", Dim::Of(Quantity::ForcePerVolume))
        .help("Equal stiffness in both in-plane directions; default zero."),
    num("normal_damping", "Normal damping per area", Dim::Of(Quantity::SurfaceDamping))
        .help("Pressure per unit velocity in this direction; default zero."),
    num("tangential_damping", "Tangential damping per area", Dim::Of(Quantity::SurfaceDamping))
        .help("Pressure per unit velocity in this direction; default zero."),
    num("loss_factor", "Foundation loss factor", Dim::None)
        .help("Dimensionless structural loss multiplying this foundation stiffness; default zero."),
];
const DESCRIPTION: &str = r#"{"name":"Bed","set":"bottom","normal":1000,"tangential":0,"normal_damping":0,"tangential_damping":0,"loss_factor":0} — reference-surface local normal and isotropic in-plane stiffness per area, viscous damping per area and dimensionless loss factor. Numbers in display units or explicit units. Shell reference translations only; not unilateral contact. edit replaces the complete definition."#;
fn put(s: &mut Session, p: &Value, replace: bool) -> Result<Value> {
    let name = super::name(p)?;
    let mut fields = p.as_object().cloned().ok_or_else(|| EngineError::Params("foundation parameters must be an object".into()))?;
    fields.remove("name");
    let foundation: Foundation = serde_json::from_value(Value::Object(fields)).map_err(|e| EngineError::Params(e.to_string()))?;
    if replace {
        s.model.replace_foundation(name, foundation)?;
    } else {
        s.model.add_foundation(name, foundation)?;
    }
    Ok(json!({"name":name}))
}
pub const CREATE: CommandSpec = CommandSpec {
    id: "model.foundation.create",
    label: "New Surface Foundation",
    params: DESCRIPTION,
    form: Some(FORM),
    kind: CommandKind::Model,
    requires: &[Requirement::Set],
    run: |s, p| put(s, p, false),
};
pub const EDIT: CommandSpec = CommandSpec {
    id: "model.foundation.edit",
    label: "Edit Surface Foundation",
    params: DESCRIPTION,
    form: Some(FORM),
    kind: CommandKind::Model,
    requires: &[Requirement::Set],
    run: |s, p| put(s, p, true),
};
pub const DELETE: CommandSpec = CommandSpec {
    id: "model.foundation.delete",
    label: "Delete Surface Foundation",
    params: r#"{"name":"Bed"}"#,
    form: Some(&[text("name", "Name").req().help("The unique foundation definition name.")]),
    kind: CommandKind::Model,
    requires: &[],
    run: |s, p| {
        let name = super::name(p)?;
        s.model.remove_foundation(name)?;
        Ok(json!({"name":name}))
    },
};

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::{ElementFace, ElementKind, Model, Set};
    #[test]
    fn commands_convert_units_roundtrip_and_reject_invalid_edit() {
        let mut s = Session::default();
        for (i, p) in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]].into_iter().enumerate() {
            s.model.set_node(i as u32 + 1, p).unwrap();
        }
        s.model.set_element(1, ElementKind::Tri3, vec![1, 2, 3]).unwrap();
        s.model.add_set("bed", Set { faces: vec![ElementFace { element: 1, face: 1 }], ..Set::default() }).unwrap();
        s.execute("model.step.create", &json!({"name":"S","procedure":{"type":"static"}})).unwrap();
        let inputs = |s: &Session| serde_json::to_value(s.model.step_inputs("S").unwrap()).unwrap();
        let empty = inputs(&s);
        assert!(empty.get("foundations").is_none(), "preserve fingerprints of models without foundations");
        let params = json!({"name":"Bed","set":"bed","normal":"2 N/mm³","tangential":"1 N/mm³","normal_damping":"3 N*s/mm³"});
        s.execute("model.foundation.create", &params).unwrap();
        let f = s.model.foundation("Bed").unwrap();
        assert!((f.normal / 2e9 - 1.).abs() < 1e-14);
        assert!((f.normal_damping / 3e9 - 1.).abs() < 1e-14);
        let defined = inputs(&s);
        assert_ne!(defined, empty);
        let mut edited = params.clone();
        edited["normal"] = json!("4 N/mm³");
        s.execute("model.foundation.edit", &edited).unwrap();
        assert_ne!(inputs(&s), defined, "editing stiffness must invalidate results");
        s.execute("model.foundation.edit", &params).unwrap();
        let text = serde_json::to_string(&s.model).unwrap();
        let restored: Model = serde_json::from_str(&text).unwrap();
        assert_eq!(restored, s.model);
        let before = s.model.clone();
        let mut bad = params.clone();
        bad["normal"] = json!(-1.);
        assert!(s.execute("model.foundation.edit", &bad).is_err());
        assert_eq!(before, s.model);
        assert!(s.execute("model.foundation.create", &params).is_err());
        s.execute("model.foundation.delete", &json!({"name":"Bed"})).unwrap();
        assert!(s.model.foundations().is_empty());
        assert_eq!(inputs(&s), empty);
    }
}
