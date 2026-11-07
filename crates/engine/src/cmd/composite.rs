//! Inspection of material definitions and laminate section properties. Constitutive and
//! failure calculations belong to `fem::laminate`; this module only resolves names and units.

use fieldcraft_fem::laminate::{Laminate, resolve_plies};
use fieldcraft_model::SectionKind;
use fieldcraft_model::units::{DisplaySystem, Quantity as Q};
use serde::Deserialize;
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Named {
    name: String,
}

fn named(p: &Value) -> Result<Named> {
    let named: Named = serde_json::from_value(p.clone()).map_err(|e| EngineError::Params(e.to_string()))?;
    fieldcraft_model::check_name("definition", &named.name)?;
    Ok(named)
}

pub const MATERIAL_GET: CommandSpec = CommandSpec {
    id: "model.material.get",
    label: "Inspect Material",
    params: r#"{"name": "Carbon"} — the named material, including its optional ply strengths, in the display units; usable as model.material.edit parameters"#,
    form: Some(super::forms::MATERIAL_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Material],
    run: |s, p| {
        let n = named(p)?.name;
        let m = s.model.material(&n).ok_or_else(|| EngineError::Params(format!("material `{n}` does not exist")))?;
        definition(s, "model.material.create", &n, m)
    },
};

pub const SECTION_GET: CommandSpec = CommandSpec {
    id: "model.section.get",
    label: "Inspect Section",
    params: r#"{"name": "Skin"} — the named section, including its material orientation and bottom-to-top layup, in the display units; usable as model.section.edit parameters"#,
    form: Some(super::forms::SECTION_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Section],
    run: |s, p| {
        let n = named(p)?.name;
        let section = s.model.section(&n).ok_or_else(|| EngineError::Params(format!("section `{n}` does not exist")))?;
        definition(s, "model.section.create", &n, section)
    },
};

fn definition(s: &Session, command: &str, name: &str, value: &impl serde::Serialize) -> Result<Value> {
    let mut v = serde_json::to_value(value).map_err(|e| EngineError::Params(e.to_string()))?;
    if let Some(o) = v.as_object_mut() {
        o.insert("name".into(), json!(name));
    }
    crate::units::shown_params(command, &v, s.units())
}

pub const SECTION_PROPERTIES: CommandSpec = CommandSpec {
    id: "model.section.properties",
    label: "Laminate Properties",
    params: r#"{"name": "Skin"} — a composite shell's ABD blocks (a: force/length, b: force, d: moment), transverse shear matrix (force/length), thickness, centre z and mass moments [∫ρ dz, ∫ρ z dz, ∫ρ z² dz], about its nodal reference surface. Matrices use engineering in-plane shear and the section reference axes before element projection; units are explicit. Missing ply density leaves mass_moments null and lists the materials in missing_density; it never assumes zero mass"#,
    form: Some(super::forms::SECTION_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Section],
    run: properties,
};

fn shown(display: &DisplaySystem, q: Q, value: f64) -> Result<f64> {
    display.express(q, value).map_err(|e| EngineError::Params(e.to_string()))
}

fn properties(s: &mut Session, p: &Value) -> Result<Value> {
    let n = named(p)?.name;
    let section = s.model.section(&n).ok_or_else(|| EngineError::Params(format!("section `{n}` does not exist")))?;
    let SectionKind::CompositeShell { layup } = &section.kind else {
        return Err(EngineError::Params(format!("section `{n}` is not a composite shell")));
    };
    let plies = resolve_plies(layup, &s.model).map_err(|e| EngineError::Params(e.to_string()))?;
    let laminate = Laminate::new(&plies, layup.offset).map_err(|e| EngineError::Params(e.to_string()))?;
    let display = s.units();
    let abd = laminate.abd();
    let block = |row: usize, col: usize, q: Q| -> Result<Vec<Vec<f64>>> {
        (0..3).map(|i| (0..3).map(|j| shown(display, q, abd[i + row][j + col])).collect()).collect()
    };
    let shear: Vec<Vec<f64>> =
        laminate.shear().iter().map(|row| row.iter().map(|&v| shown(display, Q::ForcePerLength, v)).collect()).collect::<Result<_>>()?;
    let missing_density: Vec<&str> = layup.plies.iter().zip(&plies).filter(|(_, p)| p.density.is_none()).map(|(p, _)| p.material.as_str()).collect();
    // These compound mass quantities have no independent display override. Derive their
    // factors from the model unit table, never from hand-written conversion constants.
    let length = display.to_si(Q::Length, 1.0).map_err(|e| EngineError::Params(e.to_string()))?;
    let mass_moments = if missing_density.is_empty() {
        let moments = laminate.mass_moments().map_err(|e| EngineError::Params(e.to_string()))?;
        let values = [
            shown(display, Q::Mass, moments[0])? * length * length,
            shown(display, Q::Mass, moments[1])? * length,
            shown(display, Q::Mass, moments[2])?,
        ];
        if values.iter().any(|v| !v.is_finite()) {
            return Err(EngineError::Params("laminate mass moments overflow the display units".into()));
        }
        Some(values)
    } else {
        None
    };
    Ok(json!({
        "name": n,
        "frame": "section_reference",
        "abd": {"a": block(0, 0, Q::ForcePerLength)?, "b": block(0, 3, Q::Force)?, "d": block(3, 3, Q::Moment)?},
        "shear": shear,
        "thickness": shown(display, Q::Length, laminate.thickness())?,
        "centre_z": shown(display, Q::Length, laminate.eccentricity())?,
        "offset": layup.offset,
        "mass_moments": mass_moments,
        "missing_density": missing_density,
        "units": {
            "a": display.symbol(Q::ForcePerLength), "b": display.symbol(Q::Force), "d": display.symbol(Q::Moment),
            "shear": display.symbol(Q::ForcePerLength), "thickness": display.symbol(Q::Length), "centre_z": display.symbol(Q::Length),
            "mass_moments": [format!("({})/({})²", display.symbol(Q::Mass), display.symbol(Q::Length)), format!("({})/({})", display.symbol(Q::Mass), display.symbol(Q::Length)), display.symbol(Q::Mass)],
        },
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model() -> Session {
        let mut s = Session::default();
        s.execute("model.units", &json!({"system": "mm-t-s-N-MPa", "overrides": {"angle": "deg"}})).unwrap();
        for (i, xyz) in [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]].into_iter().enumerate() {
            s.execute("mesh.node.add", &json!({"label": i + 1, "xyz": xyz})).unwrap();
        }
        s.execute("mesh.element.add", &json!({"label": 1, "kind": "quad4", "nodes": [1, 2, 3, 4]})).unwrap();
        s.execute("model.set.create", &json!({"name": "plate", "elements": [1]})).unwrap();
        s.execute(
            "model.material.create",
            &json!({
                "name": "Ply", "elastic": {"type": "isotropic", "young": "120 GPa", "poisson": 0.25},
                "density": "2000 kg/m^3",
                "ply_strength": {"tension": ["1 GPa", 40], "compression": [800, {"value": 120, "unit": "MPa"}], "shear12": "50 MPa", "shear23": 60}
            }),
        )
        .unwrap();
        s.execute("model.section.create", &json!({
            "name": "Skin", "set": "plate", "kind": {"type": "composite_shell", "layup": {"plies": [{"material": "Ply", "thickness": {"value": 2, "unit": "mm"}, "angle": "45 deg", "points": 2}], "offset": 0.25}},
            "orientation": {"frame": {"type": "cylindrical", "origin": ["2 cm", 0, 0], "axis": [0, 0, 1]}, "shell_reference_axis": 2}
        })).unwrap();
        s
    }

    fn close(value: &Value, expected: f64) {
        let actual = value.as_f64().unwrap();
        assert!((actual - expected).abs() <= 1e-12 * expected.abs().max(1e-12), "{actual} != {expected}");
    }

    #[test]
    fn composite_definition_units_round_trip_and_reject_unknown_nested_keys() {
        let mut s = model();
        let material = s.execute("model.material.get", &json!({"name": "Ply"})).unwrap();
        assert_eq!(material["ply_strength"]["tension"], json!([1000.0, 40.0]));
        assert_eq!(material["ply_strength"]["compression"], json!([800.0, 120.0]));
        close(&material["density"], 2e-9);
        let section = s.execute("model.section.get", &json!({"name": "Skin"})).unwrap();
        close(&section["kind"]["layup"]["plies"][0]["angle"], 45.0);
        close(&section["kind"]["layup"]["plies"][0]["thickness"], 2.0);
        assert_eq!(section["orientation"]["frame"]["origin"], json!([20.0, 0.0, 0.0]));
        let before = s.model.clone();
        s.execute("model.material.edit", &material).unwrap();
        s.execute("model.section.edit", &section).unwrap();
        assert_eq!(s.model.material("Ply"), before.material("Ply"));
        assert_eq!(s.model.section("Skin"), before.section("Skin"));
        for pointer in ["/kind/layup/plies/0/typo", "/kind/layup/typo", "/orientation/frame/typo", "/orientation/typo"] {
            let mut bad = section.clone();
            let (parent, key) = pointer.rsplit_once('/').unwrap();
            bad.pointer_mut(parent).unwrap()[key] = json!(1);
            assert!(s.execute("model.section.edit", &bad).is_err(), "{pointer}");
        }
        let mut bad = material.clone();
        bad["ply_strength"]["shear13"] = json!(1);
        assert!(s.execute("model.material.edit", &bad).is_err());
        let mut bad = section.clone();
        bad["kind"]["layup"]["plies"][0]["thickness"] = json!("1 MPa");
        assert!(s.execute("model.section.edit", &bad).is_err());
        assert_eq!(s.model.section("Skin"), before.section("Skin"));
        assert!(s.execute("model.section.properties", &json!({"name": "Skin", "typo": 1})).is_err());
    }

    #[test]
    fn composite_properties_follow_offset_and_display_units_and_report_missing_density() {
        let mut s = model();
        let p = s.execute("model.section.properties", &json!({"name": "Skin"})).unwrap();
        close(&p["thickness"], 2.0);
        close(&p["centre_z"], -0.5);
        close(&p["abd"]["a"][0][0], 256000.0);
        close(&p["abd"]["b"][0][0], -128000.0);
        close(&p["abd"]["d"][0][0], 256000.0 * (1.0 / 3.0 + 0.25));
        close(&p["shear"][0][0], 80000.0);
        close(&p["mass_moments"][0], 4e-9);
        close(&p["mass_moments"][1], -2e-9);
        close(&p["mass_moments"][2], 7e-9 / 3.0);
        assert_eq!(p["units"]["a"], "N/mm");
        assert_eq!(p["units"]["d"], "N·mm");
        assert_eq!(p["frame"], "section_reference");
        let mut material = s.execute("model.material.get", &json!({"name": "Ply"})).unwrap();
        material.as_object_mut().unwrap().remove("density");
        s.execute("model.material.edit", &material).unwrap();
        let p = s.execute("model.section.properties", &json!({"name": "Skin"})).unwrap();
        assert!(p["mass_moments"].is_null());
        assert_eq!(p["missing_density"], json!(["Ply"]));
        close(&p["abd"]["a"][0][0], 256000.0);
        s.execute("model.units", &json!({"system": "m-kg-s-N-Pa"})).unwrap();
        let p = s.execute("model.section.properties", &json!({"name": "Skin"})).unwrap();
        close(&p["thickness"], 0.002);
        close(&p["abd"]["a"][0][0], 256000000.0);
        close(&p["abd"]["b"][0][0], -128000.0);
    }
}
