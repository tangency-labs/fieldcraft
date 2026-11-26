//! `deck.openfoam` with a `compressible` gas: the case it writes, the conditions it refuses, and
//! the gas its marker records.

use fieldcraft_model::units::Quantity;
use serde_json::{Value, json};

use super::{MARKER, Marker, MarkerGas};
use crate::Session;

/// A value of `q` given in SI, in the session's display units, as a command parameter takes it.
fn show(s: &Session, q: Quantity, si: f64) -> f64 {
    s.units().express(q, si).unwrap()
}

/// The session of a plate's sample design, its box holding air: a pressure face at 2 bar gauge on
/// the inlet side and the ambient pressure on the outlet side, the plate as the obstacle.
struct Air {
    _guard: crate::tests::TempDir,
    s: Session,
    dir: std::path::PathBuf,
    min: Vec<f64>,
    max: Vec<f64>,
}

/// `thermal`: the material has its conductivity and specific heat. `turbulence`: the step has a
/// Smagorinsky model.
fn air(thermal: bool, turbulence: bool) -> Air {
    let (guard, design) = crate::tests::sample_design_file(&fieldcraft_geometry::sample::Plate::default());
    let dir = std::path::Path::new(&design).parent().unwrap().to_path_buf();
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": design})).unwrap();
    let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
    let [lo, hi] = [0, 1].map(|k| bodies[0]["bbox"][k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>());
    let min: Vec<f64> = lo.iter().map(|v| (v / 5.0).floor() * 5.0 - 20.0).collect();
    let max: Vec<f64> = hi.iter().map(|v| (v / 5.0).ceil() * 5.0 + 20.0).collect();
    let mut material = json!({
        "name": "Air",
        "density": show(&s, Quantity::Density, 1.2),
        "fluid": {"type": "newtonian", "viscosity": show(&s, Quantity::DynamicViscosity, 1.8e-5)},
    });
    if thermal {
        material["thermal"] = json!({
            "conductivity": show(&s, Quantity::ThermalConductivity, 0.0257),
            "specific_heat": show(&s, Quantity::SpecificHeat, 1005.0),
        });
    }
    let mut procedure = json!({"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10});
    if turbulence {
        procedure["turbulence"] = json!({"type": "smagorinsky", "constant": 0.1, "prandtl": 0.9});
    }
    let script = json!([
        {"command": "model.material.create", "params": material},
        {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Air", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"],
            "faces": {"x_min": {"type": "pressure", "pressure": show(&s, Quantity::Pressure, 2.0e5)}, "x_max": {"type": "pressure", "pressure": 0}}}},
        {"command": "model.step.create", "params": {"name": "Flow", "procedure": procedure}},
    ]);
    s.run_script(&script).unwrap();
    let dir = dir.join("case-root");
    std::fs::create_dir_all(&dir).unwrap();
    Air { _guard: guard, s, dir, min, max }
}

/// The compressible conditions of [`air`], each value in the session's display units.
fn gas_json(s: &Session) -> Value {
    json!({
        "ambient": show(s, Quantity::Pressure, 101_325.0),
        "temperature": show(s, Quantity::Temperature, 293.15),
        "gas_constant": show(s, Quantity::SpecificHeat, 287.05),
        "far_field": show(s, Quantity::Length, 0.05),
    })
}

/// The refusal `deck.openfoam` gives for `extra` added to the path of case `name`, a new directory.
fn refusal(a: &mut Air, name: &str, extra: Value) -> String {
    let mut params = json!({"path": a.dir.join(name).to_str().unwrap()});
    params.as_object_mut().unwrap().extend(extra.as_object().unwrap().clone());
    a.s.execute("deck.openfoam", &params).unwrap_err().to_string()
}

fn close(got: f64, want: f64, relative: f64) -> bool {
    (got - want).abs() <= relative * want.abs()
}

#[test]
fn a_compressible_case_is_written_with_its_gas_recorded() {
    let mut a = air(true, false);
    let gas = gas_json(&a.s);
    let case = a.dir.join("case");
    let r = a.s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "compressible": gas})).unwrap();
    assert_eq!(r["solver"], "rhoCentralFoam");
    for file in ["0/T", "constant/thermophysicalProperties", "system/fvSchemes"] {
        assert!(case.join(file).is_file(), "{file} is written");
    }
    for file in ["0/nut", "constant/transportProperties"] {
        assert!(!case.join(file).exists(), "{file} is not written for a gas");
    }
    let control = std::fs::read_to_string(case.join("system/controlDict")).unwrap();
    assert!(control.contains("maxCo 2e-1;"), "the default Courant limit is 0.2, the CFL number of Greenshields et al. 2010 §4.1: {control}");
    let gamma = r["gas"]["gamma"].as_f64().unwrap();
    assert!((gamma - 1005.0 / (1005.0 - 287.05)).abs() < 1e-9, "gamma {gamma}");
    let marker: Value = serde_json::from_str(&std::fs::read_to_string(case.join(MARKER)).unwrap()).unwrap();
    assert_eq!(marker["solver"], "rhoCentralFoam");
    for (key, want) in [("ambient", 101_325.0), ("temperature", 293.15), ("gas_constant", 287.05), ("specific_heat", 1005.0)] {
        let got = marker["gas"][key].as_f64().unwrap();
        assert!(close(got, want, 1e-12), "gas.{key} = {got}, want {want}");
    }
}

#[test]
fn nothing_about_the_gas_is_defaulted() {
    let mut a = air(true, false);
    let gas = gas_json(&a.s);
    let without = |key: &str| {
        let mut g = gas.clone();
        g.as_object_mut().unwrap().remove(key);
        json!({"compressible": g})
    };
    for key in ["gas_constant", "ambient", "temperature", "far_field"] {
        let e = refusal(&mut a, &format!("no-{key}"), without(key));
        assert!(e.contains(key), "a missing {key} is refused naming it: {e}");
    }
    let e = refusal(&mut a, "not-an-object", json!({"compressible": true}));
    assert!(e.contains("compressible"), "{e}");
    let e = refusal(&mut a, "courant", json!({"compressible": gas, "max_courant": 0.6}));
    assert!(e.contains("Courant"), "{e}");
    let zero = show(&a.s, Quantity::Temperature, 0.0);
    let mut cold = gas.clone();
    cold["temperature"] = json!(zero);
    let e = refusal(&mut a, "zero-kelvin", json!({"compressible": cold}));
    assert!(e.contains("temperature"), "0 K is refused: {e}");
    let case = a.dir.join("courant-half");
    a.s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "compressible": gas, "max_courant": 0.5})).unwrap();
    let control = std::fs::read_to_string(case.join("system/controlDict")).unwrap();
    assert!(control.contains("maxCo 5e-1;"), "a Courant limit of 1/2 is accepted: {control}");
}

#[test]
fn a_material_without_its_thermal_properties_is_refused_for_a_gas() {
    let mut a = air(false, false);
    let gas = gas_json(&a.s);
    let e = refusal(&mut a, "no-thermal", json!({"compressible": gas}));
    assert!(e.contains("specific heat"), "{e}");
}

#[test]
fn initial_pressure_boxes_run_set_fields() {
    let mut a = air(true, false);
    let gas = gas_json(&a.s);
    let inside = |f: f64| a.min.iter().zip(&a.max).map(|(l, h)| l + f * (h - l)).collect::<Vec<f64>>();
    let (lo, hi) = (inside(0.25), inside(0.75));
    let boxed = a.dir.join("boxed");
    let mut boxes = gas.clone();
    boxes["initial_pressure"] = json!([{"min": lo, "max": hi, "pressure": show(&a.s, Quantity::Pressure, 1.0e6)}]);
    a.s.execute("deck.openfoam", &json!({"path": boxed.to_str().unwrap(), "compressible": boxes})).unwrap();
    assert!(boxed.join("system/setFieldsDict").is_file(), "the boxes are set by setFields");
    assert!(std::fs::read_to_string(boxed.join("Allrun")).unwrap().contains("run setFields"));
    let plain = a.dir.join("plain");
    a.s.execute("deck.openfoam", &json!({"path": plain.to_str().unwrap(), "compressible": gas})).unwrap();
    assert!(!plain.join("system/setFieldsDict").exists(), "no boxes, no setFields");
    assert!(!std::fs::read_to_string(plain.join("Allrun")).unwrap().contains("run setFields"));
}

#[test]
fn a_turbulence_model_is_refused_by_name_in_a_compressible_case() {
    let mut a = air(true, true);
    let gas = gas_json(&a.s);
    let e = refusal(&mut a, "turbulent", json!({"compressible": gas}));
    assert!(e.contains("the smagorinsky turbulence model in a compressible case"), "{e}");
}

#[test]
fn an_incompressible_marker_has_no_gas_and_a_mismatch_is_refused() {
    let mut a = air(true, false);
    let case = a.dir.join("incompressible");
    let r = a.s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap()})).unwrap();
    assert!(r["gas"].is_null(), "an incompressible case reports no gas");
    let marker: Value = serde_json::from_str(&std::fs::read_to_string(case.join(MARKER)).unwrap()).unwrap();
    assert!(marker.get("gas").is_none(), "an incompressible marker has no gas key");

    let gas = || MarkerGas { ambient: 101_325.0, temperature: 293.15, gas_constant: 287.05, specific_heat: 1005.0 };
    let marker = |solver: &str, gas: Option<MarkerGas>| Marker {
        step: "Flow".into(),
        solver: solver.into(),
        density: 1.2,
        background: 0.005,
        spacing: 0.005,
        time_precision: 6,
        gas,
    };
    assert!(marker("rhoCentralFoam", Some(gas())).check().is_ok());
    assert!(marker("rhoCentralFoam", None).check().is_err(), "a compressible solver without its gas");
    assert!(marker("pimpleFoam", Some(gas())).check().is_err(), "an incompressible solver with a gas");
    assert!(marker("simpleFoam", None).check().is_err(), "a solver deck.openfoam does not write");
    let weak = MarkerGas { specific_heat: 200.0, gas_constant: 287.0, ..gas() };
    assert!(marker("rhoCentralFoam", Some(weak)).check().is_err(), "cp below R");
}
