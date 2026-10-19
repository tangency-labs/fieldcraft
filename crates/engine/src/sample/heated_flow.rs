//! A small transient heated-flow demonstration using the original upright bracket.
//! Air at 20 °C passes an 80 °C bracket; no buoyancy or conjugate solid conduction.

use fieldcraft_geometry::Shape;
use serde_json::{Value, json};

pub(super) fn script(bodies: &[Shape]) -> Value {
    let walls: Vec<Value> = bodies
        .iter()
        .filter(|b| b.name == "Bracket")
        .flat_map(|b| {
            b.faces().iter().map(move |f| {
                json!({
                    "face": {"body": b.name, "face": f.name},
                    "condition": {"type": "temperature", "value": 80.0}
                })
            })
        })
        .collect();
    let [density, viscosity] = super::AIR;
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "model.material.create", "params": {"name": "Air", "density": density,
            "fluid": {"type": "newtonian", "viscosity": viscosity},
            "thermal": {"conductivity": "0.026 W/(m K)", "specific_heat": "1005 J/(kg K)"}}},
        {"command": "model.fluid_domain.create", "params": {
            "name": "WarmAir", "material": "Air", "min": [-30, -20, 0], "max": [100, 60, 70], "spacing": 5,
            "faces": {"x_min": {"type": "velocity", "velocity": [1, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0},
                "y_min": {"type": "slip"}, "y_max": {"type": "slip"}, "z_min": {"type": "wall"}, "z_max": {"type": "slip"}},
            "obstacles": ["Bracket"], "initial_velocity": [1, 0, 0],
            "thermal": {"initial_temperature": 20, "surfaces": walls,
                "faces": {"x_min": {"type": "temperature", "value": 20}}}
        }},
        {"command": "model.step.create", "params": {"name": "HeatedFlow", "procedure": {
            "type": "fluid", "domain": "WarmAir", "steps": 180, "frames": 6
        }}}
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn heated_flow_sample_solves_and_exposes_temperature_and_wall_heat() {
        let mut s = crate::Session::default();
        s.execute("file.sample", &json!({"name": "heated-flow"})).unwrap();
        s.execute("solve.fluid", &json!({"backend": "cpu"})).unwrap();
        let report = s.execute("results.slice", &json!({"field": "temperature"})).unwrap();
        let max = report["contour"]["range"]["max"]["value"].as_f64().unwrap();
        let min = report["contour"]["range"]["min"]["value"].as_f64().unwrap();
        assert!(max > 20.1 && max <= 80.01 && min >= 19.99, "{report}");
        let report = s.execute("results.surface", &json!({"field": "wall_heat_flux"})).unwrap();
        assert_eq!(report["contour"]["unit"], "mW/mm²");
        assert!(report["wall_heat"]["faces"].as_array().unwrap().iter().any(|f| f["power"].as_f64().unwrap_or(0.0) < 0.0));
        assert_eq!(report["frames"].as_array().unwrap().len(), 7);
    }
}
