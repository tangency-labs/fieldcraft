//! Saved-frame measurements for cooling studies. Port integrals use the solver's existing
//! reference-density convention; wall powers use recorded exchange, without requiring geometry.

use fieldcraft_model::units::Quantity as Q;
use serde_json::{Value, json};

use fieldcraft_model::BoxFace;

use super::results::{Cell, Table, Want, fluid_frame, step_for};
use crate::results::{Pick, shown, unit_of};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult, requirement::Requirement};

pub const MEASURE: CommandSpec = CommandSpec {
    id: "results.fluid.measure",
    label: "Fluid Measurements",
    params: r#"{"step": "Flow", "frame": 10, "from": "x_min", "to": "x_max"} — read-only measurements of one saved fluid frame (default the last): every open box face's fluid area, area-mean static gauge pressure, temperature and velocity, separate inflow/outflow volume and reference-density mass rates, mass-flow-weighted temperatures and advective sensible heat relative to the mixed inflow temperature; every recorded box/obstacle wall's heat rate, positive into the wall, interpolated at the frame time without extrapolation. Optional `from` and `to` give their area-mean static pressure difference. Returns `rows` with face/body, quantity, value and unit, and the same table as `csv` including step/frame/time. Null means unavailable; a blocked face has no mean. Instantaneous saved-frame measurements, not time averages or a validation verdict; no port conduction, arbitrary surfaces, exchanger effectiveness or complete energy balance yet. Stale results are refused"#,
    form: Some(super::forms::FLUID_MEASURE),
    kind: CommandKind::Session,
    requires: &[Requirement::Result],
    run: measure,
};

/// Measures one saved frame of a current fluid result: the face, wall and pressure-drop rows.
fn measure(s: &mut Session, p: &Value) -> Result<Value> {
    let step = step_for(s, p, Want::Fluid)?;
    let result = super::path::current_result(s, &step)?;
    let StepResult::Fluid(r) = result.as_ref() else {
        return Err(EngineError::Params("measurements require a fluid result".into()));
    };
    if r.surface.is_some() {
        return Err(EngineError::Params("fluid measurements do not support free-surface results yet".into()));
    }
    let frame = crate::flow::frame_index(r, Pick { frame: fluid_frame(r, p)?, ..Pick::default() })?;
    let frame_time = r.frames.get(frame).map(|f| f.time).ok_or_else(|| EngineError::invalid("frame", "no saved fluid frame at this index"))?;
    let ports = r.ports(frame)?;
    let pressure_drop = match (super::opt::<String>(p, "from")?, super::opt::<String>(p, "to")?) {
        (Some(from), Some(to)) => {
            let pressure = |key: &str, name: &str| {
                let face = BoxFace::ALL
                    .into_iter()
                    .find(|f| f.name() == name)
                    .ok_or_else(|| EngineError::invalid(key, "use x_min, x_max, y_min, y_max, z_min or z_max"))?;
                ports
                    .iter()
                    .find(|port| port.face == face)
                    .and_then(|port| port.mean)
                    .map(|m| m.pressure)
                    .ok_or_else(|| EngineError::invalid(key, format!("face `{name}` has no open fluid area to measure")))
            };
            Some((format!("{from} -> {to}"), pressure("from", &from)? - pressure("to", &to)?))
        }
        (None, None) => None,
        _ => return Err(EngineError::Params("pressure drop requires both `from` and `to`".into())),
    };
    let display = s.units();
    let time = shown(Q::Time, frame_time, display);
    let time_unit = unit_of(Q::Time, display);
    let text = |t: Option<&str>| t.map_or(Cell::Empty, |t| Cell::Text(t.into()));
    let mut rows = Vec::new();
    let mut table = Table {
        header: ["step", "frame", "time", "time_unit", "body", "face", "quantity", "value", "unit"].map(String::from).to_vec(),
        rows: Vec::new(),
    };
    let mut add = |body: Option<&str>, face: &str, quantity: &str, q: Q, value: Option<f64>| {
        let value = value.map(|v| shown(q, v, display)).filter(|v| v.is_finite());
        let unit = unit_of(q, display);
        rows.push(json!({"body": body, "face": face, "quantity": quantity, "value": value, "unit": unit}));
        table.rows.push(vec![
            Cell::Text(step.clone()),
            Cell::Number(frame as f64),
            Cell::Number(time),
            text(time_unit.as_deref()),
            text(body),
            Cell::Text(face.into()),
            Cell::Text(quantity.into()),
            value.map_or(Cell::Empty, Cell::Number),
            text(unit.as_deref()),
        ]);
    };
    for port in &ports {
        let face = port.face.name();
        let mean = port.mean;
        add(None, face, "area", Q::Area, Some(port.area));
        add(None, face, "area_mean_pressure", Q::Pressure, mean.map(|m| m.pressure));
        add(None, face, "area_mean_temperature", Q::Temperature, mean.and_then(|m| m.temperature));
        for (axis, name) in ["area_mean_velocity_x", "area_mean_velocity_y", "area_mean_velocity_z"].into_iter().enumerate() {
            add(None, face, name, Q::Velocity, mean.and_then(|m| m.velocity.get(axis).copied()));
        }
        for (direction, stream) in [("inflow", port.inflow), ("outflow", port.outflow)] {
            for (name, q, value) in [
                ("volume_flow", Q::VolumeFlow, Some(stream.volume_flow)),
                ("mass_flow", Q::MassFlow, stream.mass_flow),
                ("mean_temperature", Q::Temperature, stream.temperature),
                ("sensible_heat", Q::Power, stream.sensible_heat),
            ] {
                add(None, face, &format!("{direction}_{name}"), q, value);
            }
        }
    }
    add(None, "all_open_faces", "heat_reference", Q::Temperature, fieldcraft_lbm::ports::inlet_temperature(&ports));
    for wall in &r.heat {
        add(wall.body.as_deref(), &wall.name, "wall_heat_rate", Q::Power, r.wall_power(wall, frame_time));
    }
    if let Some((faces, drop)) = &pressure_drop {
        add(None, faces, "pressure_drop", Q::Pressure, Some(*drop));
    }
    Ok(json!({
        "step": step, "frame": frame, "time": time, "time_unit": time_unit,
        "time_basis": "one saved frame; wall heat interpolated between recorded samples, without extrapolation",
        "pressure_drop": pressure_drop.map(|(_, drop)| shown(Q::Pressure, drop, display)),
        "pressure_unit": unit_of(Q::Pressure, display),
        "conventions": {
            "face_means": "fluid cells on the saved open boundary layer; area-weighted pressure, temperature and velocity",
            "flow": "separate nonnegative inflow and outflow; reference-density mass flow; mass-flow-weighted temperatures",
            "sensible_heat": "advective sensible heat relative to heat_reference (mixed inflow temperature); excludes port conduction",
            "wall_heat_rate": "positive into the wall; recorded thermal exchange, separate from port advection",
            "scope": "measurements only; no time averaging or complete energy balance"
        },
        "rows": rows, "csv": table.csv(),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn session() -> Session {
        let mut s = Session::default();
        let mut r = crate::flow::tests::channel();
        r.open_faces = [true, true, false, false, false, false];
        r.port_properties = Some(fieldcraft_lbm::PortProperties { density: 2.0, specific_heat: Some(5.0) });
        r.temperature_reference = fieldcraft_model::units::celsius_to_kelvin(0.0);
        let f = r.frames.last_mut().unwrap();
        f.velocity.fill([1.0, 0.0, 0.0]);
        f.temperature = (0..r.grid.len()).map(|c| if r.grid.coords(c)[0] == 0 { 20.0 } else { 30.0 }).collect();
        r.history_times = vec![0.5, 1.5];
        r.heat = vec![fieldcraft_lbm::HeatHistory { name: "heater,\"face\"".into(), body: Some("Block".into()), power: vec![-1200.0, -2000.0] }];
        s.store_result(StepResult::Fluid(Box::new(r)));
        s.results_view = None;
        s
    }

    fn value(report: &Value, face: &str, quantity: &str) -> Value {
        report["rows"].as_array().unwrap().iter().find(|r| r["face"] == face && r["quantity"] == quantity).unwrap()["value"].clone()
    }

    #[test]
    fn measurements_report_analytic_heat_pressure_units_and_csv_without_changing_view() {
        let mut s = session();
        let journal = s.journal.len();
        let report = s.execute(MEASURE.id, &json!({"from": "x_min", "to": "x_max"})).unwrap();
        // Default mm-t-s-N-MPa display: pressure 10 - 1 Pa, area 16 m²,
        // flow 16 m³/s, mass 32 kg/s, heat rho cp Q (Tout - Tin) = 1600 W.
        let display = s.units().clone();
        let shown = |q, x| json!(display.express(q, x).unwrap());
        assert_eq!(report["pressure_drop"], shown(Q::Pressure, 9.0));
        assert_eq!(value(&report, "x_min", "area"), shown(Q::Area, 16.0));
        assert_eq!(value(&report, "x_max", "outflow_mass_flow"), shown(Q::MassFlow, 32.0));
        assert_eq!(value(&report, "x_max", "outflow_sensible_heat"), shown(Q::Power, 1600.0));
        assert_eq!(value(&report, "heater,\"face\"", "wall_heat_rate"), shown(Q::Power, -1600.0));
        assert_eq!(value(&report, "x_max", "area_mean_temperature"), json!(30.0));
        assert_eq!(value(&report, "x_min", "outflow_mean_temperature"), Value::Null);
        // The solve reply and results view report the same area means.
        let StepResult::Fluid(r) = s.latest_result("Flow").unwrap().0.as_ref() else { panic!() };
        let ports = crate::results::fluid_ports(r, report["frame"].as_u64().unwrap() as usize, s.units());
        assert_eq!(ports["faces"][1]["area_mean"]["temperature"], value(&report, "x_max", "area_mean_temperature"));
        assert_eq!(ports["faces"][0]["area_mean"]["pressure"], value(&report, "x_min", "area_mean_pressure"));
        assert_eq!(ports["units"]["velocity"], json!(unit_of(Q::Velocity, s.units())));
        let csv = report["csv"].as_str().unwrap();
        assert!(csv.starts_with("step,frame,time,time_unit,body,face,quantity,value,unit\r\n"));
        assert!(csv.contains("\"heater,\"\"face\"\"\""));
        assert_eq!(csv.lines().count(), report["rows"].as_array().unwrap().len() + 1);
        assert_eq!(s.journal.len(), journal);
        assert!(s.results_view.is_none());
        let earlier = s.execute(MEASURE.id, &json!({"frame": 0, "from": "x_max", "to": "x_min"})).unwrap();
        assert_eq!(earlier["pressure_drop"], 0.0);
        assert_eq!(earlier["time"], 0.0);
        assert!(value(&earlier, "heater,\"face\"", "wall_heat_rate").is_null());
        assert!(value(&earlier, "x_min", "area_mean_temperature").is_null());
        let reversed = s.execute(MEASURE.id, &json!({"from": "x_max", "to": "x_min"})).unwrap();
        assert_eq!(reversed["pressure_drop"], shown(Q::Pressure, -9.0));
    }

    #[test]
    fn measurements_refuse_bad_parameters_incomplete_and_stale_results() {
        let mut s = session();
        for p in [
            json!({"frame": -1}),
            json!({"frame": u64::MAX}),
            json!({"frame": 2}),
            json!({"from": "x_min"}),
            json!({"from": "y_min", "to": "x_max"}),
            json!({"from": "unknown", "to": "x_max"}),
            json!({"step": "missing"}),
        ] {
            assert!(s.execute(MEASURE.id, &p).is_err(), "{p}");
        }
        let StepResult::Fluid(mut r) = s.latest_result("Flow").unwrap().0.as_ref().clone() else { panic!() };
        r.frames.last_mut().unwrap().pressure.clear();
        s.store_result(StepResult::Fluid(r));
        assert!(s.execute(MEASURE.id, &json!({})).is_err());
        let mut s = session();
        s.execute("mesh.node.add", &json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        assert!(s.execute(MEASURE.id, &json!({})).unwrap_err().to_string().contains("stale"));
    }
}
