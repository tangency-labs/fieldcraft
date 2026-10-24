//! `results.path`: exact solid shape-function interpolation along reference-coordinate paths.

use std::sync::Arc;

use serde_json::{Value, json};

use super::results::{Cell, Table, Want, step_for};
use crate::path::{MAX_POINTS, MAX_SAMPLES};
use crate::results::{Field, Kind, Pick, Quantity, field_values};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};

pub const PATH: CommandSpec = CommandSpec {
    id: "results.path",
    label: "Path Plot",
    params: r#"{"step":"Load","field":"stress_xx","points":[[20,5,20],[80,5,20]],"samples":101,"frame":null,"mode":null,"path":"top-fibre.csv"} — interpolate the existing nodal scalar field with the containing solid element's shape functions along 2 to 1000 reference-coordinate polyline points, at equally spaced arc-length stations (2 to 10000 samples, default 101). Supports tet4/tet10/wedge6/pyramid5/hex8/hex20, including curved quadratic geometry. Outside elements gives null, never nearest-node values; unresolved inverse searches are errors. Shared faces use the lowest element label. Stale results, unavailable intermediate fields and categorical contact_status are refused. Nonlinear stored frames permit displacements only; final stress/contact/temperature fields require omitted/null frame. `frame` is from 0 and `mode` from 1; omitted/null frame selects the final state, and omitted/null mode selects mode 1. Optional CSV is distance,x,y,z,value in the same display units as the reply, with empty values at gaps. Replies {step,field,distance,xyz,values,elements,units}. The index is cached per exact result/definition; caps are 500000 mesh elements, 2000000 result nodes, 96 MiB index allocations, 4096 candidates per station, 2000000 total BVH visits, 14000000 separately capped construction basis/enclosure work units and 4000000 query basis/enclosure work units. Construction and sampling report progress and support cancellation."#,
    form: Some(super::forms::RESULTS_PATH),
    kind: CommandKind::Session,
    requires: &[crate::requirement::Requirement::Result],
    run: path,
};

fn optional_index(p: &Value, key: &str) -> Result<Option<usize>> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => v
            .as_u64()
            .and_then(|v| usize::try_from(v).ok())
            .map(Some)
            .ok_or_else(|| EngineError::invalid(key, "must be a nonnegative integer or null")),
    }
}

fn pick(result: &StepResult, field: Field, p: &Value) -> Result<Pick> {
    if field == Field::ContactStatus {
        return Err(EngineError::invalid(
            "field",
            "contact_status is categorical; interpolating status codes would invent fractional states (use solve.probe at nodes)",
        ));
    }
    let frame = optional_index(p, "frame")?;
    let mode = optional_index(p, "mode")?;
    let pick = Pick { frame, mode: mode.unwrap_or(1) };
    let displacement = matches!(field, Field::Displacement | Field::DisplacementX | Field::DisplacementY | Field::DisplacementZ);
    match result.kind() {
        Kind::Frequency(f) => {
            if frame.is_some() {
                return Err(EngineError::invalid("frame", "frequency results have modes, not frames"));
            }
            if pick.mode == 0 || pick.mode > f.modes.len() {
                return Err(EngineError::invalid("mode", format!("must be from 1 to {}", f.modes.len())));
            }
        }
        Kind::Static(s) => {
            if mode.is_some() {
                return Err(EngineError::invalid("mode", "only frequency results have modes"));
            }
            if let Some(k) = frame {
                let history = s.nonlinear.as_ref().ok_or_else(|| {
                    EngineError::invalid("frame", "linear static results have no stored frames; leave frame empty for the final state")
                })?;
                history.frames.get(k).ok_or_else(|| EngineError::invalid("frame", "no such nonlinear displacement frame"))?;
                // A repeated load factor does not identify committed material/contact state,
                // and stored displacement frames can be capped before the actual final state.
                if !displacement {
                    return Err(EngineError::invalid(
                        "field",
                        "nonlinear frames store displacements only; stresses, contact and temperature require the final state (frame null)",
                    ));
                }
            }
        }
        Kind::Heat(h) => {
            if mode.is_some() {
                return Err(EngineError::invalid("mode", "heat results have no modes"));
            }
            let last = h.frames.len().checked_sub(1).ok_or_else(|| EngineError::invalid("frame", "heat result has no stored temperatures"))?;
            let k = frame.unwrap_or(last);
            if k > last {
                return Err(EngineError::invalid("frame", format!("must be from 0 to {last}")));
            }
            if field == Field::HeatFlux && k != last {
                return Err(EngineError::invalid("field", "heat flux is stored only for the final state; select the final frame"));
            }
        }
        Kind::Explicit(e) => {
            if mode.is_some() {
                return Err(EngineError::invalid("mode", "explicit results have no modes"));
            }
            let k = frame.unwrap_or_else(|| e.frames.len().saturating_sub(1));
            if e.frames.get(k).is_none() {
                return Err(EngineError::invalid("frame", "no such explicit field frame"));
            }
        }
        Kind::Fluid(_) => return Err(EngineError::invalid("step", "path plots currently interpolate solid nodal fields; fluid fields are per cell")),
    }
    if !Field::available(result).any(|f| f == field) {
        return Err(EngineError::invalid("field", "field is not available for this solved step"));
    }
    Ok(pick)
}

fn path(s: &mut Session, p: &Value) -> Result<Value> {
    fieldcraft_model::progress::check()?;
    let step = step_for(s, p, Want::Any)?;
    let result =
        s.result(&step).cloned().ok_or_else(|| EngineError::invalid("step", "result is stale; solve the step again before sampling a path"))?;
    let field_id = p.get("field").filter(|v| !v.is_null()).ok_or_else(|| EngineError::param("field", crate::schema::ProblemKind::Required))?;
    let field = field_id.as_str().and_then(Field::parse).ok_or_else(|| EngineError::invalid("field", "use a supported nodal contour field id"))?;
    let pick = pick(&result, field, p)?;
    let points = super::points(p, "points", MAX_POINTS)?;
    let samples = super::opt::<usize>(p, "samples")?.unwrap_or(101);
    // Check limits/arclength before a potentially expensive index build or field allocation.
    crate::path::stations(&points, samples)?;
    if !(2..=MAX_SAMPLES).contains(&samples) {
        return Err(EngineError::invalid("samples", "must be from 2 to 10000"));
    }
    let path = p.get("path").filter(|v| !v.is_null()).map(|_| super::file_param(s, p, "path", fieldcraft_io::roots::Access::Write)).transpose()?;
    if path.as_ref().is_some_and(|p| p.ends_with(['/', '\\'])) {
        return Err(EngineError::invalid("path", "name a CSV file, not a directory path ending in a separator"));
    }
    let definition =
        s.result_inputs.get(&step).copied().ok_or_else(|| EngineError::invalid("step", "result has no definition identity; solve the step again"))?;
    let index = s.path_index.index(&s.model, &result, definition)?;
    let values = field_values(&s.model, &result, field, pick)?;
    let plot = Arc::new(index.sample(&s.model, &result, field, pick, points, samples, &values)?);
    let shown = |q, v| s.units().express(q, v).map_err(|e| EngineError::Params(e.to_string()));
    let distance = plot.distance.iter().map(|v| shown(Quantity::Length, *v)).collect::<Result<Vec<_>>>()?;
    let xyz = plot
        .xyz
        .iter()
        .map(|p| Ok([shown(Quantity::Length, p[0])?, shown(Quantity::Length, p[1])?, shown(Quantity::Length, p[2])?]))
        .collect::<Result<Vec<_>>>()?;
    let values = plot
        .values
        .iter()
        .map(|v| {
            v.map(|v| {
                let shown = field.shown(&result, v, s.units());
                if shown.is_finite() { Ok(shown) } else { Err(EngineError::invalid("field", "displayed field value overflows the selected units")) }
            })
            .transpose()
        })
        .collect::<Result<Vec<_>>>()?;
    let shown_mode = plot.mode;
    let mut reply = json!({"step":step,"field":field.id(),"frame":pick.frame,"mode":shown_mode,"distance":distance,"xyz":xyz,"values":values,"elements":plot.elements,
        "units":{"distance":s.units().symbol(Quantity::Length),"xyz":s.units().symbol(Quantity::Length),"value":field.unit(&result,s.units())},
        "field_title":field.title(&result,s.units()),"configuration":"reference","samples":samples});
    if let Some(path) = path {
        fieldcraft_model::progress::check()?;
        let rows = distance
            .iter()
            .zip(&xyz)
            .zip(&values)
            .map(|((d, p), value)| {
                vec![Cell::Number(*d), Cell::Number(p[0]), Cell::Number(p[1]), Cell::Number(p[2]), value.map_or(Cell::Empty, Cell::Number)]
            })
            .collect();
        let table = Table { header: ["distance", "x", "y", "z", "value"].map(str::to_owned).to_vec(), rows };
        fieldcraft_io::fs::write_output(std::path::Path::new(&path), table.csv().as_bytes())?;
        reply["path"] = json!(path);
    }
    // Failed writes/searches leave the preceding successful plot intact.
    s.path_plot = Some(plot);
    Ok(reply)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::units::{DisplaySystem, Preset};

    fn solved() -> Session {
        let mut s = Session::default();
        crate::cmd::solve::tests::cantilever(&mut s);
        s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
        s
    }

    fn request() -> Value {
        json!({"step":"Bend","field":"stress_xx","points":[[2.0,0.4,0.7],[8.0,0.4,0.7]],"samples":7,"frame":null,"mode":null})
    }

    fn linear_field(s: &mut Session) {
        let mut result = s.result("Bend").unwrap().as_ref().clone();
        let StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = &mut result else { panic!("static") };
        r.stress = r
            .nodes
            .iter()
            .map(|n| {
                let p = s.model.node(*n).unwrap().xyz;
                [7.0 + 2.0 * p[0] + 3.0 * p[1] - p[2], 0.0, 0.0, 0.0, 0.0, 0.0]
            })
            .collect();
        s.store_result(result);
    }

    #[test]
    fn public_query_interpolates_a_linear_field_reuses_index_and_is_not_an_edit() {
        let mut s = solved();
        linear_field(&mut s);
        let journal = s.journal.len();
        let revision = s.revision();
        let reply = s.execute("results.path", &request()).unwrap();
        assert_eq!(reply["samples"], 7);
        assert!(reply["mode"].is_null());
        assert_eq!(reply["units"]["value"], "MPa");
        for (point, value) in reply["xyz"].as_array().unwrap().iter().zip(reply["values"].as_array().unwrap()) {
            let shown: [f64; 3] = serde_json::from_value(point.clone()).unwrap();
            let si = shown.map(|v| s.units().to_si(Quantity::Length, v).unwrap());
            let expected = 7.0 + 2.0 * si[0] + 3.0 * si[1] - si[2];
            let actual = s.units().to_si(Quantity::Pressure, value.as_f64().unwrap()).unwrap();
            assert!((actual - expected).abs() < 1e-12 * expected.abs(), "{reply}");
        }
        let result = s.result("Bend").unwrap().clone();
        let definition = s.result_inputs["Bend"];
        let index = s.path_index.index(&s.model, &result, definition).unwrap();
        s.execute("results.path", &request()).unwrap();
        let reused = s.path_index.index(&s.model, &result, definition).unwrap();
        assert!(Arc::ptr_eq(&index, &reused));
        assert_eq!(s.journal.len(), journal);
        assert_eq!(s.revision(), revision);
        assert!(s.path_plot.as_ref().unwrap().current(&s));
        // A fresh solve changes exact result identity even at the same model definition.
        s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
        assert!(!s.path_plot.as_ref().unwrap().current(&s));
        let fresh = s.result("Bend").unwrap().clone();
        let rebuilt = s.path_index.index(&s.model, &fresh, definition).unwrap();
        assert!(!Arc::ptr_eq(&index, &rebuilt));
    }

    #[test]
    fn path_gaps_are_null_and_shared_faces_choose_lowest_label() {
        let mut s = solved();
        let reply =
            s.execute("results.path", &json!({"step":"Bend","field":"displacement_z","points":[[-1,0.5,0.5],[11,0.5,0.5]],"samples":13})).unwrap();
        let values = reply["values"].as_array().unwrap();
        assert!(values[0].is_null() && values[12].is_null());
        assert!(values[1..12].iter().all(Value::is_number), "{reply}");
        assert_eq!(reply["elements"][2], 1, "x=1 is shared by elements 1 and 2");
        // Delete an interior element from the model before storing the same synthetic result:
        // both sides retain nodes, but the empty volume must remain a gap.
        s.execute("mesh.element.remove", &json!({"label":5})).unwrap();
        s.execute("model.set.edit", &json!({"name":"all","elements":[1,2,3,4,6,7,8,9,10]})).unwrap();
        let result = s.latest_result("Bend").unwrap().0.as_ref().clone();
        s.store_result(result);
        let gap =
            s.execute("results.path", &json!({"step":"Bend","field":"displacement_z","points":[[4.25,0.5,0.5],[4.75,0.5,0.5]],"samples":3})).unwrap();
        assert!(gap["values"].as_array().unwrap().iter().all(Value::is_null));
        assert!(gap["elements"].as_array().unwrap().iter().all(Value::is_null));
    }

    #[test]
    fn public_points_distances_values_and_csv_convert_once_in_every_preset() {
        let mut s = solved();
        linear_field(&mut s);
        let si_points = [[0.002, 0.0004, 0.0007], [0.008, 0.0004, 0.0007]];
        for preset in Preset::ALL {
            let display = DisplaySystem::new(preset);
            s.default_units = display.clone();
            let points = si_points.map(|p| p.map(|v| display.express(Quantity::Length, v).unwrap()));
            let file = std::env::temp_dir().join(format!("fieldcraft-path-units-{:?}-{}.csv", preset, std::process::id()));
            let reply = s
                .execute("results.path", &json!({"step":"Bend","field":"stress_xx","points":points,"samples":7,"path":file.to_str().unwrap()}))
                .unwrap();
            let lines = std::fs::read_to_string(&file).unwrap();
            let mut rows = lines.lines();
            assert_eq!(rows.next(), Some("distance,x,y,z,value"));
            for (i, row) in rows.enumerate() {
                let cells: Vec<f64> = row.split(',').map(|v| v.parse().unwrap()).collect();
                assert_eq!(cells[0], reply["distance"][i].as_f64().unwrap());
                for k in 0..3 {
                    assert_eq!(cells[k + 1], reply["xyz"][i][k].as_f64().unwrap());
                }
                assert_eq!(cells[4], reply["values"][i].as_f64().unwrap());
                let x = si_points[0][0] + 0.001 * i as f64;
                assert!((display.to_si(Quantity::Length, cells[0]).unwrap() - 0.001 * i as f64).abs() < 1e-15);
                assert!((display.to_si(Quantity::Pressure, cells[4]).unwrap() - (7.0 + 2.0 * x + 0.0012 - 0.0007)).abs() < 1e-12);
            }
            std::fs::remove_file(file).unwrap();
        }
    }

    #[test]
    fn stale_hostile_and_failed_exports_preserve_previous_successful_plot() {
        let mut s = solved();
        s.execute("results.path", &request()).unwrap();
        let previous = s.path_plot.clone().unwrap();
        let mut hostile = vec![
            json!({}),
            json!({"field":"not_a_field"}),
            json!({"field":"contact_status"}),
            json!({"samples":1}),
            json!({"samples":10001}),
            json!({"samples":-1}),
            json!({"frame":0}),
            json!({"frame":-1}),
            json!({"mode":1}),
            json!({"points":[[1,2,3]]}),
            json!({"points":[[0,0,0],[0,0,0]]}),
            json!({"points":[[0,0],[1,0,0]]}),
            json!({"points":vec![[0,0,0];MAX_POINTS+1]}),
            json!({"extra":"rejected"}),
        ];
        let directory = std::env::temp_dir();
        hostile.push(json!({"path":directory.to_str().unwrap()}));
        let missing = directory.join(format!("fieldcraft-path-no-parent-{}", std::process::id())).join("out.csv");
        hostile.push(json!({"path":missing.to_str().unwrap()}));
        for change in hostile {
            let mut p = request();
            if change.as_object().unwrap().is_empty() {
                p.as_object_mut().unwrap().remove("field");
            } else {
                for (k, v) in change.as_object().unwrap() {
                    p[k] = v.clone();
                }
            }
            let tested = s.execute("results.path", &p);
            assert!(tested.is_err(), "{p}: {tested:?}");
            assert!(Arc::ptr_eq(s.path_plot.as_ref().unwrap(), &previous));
        }
        s.execute("model.material.create", &json!({"name":"Other","elastic":{"type":"isotropic","young":1200,"poisson":0.2},"density":1})).unwrap();
        // An unrelated material can keep a result current; changing a used modulus cannot.
        s.execute("model.material.edit", &json!({"name":"Mat","elastic":{"type":"isotropic","young":1200,"poisson":0.2},"density":1})).unwrap();
        assert!(s.execute("results.path", &request()).is_err());
    }

    #[test]
    fn intermediate_fields_are_refused_without_final_state_substitution() {
        let mut s = solved();
        let mut result = s.result("Bend").unwrap().as_ref().clone();
        let StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = &mut result else { panic!("static") };
        r.nonlinear = Some(fieldcraft_solver::nonlinear::History {
            frames: vec![fieldcraft_solver::nonlinear::Frame { load_factor: 0.0, displacement: vec![[0.0; 3]; r.nodes.len()], reaction: vec![] }],
            ..Default::default()
        });
        s.store_result(result);
        let mut p = request();
        p["frame"] = json!(0);
        assert!(s.execute("results.path", &p).unwrap_err().to_string().contains("displacements only"));
        p["field"] = json!("displacement_x");
        assert!(s.execute("results.path", &p).unwrap()["values"].as_array().unwrap().iter().all(|v| v.as_f64() == Some(0.0)));
        p["frame"] = json!(1);
        assert!(s.execute("results.path", &p).is_err());
        // A stored thermal frame may be sampled, but final-only heat flux may not.
        let nodes = s.result("Bend").unwrap().nodes().to_vec();
        let count = nodes.len();
        let heat = fieldcraft_solver::HeatResult {
            step: "Bend".into(),
            nodes,
            times: vec![0.0, 1.0],
            frames: vec![vec![300.0; count], vec![310.0; count]],
            flux: vec![[2.0, 0.0, 0.0]; count],
            elements: vec![],
            reaction: vec![0.0; count],
            loads: vec![],
            applied_total: 0.0,
            reaction_total: 0.0,
            stats: Default::default(),
            warnings: vec![],
        };
        s.store_result(StepResult::Implicit(fieldcraft_solver::StepResult::Heat(heat)));
        p["frame"] = json!(0);
        p["field"] = json!("heat_flux");
        assert!(s.execute("results.path", &p).is_err());
        p["field"] = json!("temperature");
        for preset in Preset::ALL {
            for temperature in ["K", "°C", "°F"] {
                s.default_units = DisplaySystem::new(preset).with_override(Quantity::Temperature, temperature).unwrap();
                let mut query = p.clone();
                query["points"] = json!([["2 mm", "0.4 mm", "0.7 mm"], ["8 mm", "0.4 mm", "0.7 mm"]]);
                let file = std::env::temp_dir().join(format!("fieldcraft-path-temperature-{}.csv", std::process::id()));
                query["path"] = json!(file.to_str().unwrap());
                let reply = s.execute("results.path", &query).unwrap();
                let value = reply["values"][0].as_f64().unwrap();
                assert!((s.units().to_si(Quantity::Temperature, value).unwrap() - 300.0).abs() < 1e-12);
                let csv = std::fs::read_to_string(&file).unwrap();
                let row = csv.lines().nth(1).unwrap().split(',').next_back().unwrap().parse::<f64>().unwrap();
                assert_eq!(row, value);
                std::fs::remove_file(file).unwrap();
            }
        }
    }

    #[test]
    fn actual_cantilever_top_fibre_tracks_euler_bernoulli_and_exact_nodes() {
        let mut s = Session::default();
        s.run_script(&crate::sample::cantilever()).unwrap();
        s.execute("solve.static", &json!({"step":"Load"})).unwrap();
        let reply = s.execute("results.path", &json!({"step":"Load","field":"stress_xx","points":[[20,5,20],[80,5,20]],"samples":101})).unwrap();
        for (point, value) in reply["xyz"].as_array().unwrap().iter().zip(reply["values"].as_array().unwrap()) {
            let x = point[0].as_f64().unwrap();
            let expected = 1.5 * (100.0 - x);
            assert!((value.as_f64().unwrap() - expected).abs() / expected < 0.10, "{point}: {value}, expected {expected}");
        }
        let result = s.result("Load").unwrap().clone();
        let values = field_values(&s.model, &result, Field::StressXx, Pick { mode: 1, frame: None }).unwrap();
        for (i, node) in result.nodes().iter().enumerate().step_by(11) {
            let p = s.model.node(*node).unwrap().xyz.map(|v| s.units().express(Quantity::Length, v).unwrap());
            let end = [p[0] + 1e-8, p[1], p[2]];
            let sampled = s.execute("results.path", &json!({"step":"Load","field":"stress_xx","points":[p,end],"samples":2})).unwrap();
            let actual = s.units().to_si(Quantity::Pressure, sampled["values"][0].as_f64().unwrap()).unwrap();
            assert!((actual - values[i]).abs() <= 1e-12 * values[i].abs().max(1.0));
        }
    }

    #[test]
    fn roots_and_gap_csv_preserve_blank_values() {
        let mut s = solved();
        let file = std::env::temp_dir().join(format!("fieldcraft-path-gaps-{}.csv", std::process::id()));
        let p = json!({"step":"Bend","field":"displacement_z","points":[[-1,0.5,0.5],[11,0.5,0.5]],"samples":13,"path":file.to_str().unwrap()});
        s.file_access = fieldcraft_io::roots::Roots::new(None, None).unwrap();
        assert!(s.execute("results.path", &p).is_err());
        s.file_access = fieldcraft_io::roots::Roots::writes_within(&std::env::temp_dir()).unwrap();
        let reply = s.execute("results.path", &p).unwrap();
        let text = std::fs::read_to_string(&file).unwrap();
        let rows: Vec<_> = text.lines().collect();
        assert!(rows[1].ends_with(',') && rows[13].ends_with(','));
        assert!(reply["values"][0].is_null());
        std::fs::remove_file(file).unwrap();
    }

    #[test]
    fn revisited_load_factor_is_not_the_final_committed_field_state() {
        let mut s = solved();
        let mut result = s.result("Bend").unwrap().as_ref().clone();
        let StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = &mut result else { panic!("static") };
        let count = r.nodes.len();
        r.temperature = Some(vec![300.0; count]);
        r.nonlinear = Some(fieldcraft_solver::nonlinear::History {
            frames: vec![
                fieldcraft_solver::nonlinear::Frame { load_factor: 1.0, displacement: vec![[0.0; 3]; count], reaction: vec![] },
                fieldcraft_solver::nonlinear::Frame { load_factor: 1.0, displacement: vec![[0.0001, 0.0, 0.0]; count], reaction: vec![] },
            ],
            increments: vec![fieldcraft_solver::nonlinear::Increment {
                number: 3,
                load_factor: 1.0,
                size: 1.0,
                iterations: vec![],
                converged: true,
                cutback: None,
                negative_pivots: 0,
                augmentations: 0,
            }],
            ..Default::default()
        });
        s.store_result(result);
        for frame in [0, 1] {
            for field in ["stress_xx", "temperature", "contact_pressure"] {
                let mut p = request();
                p["frame"] = json!(frame);
                p["field"] = json!(field);
                assert!(s.execute("results.path", &p).unwrap_err().to_string().contains("displacements only"));
            }
            let mut p = request();
            p["frame"] = json!(frame);
            p["field"] = json!("displacement_x");
            let reply = s.execute("results.path", &p).unwrap();
            assert!((reply["values"][0].as_f64().unwrap() - 0.1 * frame as f64).abs() < 1e-12);
        }
        let mut final_query = request();
        final_query["frame"] = Value::Null;
        assert!(s.execute("results.path", &final_query).is_ok());
    }

    #[test]
    fn cancellation_preserves_the_previous_plot_and_checks_index_and_sampling() {
        let mut s = solved();
        let progress = fieldcraft_model::progress::Progress::new();
        progress.cancel();
        let stopped = fieldcraft_model::progress::scoped(Some(progress.clone()), || s.execute("results.path", &request()));
        assert!(stopped.unwrap_err().is_cancelled());
        assert!(s.path_plot.is_none());
        s.execute("results.path", &request()).unwrap();
        let previous = s.path_plot.clone().unwrap();
        let result = s.result("Bend").unwrap().clone();
        let definition = s.result_inputs["Bend"];
        let index = s.path_index.index(&s.model, &result, definition).unwrap();
        let values = field_values(&s.model, &result, Field::StressXx, Pick { mode: 1, frame: None }).unwrap();
        let stopped = fieldcraft_model::progress::scoped(Some(progress.clone()), || {
            index.sample(
                &s.model,
                &result,
                Field::StressXx,
                Pick { mode: 1, frame: None },
                vec![[0.002, 0.0004, 0.0007], [0.008, 0.0004, 0.0007]],
                101,
                &values,
            )
        });
        assert!(stopped.unwrap_err().is_cancelled());
        let mut empty = crate::path::Cache::default();
        let stopped = fieldcraft_model::progress::scoped(Some(progress), || empty.index(&s.model, &result, definition));
        assert!(stopped.unwrap_err().is_cancelled());
        assert!(Arc::ptr_eq(s.path_plot.as_ref().unwrap(), &previous));
    }

    #[test]
    fn modal_shapes_use_the_existing_normalized_field_boundary_not_length_units() {
        let mut s = solved();
        s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap();
        let result = s.result("Modes").unwrap().clone();
        let Kind::Frequency(f) = result.kind() else { panic!("modes") };
        let label = f.nodes[10];
        let point = s.model.node(label).unwrap().xyz;
        let stored = f.modes[1].shape[10][2];
        for preset in Preset::ALL {
            s.default_units = DisplaySystem::new(preset);
            let p = point.map(|v| s.units().express(Quantity::Length, v).unwrap());
            let q = [p[0] * 0.999999, p[1], p[2]];
            let reply = s.execute("results.path", &json!({"step":"Modes","field":"displacement_z","mode":2,"points":[p,q],"samples":2})).unwrap();
            assert_eq!(reply["mode"], 2);
            assert!(reply["units"]["value"].is_null());
            assert!(reply["field_title"].as_str().unwrap().contains("normalised"));
            let expected = Field::DisplacementZ.shown(&result, stored, s.units());
            assert!((reply["values"][0].as_f64().unwrap() - expected).abs() <= 1e-12 * expected.abs().max(1.0));
            let factor = 2.0 * crate::results::deformation_unit(&result, s.units());
            let placed = s.path_plot.as_ref().unwrap().displayed_point(&s, 0, "Modes", Pick { mode: 2, frame: None }, factor).unwrap();
            let expected = [0, 1, 2].map(|k| point[k] + factor * f.modes[1].shape[10][k]);
            assert!((0..3).all(|k| (placed[k] - expected[k]).abs() <= 1e-12 * expected[k].abs().max(1.0)));
        }
        for extra in [json!({"mode":0}), json!({"mode":999}), json!({"frame":0})] {
            let mut p = json!({"step":"Modes","field":"displacement_z","points":[["2 mm","0.4 mm","0.7 mm"],["8 mm","0.4 mm","0.7 mm"]]});
            for (k, v) in extra.as_object().unwrap() {
                p[k] = v.clone();
            }
            assert!(s.execute("results.path", &p).is_err());
        }
    }

    #[test]
    fn stored_station_overlay_uses_reference_or_current_deformation_and_guards_stale_steps() {
        let mut s = solved();
        let result = s.result("Bend").unwrap().clone();
        let Kind::Static(r) = result.kind() else { panic!("static") };
        let node = r.nodes[10];
        let point = s.model.node(node).unwrap().xyz;
        let p = point.map(|v| s.units().express(Quantity::Length, v).unwrap());
        s.execute("results.path", &json!({"step":"Bend","field":"stress_xx","points":[p,[p[0]*0.99999,p[1],p[2]]],"samples":2})).unwrap();
        let path = s.path_plot.clone().unwrap();
        let pick = Pick { mode: 1, frame: None };
        assert_eq!(path.displayed_point(&s, 0, "Bend", pick, 0.0), Some(point));
        let factor = 1.7;
        let actual = path.displayed_point(&s, 0, "Bend", pick, factor).unwrap();
        let expected = [0, 1, 2].map(|k| point[k] + factor * r.displacement[10][k]);
        assert!((0..3).all(|k| (actual[k] - expected[k]).abs() < 1e-12));
        assert!(path.displayed_point(&s, 0, "Modes", pick, factor).is_none());
        assert!(path.displayed_point(&s, 100, "Bend", pick, factor).is_none());
        assert!(path.displayed_point(&s, 0, "Bend", pick, f64::NAN).is_none());
        let gap = s.execute("results.path", &json!({"step":"Bend","field":"stress_xx","points":[[-1,0.5,0.5],[-0.5,0.5,0.5]],"samples":2})).unwrap();
        assert!(gap["values"][0].is_null());
        let gap = s.path_plot.clone().unwrap();
        assert!(gap.displayed_point(&s, 0, "Bend", pick, 0.0).is_some());
        assert!(gap.displayed_point(&s, 0, "Bend", pick, 1.0).is_none());
        s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
        assert!(path.displayed_point(&s, 0, "Bend", pick, 0.0).is_none());
        let mut missing = s.result("Bend").unwrap().as_ref().clone();
        let StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = &mut missing else { panic!("static") };
        r.displacement[10] = [f64::NAN; 3];
        s.store_result(missing);
        s.execute("results.path", &json!({"step":"Bend","field":"stress_xx","points":[p,[p[0]*0.99999,p[1],p[2]]],"samples":2})).unwrap();
        let undefined = s.path_plot.as_ref().unwrap();
        assert!(undefined.displayed_point(&s, 0, "Bend", pick, 0.0).is_some());
        assert!(undefined.displayed_point(&s, 0, "Bend", pick, 1.0).is_none());
    }

    #[test]
    fn documented_conductor_recipe_matches_analytic_temperature_and_flux() {
        let script: Value = serde_json::from_str(include_str!("../../../../docs/recipes/path-temperature.json")).unwrap();
        let mut s = Session::default();
        s.run_script(&script).unwrap();
        let plot = s.path_plot.as_ref().unwrap();
        assert_eq!(s.model.elements().len(), 48);
        assert_eq!(plot.values.len(), 121);
        for (xyz, value) in plot.xyz.iter().zip(&plot.values) {
            let expected = 293.15 + xyz[2] * (10.0 / 0.12);
            assert!((value.unwrap() - expected).abs() < 1e-9, "{xyz:?}: {value:?}");
        }
        let result = s.result("Conduction").unwrap();
        let Kind::Heat(heat) = result.kind() else { panic!("heat result") };
        assert!(heat.flux.iter().all(|q| q[0].abs() < 1e-7 && q[1].abs() < 1e-7 && (q[2] + 1250.0).abs() < 1e-7));
        let reply =
            s.execute("results.path", &json!({"step":"Conduction","field":"heat_flux","points":[[10,15,0],[10,15,120]],"samples":121})).unwrap();
        for value in reply["values"].as_array().unwrap() {
            let si = s.units().to_si(Quantity::HeatFlux, value.as_f64().unwrap()).unwrap();
            assert!((si - 1250.0).abs() < 1e-7, "{si}");
        }
    }
}
