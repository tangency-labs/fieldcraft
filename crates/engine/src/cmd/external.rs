//! Load external visualization grids independently of the editable model. A PVD is committed
//! to the session only after every referenced frame passes the same checks as a single VTU.

use std::path::{Path, PathBuf};

use fieldcraft_io::external::{Frame, MAX_BYTES, MAX_FRAME_BYTES, Series};
use fieldcraft_io::roots::Access;
use fieldcraft_model::units::{Quantity, Unit};
use serde_json::{Value, json};

use crate::external_view::{Association, Component, ExternalView};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const IMPORT: CommandSpec = CommandSpec {
    id: "results.import.vtu",
    label: "Import VTU Results",
    params: r#"{"path":"solution.vtu" | "solution.pvtu" | "series.pvd","unit"?:"m"} — replace the imported results with a visualization grid or time series. VTU carries no units. Set the length unit of coordinates and displacement vectors (default: meters). The editable model and its solved steps are kept. Saving the project embeds the imported frames by default."#,
    form: Some(super::forms::IMPORT_VTU),
    kind: CommandKind::Session,
    requires: &[],
    run: import,
};

pub const VIEW: CommandSpec = CommandSpec {
    id: "results.external.view",
    label: "External Result View",
    params: r#"{"frame"?:0,"association"?:"point"|"cell","field"?:"temperature","component"?:"value"|"magnitude"|"x"|"y"|"z","displacement"?:"u","scale"?:1,"show"?:true} — select an imported frame, contour or point displacement. An empty field or displacement name clears that selection. Returns the selected frame's time and available arrays."#,
    form: Some(super::forms::EXTERNAL_VIEW),
    kind: CommandKind::Session,
    requires: &[Requirement::ExternalResults],
    run: view,
};

fn import(s: &mut Session, p: &Value) -> Result<Value> {
    let symbol = super::opt::<String>(p, "unit")?;
    let unit = symbol.as_deref().map(Unit::parse).transpose().map_err(|e| EngineError::invalid("unit", e.to_string()))?;
    if let Some(unit) = &unit {
        Quantity::Length.to_si(0.0, unit).map_err(|e| EngineError::invalid("unit", e.to_string()))?;
    }
    let unit = unit.as_ref().filter(|u| u.symbol() != "m");
    let path = super::file_param(s, p, "path", Access::Read)?;
    let source = Path::new(&path);
    let mut remaining_input = MAX_BYTES;
    let frames = if extension(source, "pvd") {
        read_collection(s, source, &mut remaining_input, unit)?
    } else if extension(source, "vtu") || extension(source, "pvtu") {
        vec![read_grid(s, source, extension(source, "pvtu"), 0.0, &mut remaining_input, MAX_FRAME_BYTES, unit)?]
    } else {
        return Err(EngineError::invalid("path", "expected a .vtu grid, a single-piece .pvtu, or a .pvd time series"));
    };
    let series = Series { frames };
    let initial = ExternalView::initial(&series);
    let result = describe(&series, &initial)?;
    s.set_external_results(series)?;
    Ok(result)
}

fn extension(path: &Path, expected: &str) -> bool {
    path.extension().is_some_and(|e| e.eq_ignore_ascii_case(expected))
}

fn utf8_path(path: &Path) -> Result<&str> {
    path.to_str().ok_or_else(|| EngineError::invalid("path", "the resolved VTK path is not valid UTF-8"))
}

fn read_bytes(path: &Path, limit: usize, remaining: &mut usize) -> Result<Vec<u8>> {
    let bytes = super::file::read_file(utf8_path(path)?, (*remaining).min(limit) as u64)?;
    *remaining = remaining.checked_sub(bytes.len()).ok_or_else(|| EngineError::invalid("path", "VTK input exceeds the series byte limit"))?;
    Ok(bytes)
}

fn relative_file(s: &Session, index: &Path, name: &str) -> Result<PathBuf> {
    let base = index.parent().unwrap_or(Path::new(""));
    let mut child = base.to_owned();
    child.extend(Path::new(name).components().filter(|c| *c != std::path::Component::CurDir));
    // A metadata file cannot grant additional filesystem access, including through symlinks.
    #[cfg(not(target_arch = "wasm32"))]
    let child = {
        let dir = if base.as_os_str().is_empty() { Path::new(".") } else { base };
        fieldcraft_io::roots::Roots::new(Some(dir), None)?.check(&child, Access::Read)?
    };
    Ok(s.checked_path(utf8_path(&child)?, Access::Read)?.into())
}

fn read_grid(
    s: &Session,
    source: &Path,
    is_wrapper: bool,
    time: f64,
    remaining_input: &mut usize,
    remaining_decoded: usize,
    unit: Option<&Unit>,
) -> Result<Frame> {
    let piece;
    let source = if is_wrapper {
        let xml = read_bytes(source, fieldcraft_io::vtu::MAX_INDEX_BYTES, remaining_input)?;
        let name = fieldcraft_io::vtu::single_piece(&xml)?;
        piece = relative_file(s, source, &name)?;
        piece.as_path()
    } else {
        source
    };
    let bytes = read_bytes(source, MAX_FRAME_BYTES, remaining_input)?;
    Ok(fieldcraft_io::vtu::read_with_unit(&bytes, time, remaining_decoded, unit)?)
}

fn read_collection(s: &Session, source: &Path, remaining_input: &mut usize, unit: Option<&Unit>) -> Result<Vec<Frame>> {
    let index = read_bytes(source, fieldcraft_io::vtu::MAX_INDEX_BYTES, remaining_input)?;
    let entries = fieldcraft_io::vtu::collection(&index)?;
    let mut frames = Vec::new();
    let mut decoded = 0usize;
    for (time, name) in entries {
        let path = relative_file(s, source, &name)?;
        let frame = read_grid(s, &path, extension(Path::new(&name), "pvtu"), time, remaining_input, MAX_BYTES - decoded, unit)
            .map_err(|e| EngineError::invalid("path", format!("{name} at time {time}: {e}")))?;
        decoded = decoded.saturating_add(frame.storage_bytes());
        if decoded > MAX_BYTES {
            return Err(EngineError::invalid("path", format!("imported series exceeds {MAX_BYTES} bytes")));
        }
        frames.push(frame);
    }
    Ok(frames)
}

fn view(s: &mut Session, p: &Value) -> Result<Value> {
    let series = s.external_results().ok_or_else(|| EngineError::Params("import a VTU or PVD result first".into()))?;
    let mut selected = s.external_view.clone().unwrap_or_else(|| ExternalView::initial(series));
    let frame = super::opt(p, "frame")?.unwrap_or(selected.frame);
    let data = series.frames.get(frame).ok_or_else(|| EngineError::invalid("frame", "outside the imported time series"))?;
    if p.get("frame").is_some_and(|f| !f.is_null()) {
        selected.animate = false;
    }
    let association = super::opt::<Association>(p, "association")?;
    let field = super::opt::<String>(p, "field")?;
    if let Some(association) = association {
        selected.association = association;
    }
    if selected.frame != frame {
        selected.frame = frame;
        if field.is_none() {
            selected.inherit_contour(data, association);
        }
        selected.inherit_displacement(data, p.get("scale").is_some_and(|v| !v.is_null()));
    }
    if let Some(field) = field {
        selected.field = (!field.is_empty()).then_some(field);
        let fields = selected.association.fields(data);
        selected.component = if selected.field.as_ref().and_then(|n| fields.get(n)).is_some_and(|f| f.components > 1) {
            Component::Magnitude
        } else {
            Component::Value
        };
    }
    if let Some(component) = super::opt(p, "component")? {
        if selected.field.is_none() && component != Component::Value {
            return Err(EngineError::invalid("component", "select a contour field first"));
        }
        selected.component = component;
    }
    if let Some(name) = super::opt::<String>(p, "displacement")? {
        selected.displacement = (!name.is_empty()).then_some(name);
    }
    selected.scale = super::opt(p, "scale")?.unwrap_or(selected.scale);
    selected.show = super::opt(p, "show")?.unwrap_or(selected.show);
    if !selected.show {
        selected.animate = false;
    }
    selected.validate(series).map_err(EngineError::Params)?;
    let result = describe(series, &selected)?;
    s.external_view = Some(selected);
    Ok(result)
}

fn describe(series: &Series, selected: &ExternalView) -> Result<Value> {
    let data = series.frames.get(selected.frame).ok_or_else(|| EngineError::Params("empty imported result series".into()))?;
    let fields = |values: &std::collections::BTreeMap<String, fieldcraft_io::external::Field>| {
        values.iter().map(|(name, field)| json!({"name":name,"components":field.components})).collect::<Vec<_>>()
    };
    Ok(json!({
        "frames":series.frames.len(), "frame":selected.frame, "time":data.time, "unit":data.displacement_unit.as_ref().map_or("m", Unit::symbol),
        "points":data.grid.points.len(), "cells":data.grid.cells.len(),
        "point_fields":fields(&data.point_data), "cell_fields":fields(&data.cell_data), "view":selected
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    // Two triangles share a geometric point with different temperatures. Each VTU point
    // remains a separate sample. Moving every point also exercises frame-specific grids.
    fn grid(offset: f64) -> String {
        format!(
            r#"<VTKFile type="UnstructuredGrid" byte_order="LittleEndian" version="1.0">
<UnstructuredGrid><Piece NumberOfPoints="6" NumberOfCells="2">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{offset} 0 0  1 0 0  0 1 0  {offset} 0 0  0 -1 0  1 0 0</DataArray></Points>
<Cells><DataArray type="Int32" Name="connectivity" format="ascii">0 1 2 3 4 5</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">3 6</DataArray><DataArray type="UInt8" Name="types" format="ascii">5 5</DataArray></Cells>
<PointData><DataArray type="Float64" Name="T" format="ascii">21 22 23 31 32 33</DataArray>
<DataArray type="Float64" Name="u" NumberOfComponents="2" format="ascii">0 0.1 0 0 0 0 0 0 0 0 0 0</DataArray></PointData>
<CellData><DataArray type="Int32" Name="region" format="ascii">3 8</DataArray></CellData>
</Piece></UnstructuredGrid></VTKFile>"#
        )
    }

    fn collection(second: &str) -> String {
        format!(
            r#"<VTKFile type="Collection" version="1.0" byte_order="LittleEndian"><Collection>
<DataSet timestep="0" file="start.vtu"/><DataSet timestep="0.2" file="{second}"/></Collection></VTKFile>"#
        )
    }

    #[test]
    fn source_length_unit_keeps_coordinates_and_displacement_in_si_after_reopening() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-vtu-length-unit-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("source.vtu");
        let project = dir.join("embedded.fieldcraft");
        let mut s = Session::default();
        let displayed = |s: &Session| {
            let (frame, view) = crate::external_view::visible(s).unwrap();
            let (_, mesh, _) = crate::external_view::scene(
                frame,
                view,
                crate::view::Style::Shaded,
                crate::view::Show::DEFAULT,
                &crate::view::Palette::LIGHT,
                None,
                None,
            );
            (frame.grid.points.clone(), mesh.nodes.iter().map(|n| n.p).collect::<Vec<_>>(), frame.point_data["T"].values.clone())
        };
        let meters = grid(0.0);
        std::fs::write(&path, &meters).unwrap();
        s.execute(IMPORT.id, &json!({"path":path})).unwrap();
        s.execute(VIEW.id, &json!({"displacement":"u"})).unwrap();
        let expected = displayed(&s);
        assert_eq!(expected.1[0].y, 0.1);
        let millimeters =
            meters.replace("1 0 0  0 1 0", "1000 0 0  0 1000 0").replace("0 -1 0  1 0 0", "0 -1000 0  1000 0 0").replace("0 0.1 0", "0 100 0");
        std::fs::write(&path, millimeters).unwrap();
        let reply = s.execute(IMPORT.id, &json!({"path":path,"unit":"mm"})).unwrap();
        assert_eq!(reply["unit"], "mm");
        s.execute(VIEW.id, &json!({"displacement":"u"})).unwrap();
        assert_eq!(displayed(&s), expected);
        assert!(s.execute(IMPORT.id, &json!({"path":path,"unit":"s"})).is_err());
        assert_eq!(displayed(&s), expected);
        s.execute("file.save_as", &json!({"path":project,"results":true})).unwrap();
        std::fs::remove_file(&path).unwrap();
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":project})).unwrap();
        assert_eq!(displayed(&reopened), expected);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn checked_paths_must_be_representable_without_replacement() {
        use std::os::unix::ffi::OsStrExt;
        let name = std::ffi::OsStr::from_bytes(b"temperature-\xff.vtu");
        assert!(utf8_path(Path::new(name)).unwrap_err().to_string().contains("not valid UTF-8"));
        assert_eq!(utf8_path(Path::new("temperature-20°C.vtu")).unwrap(), "temperature-20°C.vtu");
    }

    #[test]
    fn single_piece_wrappers_open_directly_and_within_a_collection() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-pvtu-{}", std::process::id()));
        let folder = dir.join("step");
        fs::create_dir_all(&folder).unwrap();
        let leaf = folder.join("rank0.vtu");
        let wrapper = folder.join("field.pvtu");
        let index = dir.join("series.pvd");
        let xml = r#"<VTKFile type="PUnstructuredGrid"><PUnstructuredGrid GhostLevel="0"><Piece Source="rank0.vtu"/></PUnstructuredGrid></VTKFile>"#;
        fs::write(&leaf, grid(0.1)).unwrap();
        fs::write(&wrapper, xml).unwrap();
        fs::write(dir.join("start.vtu"), grid(0.)).unwrap();
        fs::write(&index, collection("step/field.pvtu")).unwrap();
        let mut s = Session::default();
        assert_eq!(s.execute(IMPORT.id, &json!({"path":wrapper})).unwrap()["frames"], 1);
        assert_eq!(s.execute(IMPORT.id, &json!({"path":index})).unwrap()["frames"], 2);
        assert_eq!(s.external_results().unwrap().frames[1].time, 0.2);
        let exact = fs::metadata(&leaf).unwrap().len() as usize + xml.len();
        let mut remaining = exact - 1;
        assert!(read_grid(&s, &wrapper, true, 0., &mut remaining, MAX_FRAME_BYTES, None).is_err());
        let second = xml.replace("</PUnstructuredGrid>", r#"<Piece Source="rank1.vtu"/></PUnstructuredGrid>"#);
        fs::write(&wrapper, second).unwrap();
        assert!(s.execute(IMPORT.id, &json!({"path":index})).unwrap_err().to_string().contains("multi-piece"));
        assert_eq!(s.external_results().unwrap().frames.len(), 2);
        #[cfg(unix)]
        {
            // A symlink that stays in the collection but escapes the wrapper's directory.
            let escape = folder.join("linked.vtu");
            std::os::unix::fs::symlink(dir.join("start.vtu"), &escape).unwrap();
            fs::write(&wrapper, xml.replace("rank0.vtu", "linked.vtu")).unwrap();
            assert!(s.execute(IMPORT.id, &json!({"path":index})).is_err());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn frame_changes_inherit_contours_and_deformation_independently() {
        use fieldcraft_io::external::{Cell, Field, Grid};
        let mut frames = [0., 1., 2., 3.].map(|time| Frame {
            displacement_unit: None,
            time,
            grid: Grid { points: vec![[1., 0., 0.]], cells: vec![Cell { vtk_type: 1, points: vec![0] }] },
            point_data: [
                ("A".into(), Field { components: 1, values: vec![1.] }),
                ("T".into(), Field { components: 3, values: vec![1., 2., 3.] }),
                ("u".into(), Field { components: 2, values: vec![0., 0.] }),
            ]
            .into(),
            cell_data: Default::default(),
        });
        frames[1].point_data.remove("u");
        frames[2].point_data.get_mut("u").unwrap().values[0] = fieldcraft_model::MAX_COORD;
        frames[3].point_data.remove("T");
        let mut s = Session::default();
        s.set_external_results(Series { frames: frames.into() }).unwrap();
        let selection = json!({"frame":0,"field":"T","component":"y","displacement":"u","scale":2.0});
        for request in [json!({"frame":1}), json!({"frame":1,"displacement":""}), json!({"frame":2})] {
            s.execute(VIEW.id, &selection).unwrap();
            let view = &s.execute(VIEW.id, &request).unwrap()["view"];
            assert_eq!((&view["field"], &view["component"], &view["displacement"]), (&json!("T"), &json!("y"), &Value::Null));
        }
        s.execute(VIEW.id, &selection).unwrap();
        for invalid in [
            json!({"frame":1,"displacement":"u"}),
            json!({"frame":1,"field":"missing"}),
            json!({"frame":3,"component":"z"}),
            json!({"frame":2,"scale":2.0}),
            json!({"frame":2,"displacement":"u"}),
        ] {
            assert!(s.execute(VIEW.id, &invalid).is_err(), "{invalid}");
            assert_eq!(s.external_view.as_ref().unwrap().frame, 0);
        }
        let view = &s.execute(VIEW.id, &json!({"frame":2,"scale":0.25})).unwrap()["view"];
        assert_eq!(view["displacement"], "u");
        let view = &s.execute(VIEW.id, &json!({"frame":3})).unwrap()["view"];
        assert_eq!(view["field"], "A");
        assert_eq!(view["displacement"], "u");
    }

    #[test]
    fn frame_and_association_use_the_inherited_name_at_its_requested_location() {
        use fieldcraft_io::external::{Cell, Field, Grid};
        let mut frames = [0., 1.].map(|time| Frame {
            displacement_unit: None,
            time,
            grid: Grid { points: vec![[0.; 3]], cells: vec![Cell { vtk_type: 1, points: vec![0] }] },
            point_data: [("T".into(), Field { components: 1, values: vec![0.] })].into(),
            cell_data: Default::default(),
        });
        let scalar = frames[1].point_data.remove("T").unwrap();
        frames[1].point_data.insert("A".into(), scalar.clone());
        frames[1].cell_data.insert("T".into(), scalar);
        let mut session = Session::default();
        session.set_external_results(Series { frames: frames.into() }).unwrap();
        let original = session.external_view.clone();
        for invalid in [json!({"frame":1,"association":"cell","field":"A"}), json!({"frame":1,"association":"cell","component":"y"})] {
            assert!(session.execute(VIEW.id, &invalid).is_err());
            assert_eq!(session.external_view, original);
        }
        let view = &session.execute(VIEW.id, &json!({"frame":1,"association":"cell"})).unwrap()["view"];
        assert_eq!(view["field"], "T");
        assert_eq!(view["association"], "cell");
        assert_eq!(view["component"], "value");
        // Frame zero has no cell arrays. An explicit cell location must not select a point array.
        let view = &session.execute(VIEW.id, &json!({"frame":0,"association":"cell"})).unwrap()["view"];
        assert_eq!(view["association"], "cell");
        assert!(view["field"].is_null());
    }

    #[test]
    fn moving_series_survives_save_without_sources_and_failed_replacement() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-external-series-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (pvd, first, last, project) = (dir.join("series.pvd"), dir.join("start.vtu"), dir.join("end.vtu"), dir.join("model.fieldcraft"));
        fs::write(&first, grid(0.0)).unwrap();
        fs::write(&last, grid(0.1)).unwrap();
        fs::write(&pvd, collection("end.vtu")).unwrap();
        let mut s = Session::default();
        s.execute("mesh.node.add", &json!({"label":17,"xyz":[1,2,3]})).unwrap();
        s.execute(IMPORT.id, &json!({"path":pvd})).unwrap();
        let select = json!({"frame":1,"association":"cell","field":"region","displacement":"u","scale":5.0});
        let r = s.execute(VIEW.id, &select).unwrap();
        assert_eq!((r["time"].as_f64(), r["points"].as_u64(), r["cells"].as_u64()), (Some(0.2), Some(6), Some(2)));
        assert_eq!(r["view"]["scale"], 5.0);
        let original = s.external_results().unwrap().clone();
        assert_eq!(original.frames[1].grid.points[0], [0.1, 0.0, 0.0]);
        assert_eq!(original.frames[1].point_data["T"].values[3], 31.0);
        fs::write(&pvd, collection("missing.vtu")).unwrap();
        assert!(s.execute(IMPORT.id, &json!({"path":pvd})).is_err());
        assert!(s.execute(VIEW.id, &json!({"frame":0,"displacement":"T"})).is_err());
        assert_eq!(s.external_results(), Some(&original));
        assert_eq!(s.external_view.as_ref().unwrap().frame, 1);
        s.execute("file.save_as", &json!({"path":project})).unwrap();
        for path in [&pvd, &first, &last] {
            fs::remove_file(path).unwrap();
        }
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":project})).unwrap();
        assert_eq!(reopened.external_results(), Some(&original));
        assert!(reopened.model.node(17).is_some());
        assert_eq!(reopened.execute(VIEW.id, &json!({})).unwrap(), r);
        fs::remove_dir_all(dir).unwrap();
    }
}
