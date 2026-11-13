//! OpenFOAM: write a fluid step as a case directory (`deck.openfoam`), and read the case's
//! results back onto the step's lattice (`results.openfoam`), so FieldCraft's fluid results
//! (slices, probes, lines, measurements) show what OpenFOAM computed.
//!
//! The case is written from OpenFOAM's public user guide; OpenFOAM runs as its own program.

use std::path::{Path, PathBuf};

use fieldcraft_decks::openfoam::{CaseOptions, FieldValues, Gas, ObstacleSurface, PressureBox, RefineBox, read_field, write_case};
use fieldcraft_io::roots::Access;
use fieldcraft_lbm::solve::{FluidResult, FluidStats, Frame};
use fieldcraft_lbm::voxel::Grid;
use fieldcraft_model::units::Quantity;
use fieldcraft_model::{BoxFace, FaceCondition, FluidDomain};
use serde_json::{Value, json};

use super::opt;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};

pub const WRITE: CommandSpec = CommandSpec {
    id: "deck.openfoam",
    label: "Write OpenFOAM Case",
    params: r#"{"path": "case_dir", "step"?: "Flow", "cell_factor"?: 2, "processes"?: 1, "max_courant"?: 0.8, "refine_boxes"?: [{"min": [x,y,z], "max": [x,y,z], "spacing"?: 0.1}], "location"?: [x,y,z], "initial_pressure"?: [{"min": [x,y,z], "max": [x,y,z], "pressure": 29.9}], "compressible"?: true | {"ambient"?: 0.101325, "temperature"?: 20, "molar_mass"?: 0.0289647, "specific_heat"?: 1004.5, "prandtl"?: 0.71}} — writes a fluid step (the only one when `step` is left out) as an OpenFOAM case directory: blockMesh background of spacing × cell_factor (1, 2, 4, 8), snappyHexMesh around the obstacles (one patch per named face, refined back to the spacing; `refine_boxes` to theirs, rounded to the background over a power of two), ASCII fields, decomposeParDict for `processes` ranks and an Allrun script (mesh, cell centres, solve, reconstruct). Incompressible by default: pimpleFoam with Smagorinsky LES, an adaptive time step (max_courant, 0.8; compressible 0.2). `compressible`: rhoCentralFoam, a laminar ideal gas (air by default): pressure faces above 0 are total-pressure inlets at `ambient` + their gauge pressure and the total `temperature`, the others non-reflecting outlets at `ambient`; walls adiabatic; `specific_heat` from the material's thermal data when left out, else air's; molar mass in kg/mol. Lengths, pressures and temperatures in the input units. `location` is a point in the fluid (found when left out). `initial_pressure` (compressible): boxes that start at that gauge pressure (setFields), such as the bore behind a high-pressure inlet, which would otherwise start as a shock at the inlet face. Fails naming what OpenFOAM cannot be given yet (heat, porous zones, fans, baffles, free surface, periodic faces, parabolic inlets, non-Newtonian fluids). Writes SI; nothing of the model changes; `results.openfoam` reads the results back"#,
    form: None,
    kind: CommandKind::Session,
    requires: &[Requirement::FluidDomain, Requirement::Step],
    run: write,
};

pub const READ: CommandSpec = CommandSpec {
    id: "results.openfoam",
    label: "Read OpenFOAM Results",
    params: r#"{"path": "case_dir", "spacing"?: 0.4, "max_cells"?: 20000000} — reads the results of a case `deck.openfoam` wrote and OpenFOAM ran (reconstructed, ASCII) as the fluid results of its step: every written time a frame, sampled at the centres of a lattice over the step's domain box (its spacing, or `spacing` in the input units) from the nearest OpenFOAM cell; the obstacles' cells are solid. Velocity, gauge pressure (compressible: above the ambient; incompressible: the kinematic pressure times the density) and, compressible, temperature. The reply states the frames, cells, peak speed and (compressible) peak Mach number. Replaces the step's results; the model is unchanged"#,
    form: None,
    kind: CommandKind::Session,
    requires: &[Requirement::FluidDomain, Requirement::Step],
    run: read,
};

/// The Courant number a compressible case adapts its time step to unless told otherwise.
const COMPRESSIBLE_COURANT: f64 = 0.2;

/// What FieldCraft wrote beside the case, for reading it back.
const MARKER: &str = "fieldcraft.json";

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct BoxParam {
    min: [f64; 3],
    max: [f64; 3],
    #[serde(default)]
    spacing: Option<f64>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct PressureParam {
    min: [f64; 3],
    max: [f64; 3],
    pressure: f64,
}

#[derive(serde::Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct GasParam {
    ambient: Option<f64>,
    temperature: Option<f64>,
    molar_mass: Option<f64>,
    specific_heat: Option<f64>,
    prandtl: Option<f64>,
}

fn si(s: &Session, key: &str, q: Quantity, v: f64) -> Result<f64> {
    s.input_units().to_si(q, v).map_err(|e| EngineError::invalid(key, e.to_string()))
}

/// The step `deck.openfoam` / `results.openfoam` act on: named, or the only fluid step.
fn fluid_step<'s>(s: &'s Session, name: Option<&str>) -> Result<(&'s str, &'s FluidDomain)> {
    let (step, domain) = fieldcraft_lbm::fluid_step(&s.model, name)?;
    Ok((step.name.as_str(), domain))
}

fn write(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = super::file_param(s, p, "path", Access::Write)?;
    let step: Option<String> = opt(p, "step")?;
    let length = |key: &str, v: [f64; 3]| -> Result<[f64; 3]> {
        let mut out = [0.0; 3];
        for (o, x) in out.iter_mut().zip(v) {
            *o = si(s, key, Quantity::Length, x)?;
        }
        Ok(out)
    };
    let defaults = CaseOptions::default();
    let mut options = CaseOptions {
        cell_factor: opt(p, "cell_factor")?.unwrap_or(defaults.cell_factor),
        processes: opt(p, "processes")?.unwrap_or(defaults.processes),
        max_courant: opt(p, "max_courant")?.unwrap_or(defaults.max_courant),
        ..defaults
    };
    for b in opt::<Vec<BoxParam>>(p, "refine_boxes")?.unwrap_or_default() {
        let spacing = b.spacing.map(|h| si(s, "refine_boxes", Quantity::Length, h)).transpose()?;
        options.refine_boxes.push(RefineBox { min: length("refine_boxes", b.min)?, max: length("refine_boxes", b.max)?, spacing });
    }
    for b in opt::<Vec<PressureParam>>(p, "initial_pressure")?.unwrap_or_default() {
        options.initial_pressure.push(PressureBox {
            min: length("initial_pressure", b.min)?,
            max: length("initial_pressure", b.max)?,
            pressure: si(s, "initial_pressure", Quantity::Pressure, b.pressure)?,
        });
    }
    if let Some(l) = opt::<[f64; 3]>(p, "location")? {
        options.location = Some(length("location", l)?);
    }
    let (step_name, domain) = fluid_step(s, step.as_deref())?;
    let gas = match p.get("compressible") {
        None | Some(Value::Bool(false)) | Some(Value::Null) => None,
        Some(Value::Bool(true)) => Some(GasParam::default()),
        Some(_) => opt::<GasParam>(p, "compressible")?,
    };
    if let Some(g) = gas {
        let air = Gas::AIR;
        let material_cp = s.model.materials().get(&domain.material).and_then(|m| m.thermal.as_ref()).and_then(|t| t.specific_heat);
        options.gas = Some(Gas {
            ambient: g.ambient.map(|v| si(s, "compressible.ambient", Quantity::Pressure, v)).transpose()?.unwrap_or(air.ambient),
            temperature: g.temperature.map(|v| si(s, "compressible.temperature", Quantity::Temperature, v)).transpose()?.unwrap_or(air.temperature),
            molar_mass: g.molar_mass.unwrap_or(air.molar_mass),
            specific_heat: g
                .specific_heat
                .map(|v| si(s, "compressible.specific_heat", Quantity::SpecificHeat, v))
                .transpose()?
                .or(material_cp)
                .unwrap_or(air.specific_heat),
            prandtl: g.prandtl.unwrap_or(air.prandtl),
        });
        // The central scheme through strong expansions wants a smaller Courant number.
        if p.get("max_courant").is_none() {
            options.max_courant = COMPRESSIBLE_COURANT;
        }
    }

    // The obstacles, tessellated as the lattice solver does, in metres.
    let tolerance = super::solve::obstacle_tolerance(domain);
    let mut meshes = Vec::new();
    for body in &domain.obstacles {
        let (mesh, faces) = super::solve::obstacle_mesh(s, body, tolerance)?;
        let positions: Vec<[f64; 3]> = mesh.positions.iter().map(|q| crate::units::geometry_point(*q)).collect();
        meshes.push((body.as_str(), positions, mesh.into_owned(), faces));
    }
    let surfaces: Vec<ObstacleSurface<'_>> = meshes
        .iter()
        .map(|(body, positions, mesh, faces)| ObstacleSurface { body, positions, triangles: &mesh.triangles, tri_face: &mesh.tri_face, faces })
        .collect();
    let case = write_case(&s.model, Some(step_name), &surfaces, &options)?;
    let density = s.model.materials().get(&domain.material).and_then(|m| m.density);
    let marker = json!({"step": step_name, "solver": case.solver, "density": density,
        "ambient": options.gas.map(|g| g.ambient), "temperature": options.gas.map(|g| g.temperature),
        "gamma": options.gas.map(|g| g.gamma()), "gas_constant": options.gas.map(|g| fieldcraft_decks::openfoam::GAS_CONSTANT / g.molar_mass),
        "specific_heat": options.gas.map(|g| g.specific_heat)});

    let root = Path::new(&dir);
    let mut bytes = 0;
    let marker_text = serde_json::to_string_pretty(&marker).map_err(|e| EngineError::Params(e.to_string()))?;
    for (rel, text) in case.files.iter().map(|(k, v)| (k.as_str(), v.as_str())).chain([(MARKER, marker_text.as_str())]) {
        let path = root.join(rel);
        if let Some(parent) = path.parent() {
            fieldcraft_io::fs::create_dir_all(parent)?;
        }
        fieldcraft_io::fs::write_output(&path, text.as_bytes())?;
        bytes += text.len();
    }
    Ok(json!({"format": "openfoam", "path": dir, "step": step_name, "solver": case.solver, "files": case.files.keys().collect::<Vec<_>>(),
        "bytes": bytes, "background_cells": case.background_cells, "location": case.location, "patches": case.patches, "notes": case.notes,
        "run": format!("cd {dir} && sh Allrun"), "read_back": {"command": "results.openfoam", "params": {"path": dir}}}))
}

/// The most bytes one field file may have.
const MAX_FIELD_BYTES: u64 = 4 << 30;
/// The most bytes the frames read back may take.
const MAX_FRAME_BYTES: u64 = 8 << 30;

#[derive(serde::Deserialize)]
struct Marker {
    step: String,
    solver: String,
    density: Option<f64>,
    ambient: Option<f64>,
    temperature: Option<f64>,
    gamma: Option<f64>,
    gas_constant: Option<f64>,
    #[serde(default)]
    specific_heat: Option<f64>,
}

fn field(path: &Path) -> Result<FieldValues> {
    let bytes = super::file::read_file(&path.to_string_lossy(), MAX_FIELD_BYTES)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::Params(format!("{} is not text", path.display())))?;
    read_field(text).map_err(|e| EngineError::context(format!("{}", path.display()), e))
}

/// The written times of the case, in order (`0`, the initial state, left out).
fn times(root: &Path) -> Vec<(f64, PathBuf)> {
    let mut out: Vec<(f64, PathBuf)> = fieldcraft_io::fs::list_dirs(root)
        .into_iter()
        .filter_map(|d| {
            let t: f64 = d.file_name()?.to_str()?.parse().ok()?;
            (t.is_finite() && t > 0.0).then_some((t, d))
        })
        .collect();
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

/// For each lattice cell, the OpenFOAM cell whose centre is nearest (`None` in a solid cell):
/// the centres are bucketed by lattice cell and searched ring by ring.
fn nearest_cells(grid: &Grid, solid: &[bool], centres: &[[f64; 3]]) -> Result<Vec<Option<u32>>> {
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); grid.len()];
    for (i, c) in centres.iter().enumerate() {
        if let Some(cell) = grid.cell_at(*c) {
            let i = u32::try_from(i).map_err(|_| EngineError::Params("more than 4 billion OpenFOAM cells".into()))?;
            if let Some(b) = buckets.get_mut(cell) {
                b.push(i);
            }
        }
    }
    let n = grid.cells;
    let reach = n.iter().copied().max().unwrap_or(0);
    let mut out = vec![None; grid.len()];
    for (cell, slot) in out.iter_mut().enumerate() {
        if solid.get(cell).copied().unwrap_or(true) {
            continue;
        }
        let x = grid.centre(cell);
        let ijk = grid.coords(cell);
        let mut best: Option<(f64, u32)> = None;
        for r in 0..=reach {
            // Every centre in ring r is at least (r − 1) spacings away.
            if let Some((d, _)) = best
                && d.sqrt() < (r as f64 - 1.0) * grid.spacing
            {
                break;
            }
            let lo = ijk.map(|v| v.saturating_sub(r));
            let hi = [0, 1, 2].map(|a| (ijk[a] + r).min(n[a].saturating_sub(1)));
            for k in lo[2]..=hi[2] {
                for j in lo[1]..=hi[1] {
                    for i in lo[0]..=hi[0] {
                        let on_ring = [i, j, k].iter().zip(ijk).any(|(v, c)| v.abs_diff(c) == r);
                        if !on_ring {
                            continue;
                        }
                        for &id in buckets.get(grid.index([i, j, k])).map(Vec::as_slice).unwrap_or_default() {
                            let c = centres.get(id as usize).copied().unwrap_or([f64::INFINITY; 3]);
                            let d = (c[0] - x[0]).powi(2) + (c[1] - x[1]).powi(2) + (c[2] - x[2]).powi(2);
                            if best.is_none_or(|(bd, bid)| d < bd || (d == bd && id < bid)) {
                                best = Some((d, id));
                            }
                        }
                    }
                }
            }
        }
        *slot = best.map(|(_, id)| id);
    }
    Ok(out)
}

fn read(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = super::file_param(s, p, "path", Access::Read)?;
    let root = PathBuf::from(&dir);
    let marker_path = root.join(MARKER);
    let marker_bytes = super::file::read_file(&marker_path.to_string_lossy(), 1 << 20)
        .map_err(|e| EngineError::context(format!("{} is not a case deck.openfoam wrote (no {MARKER})", root.display()), e))?;
    let marker: Marker = serde_json::from_slice(&marker_bytes).map_err(|e| EngineError::invalid("path", format!("{MARKER}: {e}")))?;
    let (step_name, domain) = fluid_step(s, Some(&marker.step))?;
    let step_name = step_name.to_string();
    if !fieldcraft_io::fs::list_dirs(&root).iter().all(|d| !d.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.starts_with("processor")))
        && times(&root).is_empty()
    {
        return Err(EngineError::invalid("path", "the case is decomposed and not reconstructed: run reconstructPar in it first"));
    }
    let compressible = marker.solver == "rhoCentralFoam";

    // The lattice over the domain box.
    let spacing = match opt::<f64>(p, "spacing")? {
        Some(h) => si(s, "spacing", Quantity::Length, h)?,
        None => domain.spacing,
    };
    if !(spacing.is_finite() && spacing > 0.0) {
        return Err(EngineError::invalid("spacing", "must be a positive length"));
    }
    let max_cells: u64 = opt(p, "max_cells")?.unwrap_or(20_000_000);
    let mut cells = [0usize; 3];
    for (a, c) in cells.iter_mut().enumerate() {
        let n = ((domain.max[a] - domain.min[a]) / spacing).round();
        if !(1.0..=1e6).contains(&n) {
            return Err(EngineError::invalid("spacing", "gives a lattice of 1 to a million cells along each axis"));
        }
        *c = n as usize;
    }
    let total = cells.iter().map(|c| *c as u64).product::<u64>();
    if total > max_cells {
        return Err(EngineError::invalid(
            "spacing",
            format!("gives {total} lattice cells, above max_cells {max_cells}: coarsen it or raise max_cells"),
        ));
    }
    let grid = Grid { origin: domain.min, spacing, cells };
    let mut solid = vec![false; grid.len()];
    for o in super::solve::obstacles(s, domain)? {
        fieldcraft_lbm::voxel::voxelise(&o.surface, &grid, &mut solid)?;
    }

    // The OpenFOAM cells' centres and the nearest one to each lattice cell.
    let centres = match field(&root.join("0/C"))? {
        FieldValues::Vectors(v) => v,
        _ => {
            return Err(EngineError::invalid(
                "path",
                "0/C holds no cell centres (Allrun writes them with `postProcess -func writeCellCentres -time 0`)",
            ));
        }
    };
    let nearest = nearest_cells(&grid, &solid, &centres)?;
    let written = times(&root);
    if written.is_empty() {
        return Err(EngineError::invalid("path", "the case has no written times yet: run it (sh Allrun)"));
    }
    let frame_bytes = grid.len() as u64 * if compressible { 20 } else { 16 };
    if frame_bytes.saturating_mul(written.len() as u64) > MAX_FRAME_BYTES {
        return Err(EngineError::invalid(
            "spacing",
            format!("{} frames of {} cells take more than {} GB: coarsen it", written.len(), grid.len(), MAX_FRAME_BYTES >> 30),
        ));
    }

    let reference = marker.temperature.unwrap_or(293.15);
    let ambient = marker.ambient.unwrap_or(0.0);
    let density = marker.density.unwrap_or(1.0);
    let gas_constant = marker.gas_constant.unwrap_or(287.0);
    let gamma = marker.gamma.unwrap_or(1.4);
    let volume = spacing.powi(3);
    let mut frames = Vec::with_capacity(written.len());
    let (mut max_speed, mut max_mach) = (0.0f64, 0.0f64);
    for (k, (time, path)) in written.iter().enumerate() {
        let u = field(&path.join("U"))?;
        let pf = field(&path.join("p"))?;
        let tf = if compressible { Some(field(&path.join("T"))?) } else { None };
        for f in [Some(&u), Some(&pf), tf.as_ref()].into_iter().flatten() {
            if f.count().is_some_and(|n| n != centres.len()) {
                return Err(EngineError::invalid("path", format!("{}: a field of another mesh than 0/C", path.display())));
            }
        }
        let mut frame = Frame { time: *time, step: k as u64, ..Frame::default() };
        frame.velocity = Vec::with_capacity(grid.len());
        frame.pressure = Vec::with_capacity(grid.len());
        let mut mass = 0.0;
        for slot in &nearest {
            let Some(id) = slot.map(|i| i as usize) else {
                frame.velocity.push([0.0; 3]);
                frame.pressure.push(0.0);
                if compressible {
                    frame.temperature.push(0.0);
                    frame.density.push(0.0);
                }
                continue;
            };
            let v = u.vector(id).unwrap_or([0.0; 3]);
            let pv = pf.scalar(id).unwrap_or(0.0);
            let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            max_speed = max_speed.max(speed);
            frame.velocity.push(v.map(|c| c as f32));
            if let Some(t) = tf.as_ref().and_then(|t| t.scalar(id)) {
                // The ideal gas: ρ = p / (R T).
                let rho = if t > 0.0 && pv > 0.0 { pv / (gas_constant * t) } else { 0.0 };
                frame.pressure.push((pv - ambient) as f32);
                frame.temperature.push((t - reference) as f32);
                frame.density.push(rho as f32);
                if t > 0.0 {
                    max_mach = max_mach.max(speed / (gamma * gas_constant * t).sqrt());
                }
                mass += rho * volume;
            } else {
                frame.pressure.push((pv * density) as f32);
                mass += density * volume;
            }
        }
        frame.mass = mass;
        frames.push(frame);
    }

    let fluid_cells = nearest.iter().filter(|c| c.is_some()).count();
    let open = |f: BoxFace| matches!(domain.faces.get(f), FaceCondition::Velocity { .. } | FaceCondition::Pressure { .. });
    let end_time = written.last().map_or(0.0, |w| w.0);
    let stats = FluidStats {
        validity: None,
        cells: grid.len(),
        fluid_cells,
        steps: 0,
        time_step: 0.0,
        end_time,
        tau_plus: 0.0,
        tau_minus: 0.0,
        magic: 0.0,
        collision: fieldcraft_model::Collision::default(),
        shear_dependent_tau: None,
        velocity: max_speed,
        lattice_velocity: 0.0,
        mach: max_mach,
        cell_updates: 0,
        backend: format!("openfoam ({})", marker.solver),
        tau_thermal: None,
        eddy_viscosity_ratio: None,
        frozen: None,
        free_surface: None,
        wall_links: None,
    };
    let history_times = frames.iter().map(|f| f.time).collect();
    let n_frames = frames.len();
    // Mass flows: the reference density (incompressible) or each cell's (compressible, its
    // frames' density); heat needs the specific heat.
    let port_properties = if compressible {
        let rho0 = ambient / (gas_constant * reference);
        (rho0.is_finite() && rho0 > 0.0).then_some(fieldcraft_lbm::PortProperties { density: rho0, specific_heat: marker.specific_heat })
    } else {
        marker.density.filter(|d| d.is_finite() && *d > 0.0).map(|density| fieldcraft_lbm::PortProperties { density, specific_heat: None })
    };
    let result = FluidResult {
        port_properties,
        step: step_name.clone(),
        grid,
        solid,
        open_faces: BoxFace::ALL.map(open),
        frames,
        history_times,
        forces: Vec::new(),
        heat: Vec::new(),
        baffles: Box::default(),
        zones: Vec::new(),
        stats,
        warnings: vec![format!("read from the OpenFOAM case {} ({}), not solved on the lattice", root.display(), marker.solver)],
        steady: None,
        temperature_reference: reference,
        solids: Vec::new(),
        surface: None,
    };
    s.store_result(StepResult::Fluid(Box::new(result)));
    let mut reply = json!({"step": step_name, "path": dir, "solver": marker.solver, "frames": n_frames, "end_time": end_time,
        "cells": cells, "spacing": spacing, "fluid_cells": fluid_cells, "openfoam_cells": centres.len(), "max_speed": max_speed});
    if compressible {
        reply["max_mach"] = json!(max_mach);
    }
    Ok(reply)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sample plate in a box of water flowing along x: its step written as a case, then a
    /// hand-made "run" (a uniform flow at two times) read back as the step's fluid results.
    #[test]
    fn writes_a_case_and_reads_its_results_back() {
        use fieldcraft_geometry::sample::{Plate, plate_design};
        let dir = std::env::temp_dir().join(format!("fieldcraft-openfoam-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, plate_design(&Plate::default()).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": design.to_str().unwrap()})).unwrap();
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let [lo, hi] = [0, 1].map(|k| bodies[0]["bbox"][k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>());
        let min: Vec<f64> = lo.iter().map(|v| (v / 5.0).floor() * 5.0 - 20.0).collect();
        let max: Vec<f64> = hi.iter().map(|v| (v / 5.0).ceil() * 5.0 + 20.0).collect();
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-9}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [1000, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0},
                          "z_min": {"type": "slip"}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10}}},
        ]);
        s.run_script(&script).unwrap();
        let journal = s.journal.len();
        let case = dir.join("case");
        let out = s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "processes": 4, "cell_factor": 1})).unwrap();
        assert_eq!(s.journal.len(), journal, "writing a case does not change the model");
        assert_eq!(out["solver"], "pimpleFoam");
        let read_text = |rel: &str| std::fs::read_to_string(case.join(rel)).unwrap();
        // Lengths in metres: the 5 mm spacing is the background cell.
        let extent: Vec<f64> = (0..3).map(|i| (max[i] - min[i]) / 5.0).collect();
        assert_eq!(out["background_cells"], json!(extent.iter().map(|n| n.round() as u64).collect::<Vec<_>>()), "{out}");
        assert!(read_text("system/blockMeshDict").contains(&format!("({:e} ", min[0] / 1000.0)));
        // A solid and a patch per named face of the plate; m/s and m²/s.
        let faces = bodies[0]["faces"].as_array().unwrap();
        let stl = read_text("constant/triSurface/Body1.stl");
        assert_eq!(stl.matches("endsolid").count(), faces.len(), "{out}");
        assert_eq!(out["patches"].as_object().unwrap().len(), 6 + faces.len());
        assert!(read_text("0/U").contains("x_min { type fixedValue; value uniform (1e0 0 0); }"));
        assert!(read_text("constant/transportProperties").contains("nu 1e-6;"));
        assert!(read_text("system/decomposeParDict").contains("numberOfSubdomains 4;"));
        assert!(read_text(MARKER).contains("\"step\": \"Flow\""));

        // Not run yet: nothing to read.
        assert!(s.execute("results.openfoam", &json!({"path": case.to_str().unwrap()})).is_err());
        // A "run": eight OpenFOAM cells at the corners of the box's inner half, and a flow at two times.
        let lo_m: Vec<f64> = min.iter().map(|v| v / 1000.0).collect();
        let hi_m: Vec<f64> = max.iter().map(|v| v / 1000.0).collect();
        let mut centres = Vec::new();
        for k in 0..8 {
            let at = |a: usize| if k >> a & 1 == 1 { lo_m[a] + 0.75 * (hi_m[a] - lo_m[a]) } else { lo_m[a] + 0.25 * (hi_m[a] - lo_m[a]) };
            centres.push(format!("({} {} {})", at(0), at(1), at(2)));
        }
        let vectors = |items: &[String]| format!("internalField nonuniform List<vector> {}({});\nboundaryField {{}}\n", items.len(), items.join(" "));
        std::fs::write(case.join("0/C"), vectors(&centres)).unwrap();
        for (t, ux) in [("0.25", 0.5), ("0.5", 1.0)] {
            std::fs::create_dir_all(case.join(t)).unwrap();
            let u: Vec<String> = (0..8).map(|k| format!("({} 0 0)", if k & 1 == 1 { 2.0 * ux } else { ux })).collect();
            std::fs::write(case.join(t).join("U"), vectors(&u)).unwrap();
            std::fs::write(case.join(t).join("p"), "internalField uniform 2;\n").unwrap();
        }
        let r = s.execute("results.openfoam", &json!({"path": case.to_str().unwrap(), "spacing": 5})).unwrap();
        assert_eq!(r["frames"], 2, "{r}");
        assert_eq!(r["max_speed"], 2.0);
        assert!(
            r["fluid_cells"].as_u64().unwrap() < r["cells"].as_array().unwrap().iter().map(|c| c.as_u64().unwrap()).product::<u64>(),
            "the plate is solid: {r}"
        );
        // Probed through the ordinary fluid commands: the nearest cell's flow, p × ρ as Pa.
        let probe_at = |s: &mut Session, x: f64, frame: u64| {
            let pt = [x, (min[1] + max[1]) / 2.0 - 0.3 * (max[1] - min[1]), (min[2] + max[2]) / 2.0 - 0.3 * (max[2] - min[2])];
            s.execute("solve.probe", &json!({"step": "Flow", "points": [pt], "frame": frame})).unwrap()
        };
        let low = probe_at(&mut s, min[0] + 0.1 * (max[0] - min[0]), 1);
        let high = probe_at(&mut s, min[0] + 0.9 * (max[0] - min[0]), 0);
        assert!(low.to_string().contains("1000") || low.to_string().contains("1e3") || low.to_string().contains("1000.0"), "{low}");
        assert!(high.to_string().contains("1000"), "frame 0 at 0.5 m/s × 2 = 1 m/s = 1000 mm/s: {high}");
        assert_eq!(s.journal.len(), journal, "reading results does not change the model");

        // Compressible: each cell's density from the ideal gas, ρ = p / (R T), carried by the
        // frames; probes show it and the ports' mass flow uses it, not the material's.
        let gas = dir.join("gas");
        let out = s.execute("deck.openfoam", &json!({"path": gas.to_str().unwrap(), "cell_factor": 1, "compressible": true})).unwrap();
        assert_eq!(out["solver"], "rhoCentralFoam");
        std::fs::write(gas.join("0/C"), vectors(&centres)).unwrap();
        std::fs::create_dir_all(gas.join("0.5")).unwrap();
        let u: Vec<String> = (0..8).map(|_| "(2 0 0)".to_string()).collect();
        std::fs::write(gas.join("0.5/U"), vectors(&u)).unwrap();
        std::fs::write(gas.join("0.5/p"), "internalField uniform 2e5;\n").unwrap();
        std::fs::write(gas.join("0.5/T"), "internalField uniform 300;\n").unwrap();
        s.execute("results.openfoam", &json!({"path": gas.to_str().unwrap(), "spacing": 5})).unwrap();
        let rho = 2e5 / (fieldcraft_decks::openfoam::GAS_CONSTANT / fieldcraft_decks::openfoam::Gas::AIR.molar_mass * 300.0);
        let display = s.units().clone();
        use fieldcraft_model::units::Quantity as Q;
        let probe = probe_at(&mut s, min[0] + 0.9 * (max[0] - min[0]), 0);
        let shown = probe["points"][0]["density"].as_f64().unwrap();
        let expected = display.express(Q::Density, rho).unwrap();
        assert!((shown / expected - 1.0).abs() < 1e-5, "{probe}");
        let report = s.execute("results.fluid.measure", &json!({})).unwrap();
        let row = |q: &str| {
            report["rows"].as_array().unwrap().iter().find(|r| r["face"] == "x_max" && r["quantity"] == q).unwrap()["value"].as_f64().unwrap()
        };
        let si = |q, v| display.to_si(q, v).unwrap();
        let ratio = si(Q::MassFlow, row("outflow_mass_flow")) / si(Q::VolumeFlow, row("outflow_volume_flow"));
        assert!((ratio / rho - 1.0).abs() < 1e-5, "mass flow at ρ = {rho}, not the water's: {ratio}");

        // What OpenFOAM cannot be given yet is named, and nothing is written.
        s.execute(
            "model.fluid_domain.edit",
            &json!({"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"], "acceleration": [0, 0, -9810]}),
        )
        .unwrap();
        let other = dir.join("other");
        let e = s.execute("deck.openfoam", &json!({"path": other.to_str().unwrap()})).unwrap_err().to_string();
        assert!(e.contains("acceleration"), "{e}");
        assert!(!other.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn nearest_cells_search_beyond_empty_buckets() {
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [5, 1, 1] };
        let solid = [false, false, true, false, false];
        let centres = [[0.4, 0.5, 0.5], [4.9, 0.5, 0.5]];
        let n = nearest_cells(&grid, &solid, &centres).unwrap();
        assert_eq!(n, vec![Some(0), Some(0), None, Some(1), Some(1)]);
    }
}
