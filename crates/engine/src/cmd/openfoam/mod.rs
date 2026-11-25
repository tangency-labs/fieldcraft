//! OpenFOAM: write a fluid step as a case directory (`deck.openfoam`), an export for running the
//! case in OpenFOAM, which runs as its own program. The case is written from OpenFOAM's public
//! user guide (`fieldcraft_decks::openfoam`); FieldCraft's own fluid results are checked against
//! theory and published benchmarks (docs/verification.md), never against OpenFOAM's.

mod location;
pub mod read;

/// The record `deck.openfoam` writes beside the case, from which `results.openfoam` reads it
/// back: what the case was written with.
pub(crate) const MARKER: &str = "fieldcraft.json";

/// The marker's content. Every field is required: a marker that lacks one is refused, never
/// filled in.
#[derive(serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Marker {
    /// The fluid step the case was written from.
    pub step: String,
    /// The OpenFOAM solver application.
    pub solver: String,
    /// The fluid's density (kg/m³): the case's pressure is kinematic, p/ρ.
    pub density: f64,
    /// The background cell (m).
    pub background: f64,
    /// The lattice spacing the surfaces are refined to (m).
    pub spacing: f64,
    /// Significant digits of the time directories' names.
    pub time_precision: u32,
    /// A compressible case's gas, in SI: what reading its results back needs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gas: Option<MarkerGas>,
}

/// A compressible case's gas, in SI: what reading its results back needs.
#[derive(serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct MarkerGas {
    /// The ambient (absolute) pressure the gauge pressures count from (Pa).
    pub ambient: f64,
    /// The initial, inflow and reservoir temperature (K).
    pub temperature: f64,
    /// The specific gas constant R/M (J/(kg·K)).
    pub gas_constant: f64,
    /// The specific heat at constant pressure (J/(kg·K)).
    pub specific_heat: f64,
}

impl MarkerGas {
    /// The ratio of specific heats, γ = c_p / (c_p − R_s), derived rather than stored so it
    /// cannot disagree with them.
    pub fn gamma(&self) -> f64 {
        self.specific_heat / (self.specific_heat - self.gas_constant)
    }
}

impl Marker {
    /// The most bytes a marker may have: it holds a step name (at most
    /// [`fieldcraft_model::MAX_NAME_LEN`] characters), a solver name and at most eight numbers, a
    /// few hundred bytes.
    pub const MAX_BYTES: u64 = 4096;

    /// Refuse a marker whose numbers cannot be what `deck.openfoam` writes.
    pub fn check(&self) -> Result<()> {
        let positive = |v: f64| v.is_finite() && v > 0.0;
        // An f64 has at most 17 significant decimal digits (IEEE 754-2008 §5.12.2).
        const F64_DIGITS: u32 = 17;
        if !(positive(self.density) && positive(self.background) && positive(self.spacing)) || !(1..=F64_DIGITS).contains(&self.time_precision) {
            return Err(crate::EngineError::invalid(
                "path",
                format!("{MARKER} holds a density, cell size or time precision deck.openfoam does not write"),
            ));
        }
        match (self.solver.as_str(), &self.gas) {
            (fieldcraft_decks::openfoam::INCOMPRESSIBLE_SOLVER, None) => Ok(()),
            (fieldcraft_decks::openfoam::gas::SOLVER, Some(g)) => {
                if positive(g.ambient)
                    && positive(g.temperature)
                    && positive(g.gas_constant)
                    && positive(g.specific_heat)
                    && g.specific_heat > g.gas_constant
                {
                    Ok(())
                } else {
                    Err(crate::EngineError::invalid(
                        "path",
                        format!(
                            "{MARKER} holds a gas deck.openfoam does not write (its pressure, temperature, gas constant and specific heat, with cp above R)"
                        ),
                    ))
                }
            }
            (fieldcraft_decks::openfoam::INCOMPRESSIBLE_SOLVER, Some(_)) => {
                Err(crate::EngineError::invalid("path", format!("{MARKER} records a gas for the incompressible solver pimpleFoam")))
            }
            (fieldcraft_decks::openfoam::gas::SOLVER, None) => {
                Err(crate::EngineError::invalid("path", format!("{MARKER} names the compressible solver rhoCentralFoam but records no gas")))
            }
            (other, _) => {
                Err(crate::EngineError::invalid("path", format!("{MARKER} names the solver `{other}`, which deck.openfoam does not write")))
            }
        }
    }
}
mod nearest;

use std::path::Path;

use fieldcraft_decks::openfoam::{CaseOptions, GasConditions, MeshOptions, ObstacleSurface, PressureBox, RefineBox, gas, write_case};
use fieldcraft_io::roots::Access;
use fieldcraft_lbm::Obstacle;
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use serde_json::{Value, json};

use super::opt;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, Result, Session};

pub const WRITE: CommandSpec = CommandSpec {
    id: "deck.openfoam",
    label: "Write OpenFOAM Case",
    params: r#"{"path": "case_dir", "step"?: "Flow", "cell_factor"?: 1, "processes"?: 1, "max_courant"?: 1, "compressible"?: {"ambient": p, "temperature": T, "gas_constant": R, "far_field": L, "initial_pressure"?: [{"min": [x,y,z], "max": [x,y,z], "pressure": p}]}, "refine_boxes"?: [{"min": [x,y,z], "max": [x,y,z], "spacing"?: 0.1}], "location"?: [x,y,z]} — writes a fluid step (the only one when `step` is left out) as an OpenFOAM case directory, for running it in OpenFOAM: a blockMesh background of spacing × cell_factor (1, 2, 4, 8), snappyHexMesh around the obstacles (one patch per named face, refined back to the spacing; `refine_boxes` to theirs, rounded to the background over a power of two), ASCII fields, pimpleFoam with the step's turbulence model and a time step that adapts to max_courant (at most 1), decomposeParDict for `processes` ranks and an Allrun script (mesh, cell centres, solve, reconstruct). Every setting is written into the case. `location` is a point in the fluid, the side snappyHexMesh keeps: found when left out; when given, it is moved to the centre of the finest cell holding it (on no cell face) and refused inside an obstacle. The case directory must be new or empty. With `compressible`: rhoCentralFoam, a laminar ideal gas: the viscosity, specific heat and conductivity are the material's (each required, constant, the conductivity one value), `gas_constant` its specific gas constant R/M, `ambient` the absolute pressure the model's gauge pressures count from, `temperature` the initial, inflow and reservoir temperature and `far_field` the far-field relaxation length of the non-reflecting outlets (OpenFOAM's waveTransmissive lInf: the distance over which an outlet's pressure relaxes towards the ambient), none with a default; a pressure face above 0 is a reservoir at that total pressure and temperature, the others non-reflecting outlets, walls adiabatic; max_courant then defaults to 0.25 (at most 0.5); `initial_pressure` boxes start at their gauge pressure (setFields). A turbulence model with `compressible` is refused. Lengths in the input units; the case is SI. Fails naming what the case cannot express yet (heat, porous zones, fans, baffles, free surface, periodic faces, parabolic inlets, non-Newtonian or temperature-dependent fluids) and any fluid property the material lacks. Nothing of the model changes"#,
    form: Some(super::forms::DECK_OPENFOAM),
    kind: CommandKind::Session,
    requires: &[Requirement::FluidDomain, Requirement::Step],
    run: write,
};

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct BoxParam {
    min: [f64; 3],
    max: [f64; 3],
    #[serde(default)]
    spacing: Option<f64>,
}

/// A box of the initial pressure field: its gauge pressure inside it.
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct PressureParam {
    min: [f64; 3],
    max: [f64; 3],
    pressure: f64,
}

/// The gas of a compressible case. Every number is required: a missing one is refused, never
/// filled in.
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct GasParam {
    ambient: f64,
    temperature: f64,
    gas_constant: f64,
    far_field: f64,
    #[serde(default)]
    initial_pressure: Vec<PressureParam>,
}

/// The gas conditions of a `compressible` parameter; `None` for an incompressible case.
fn gas_conditions(p: &Value) -> Result<Option<GasConditions>> {
    let Some(g) = opt::<GasParam>(p, "compressible")? else {
        return Ok(None);
    };
    if g.initial_pressure.len() > gas::MAX_PRESSURE_BOXES {
        return Err(crate::EngineError::invalid("compressible.initial_pressure", format!("at most {} boxes", gas::MAX_PRESSURE_BOXES)));
    }
    Ok(Some(GasConditions {
        ambient: g.ambient,
        temperature: g.temperature,
        gas_constant: g.gas_constant,
        far_field: g.far_field,
        initial_pressure: g.initial_pressure.into_iter().map(|b| PressureBox { min: b.min, max: b.max, pressure: b.pressure }).collect(),
    }))
}

/// `text` as one POSIX shell word: in single quotes, each `'` written as `'\''`.
fn shell_quote(text: &str) -> String {
    format!("'{}'", text.replace('\'', "'\\''"))
}

fn write(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = super::file_param(s, p, "path", Access::Write)?;
    let root = Path::new(&dir);
    // Nothing of an earlier case is overwritten or mixed in: a directory holding anything is refused.
    if !fieldcraft_io::fs::is_empty_dir(root)? {
        return Err(crate::EngineError::invalid("path", "the case directory is not empty: choose a new or empty directory"));
    }
    let step: Option<String> = opt(p, "step")?;
    let (step, domain) = fieldcraft_lbm::fluid_step(&s.model, step.as_deref())?;
    let step_name = step.name.clone();
    let gas = gas_conditions(p)?;
    // What the case cannot express is refused before the obstacles are tessellated and searched.
    fieldcraft_decks::openfoam::check(&s.model, Some(&step_name), gas.as_ref())?;
    let mut mesh = MeshOptions::default();
    if let Some(f) = opt(p, "cell_factor")? {
        mesh.cell_factor = f;
    }
    let boxes = opt::<Vec<BoxParam>>(p, "refine_boxes")?.unwrap_or_default();
    if boxes.len() > fieldcraft_decks::openfoam::MAX_REFINE_BOXES {
        return Err(crate::EngineError::invalid("refine_boxes", format!("at most {} boxes", fieldcraft_decks::openfoam::MAX_REFINE_BOXES)));
    }
    mesh.refine_boxes = boxes.into_iter().map(|b| RefineBox { min: b.min, max: b.max, spacing: b.spacing }).collect();

    // The obstacles, tessellated once as the lattice solver does, in metres: their triangles go
    // into the STL, and their closed surfaces decide which side of them `location` is on.
    let tolerance = super::solve::obstacle_tolerance(domain);
    let mut meshes = Vec::with_capacity(domain.obstacles.len());
    let mut solids = Vec::with_capacity(domain.obstacles.len());
    for body in &domain.obstacles {
        let (mesh, faces) = super::solve::obstacle_mesh(s, body, tolerance)?;
        let positions: Vec<[f64; 3]> = mesh.positions.iter().map(|q| crate::units::geometry_point(*q)).collect();
        let weld = SurfaceOptions { weld_distance: crate::units::geometry_length(mesh.weld_distance), ..SurfaceOptions::default() };
        let surface = BoundarySurface::new(&positions, &mesh.triangles, &mesh.tri_face, &weld)
            .map_err(|e| crate::EngineError::context(format!("obstacle `{body}`: its surface cannot bound a solid"), e))?;
        solids.push(Obstacle { body: body.clone(), surface, faces: faces.clone(), exact: Vec::new() });
        meshes.push((body.as_str(), positions, mesh.into_owned(), faces));
    }
    let given = opt::<[f64; 3]>(p, "location")?;
    let location = match given {
        Some(l) => location::given(domain, &solids, &mesh, l)?,
        None => location::find(domain, &solids, &mesh)?,
    };
    let mut options = CaseOptions::new(mesh, location);
    if let Some(n) = opt(p, "processes")? {
        options.processes = n;
    }
    options.gas = gas;
    if let Some(c) = opt(p, "max_courant")? {
        options.max_courant = c;
    } else if options.gas.is_some() {
        options.max_courant = gas::COURANT;
    }
    let surfaces: Vec<ObstacleSurface<'_>> = meshes
        .iter()
        .map(|(body, positions, mesh, faces)| ObstacleSurface { body, positions, triangles: &mesh.triangles, tri_face: &mesh.tri_face, faces })
        .collect();
    let mut case = write_case(&s.model, Some(&step_name), &surfaces, &options)?;
    if given.is_some_and(|l| l != location) {
        case.notes.push("`location` moved to the centre of the finest cell holding it, which lies on no cell face".into());
    }
    let marker = Marker {
        step: step_name.clone(),
        solver: case.solver.to_string(),
        density: case.density,
        background: case.background,
        spacing: domain.spacing,
        time_precision: case.time_precision,
        gas: case.gas.map(|g| MarkerGas {
            ambient: g.ambient,
            temperature: g.temperature,
            gas_constant: g.gas_constant,
            specific_heat: g.specific_heat,
        }),
    };
    let marker_text = serde_json::to_string_pretty(&marker).map_err(|e| crate::EngineError::Params(e.to_string()))?;
    let mut bytes = 0;
    for (rel, text) in case.files.iter().map(|(k, v)| (k.as_str(), v.as_str())).chain([(MARKER, marker_text.as_str())]) {
        let path = root.join(rel);
        if let Some(parent) = path.parent() {
            fieldcraft_io::fs::create_dir_all(parent)?;
        }
        fieldcraft_io::fs::write_output(&path, text.as_bytes())?;
        bytes += text.len();
    }
    Ok(json!({"format": "openfoam", "path": dir, "step": step_name, "solver": case.solver, "files": case.files.keys().collect::<Vec<_>>(),
        "bytes": bytes, "background_cells": case.background_cells, "location": options.location, "patches": case.patches,
        "gas": case.gas.map(|g| json!({"gamma": g.gamma(), "prandtl": g.prandtl()})),
        "notes": case.notes, "run": format!("cd {} && sh Allrun", shell_quote(&dir))}))
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_gas_read;
#[cfg(test)]
mod tests_gas_write;
