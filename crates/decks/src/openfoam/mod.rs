//! A fluid step written as an OpenFOAM case, for running the case in OpenFOAM: the
//! dictionaries of a `blockMesh` + `snappyHexMesh` mesh of the fluid domain around its
//! obstacles, the initial and boundary fields, and the controls of an incompressible transient
//! run (`pimpleFoam`), as a map of case-relative paths to text, for [`TARGET`]. The formats come
//! from OpenFOAM's public documentation: the ESI user guide and tutorial guide (openfoam.com) and
//! the documentation site doc.openfoam.com (case structure, dictionary syntax, boundary
//! conditions, each entry's section cited where it is written); nothing of OpenFOAM's source or
//! tutorial cases is read or copied. Every setting the case depends on is
//! written into it (no `#include` of the installation's files), so the same case meshes and runs
//! the same on every OpenFOAM installation that reads it.
//!
//! - The domain box is the background block (cells of `spacing × cell_factor`); its six faces
//!   are the patches `x_min` … `z_max`.
//! - Each obstacle body is a multi-solid ASCII STL in `constant/triSurface`, one solid per named
//!   face, so `snappyHexMesh` makes a patch per face (`<body>_<face>`, sanitised) and the
//!   obstacle's surface is refined back to the lattice spacing.
//! - Face conditions: velocity → fixed velocity; pressure → fixed kinematic pressure (p/ρ) with
//!   an inflow/outflow velocity; slip → a slip wall; symmetry → `symmetryPlane` (box faces are
//!   planar); wall → no-slip; moving wall → fixed velocity. Obstacle faces are no-slip unless the
//!   domain names them.
//! - The step's length and frames give `endTime` and the write interval; the time step adapts
//!   to a Courant limit (the lattice's time step is the LBM's own). The step's turbulence model
//!   is written as OpenFOAM's ([`dicts`]); a running mean starts where the step's own averaging
//!   does.
//!
//! Everything is SI, as the model stores it and OpenFOAM reads it: no unit conversion happens
//! here. What the case cannot express yet is refused by name ([`crate::Error::Unsupported`]); a
//! material property the case needs and the model lacks is an error, never a stand-in value.

mod dicts;
pub mod fields;
pub mod gas;
pub use dicts::SOLVER as INCOMPRESSIBLE_SOLVER;
pub use dicts::sample_reach;
pub use gas::{Gas, GasConditions, PressureBox};
#[cfg(test)]
mod grammar;
#[cfg(test)]
mod tests;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use fieldcraft_model::{BoxFace, FaceCondition, Fluid, FluidDomain, Model, Procedure, Profile, Step, SurfaceCondition};

use crate::{Error, Result};

const FORMAT: &str = "OpenFOAM";

/// The OpenFOAM line and release the case is written for, whose public user guide
/// (doc.openfoam.com) the dictionaries follow; stated in every file's header.
pub const TARGET: &str = "ESI OpenFOAM v2406 or later (openfoam.com)";

/// The most triangles the obstacle surfaces may have together: a cap on a hostile input, sized by
/// the text it makes. One ASCII facet is [`FACET_TEXT_BYTES`] bytes at most at the precision
/// written (a `facet normal` line, `outer loop`, three `vertex` lines, `endloop`, `endfacet`;
/// nine numbers of at most 24 characters each), so this many triangles are about 1.5 GB of text
/// held in memory before it is written.
pub const MAX_TRIANGLES: usize = 5_000_000;

/// The most bytes one facet of the STL text takes (see [`MAX_TRIANGLES`]).
#[cfg(test)]
const FACET_TEXT_BYTES: usize = 300;

/// The most MPI ranks a case is decomposed for: a cap on a hostile input, far above any one
/// machine's or cluster partition's core count a case of at most [`MAX_TRIANGLES`] meshes on.
pub const MAX_PROCESSES: u32 = 65_536;

/// The cell factors a case takes: background cells of the spacing times a power of two, so the
/// refinement levels back to the spacing are whole (each `snappyHexMesh` level halves a cell).
pub const CELL_FACTORS: [u32; 4] = [1, 2, 4, 8];

/// The most refinement boxes a case takes: a cap on a hostile input (each box is a geometry
/// entry and a refinement region of `snappyHexMeshDict`).
pub const MAX_REFINE_BOXES: usize = 1024;

/// The most refinement levels below the background a refinement box may ask for: one background
/// cell refined to level L makes 8^L cells, and the largest L with 8^L within the 32-bit label
/// cap ([`dicts::MAX_LABEL`], 2³¹ − 1) is `⌊log2(MAX_LABEL) / 3⌋`.
pub const MAX_BOX_LEVEL: u32 = dicts::MAX_LABEL.ilog2() / 3;

/// One obstacle body's surface, in metres: its triangles, the face each triangle belongs to
/// (an index into `faces`) and the faces' persistent names.
#[derive(Debug, Clone, Copy)]
pub struct ObstacleSurface<'a> {
    pub body: &'a str,
    pub positions: &'a [[f64; 3]],
    pub triangles: &'a [[u32; 3]],
    pub tri_face: &'a [u32],
    pub faces: &'a [String],
}

/// A box refined to `spacing` (metres; the lattice spacing when left out, rounded to the
/// background cell over a power of two).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
    pub spacing: Option<f64>,
}

/// How the case is meshed.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshOptions {
    /// The background cells are the domain's spacing times this (one of [`CELL_FACTORS`]);
    /// obstacle surfaces and refinement boxes are refined back to the spacing.
    pub cell_factor: u32,
    pub refine_boxes: Vec<RefineBox>,
}

impl Default for MeshOptions {
    fn default() -> Self {
        MeshOptions { cell_factor: 1, refine_boxes: Vec::new() }
    }
}

/// How the case is meshed and run.
#[derive(Debug, Clone, PartialEq)]
pub struct CaseOptions {
    pub mesh: MeshOptions,
    /// A point in the fluid (metres) that tells `snappyHexMesh` which side to keep. It must lie
    /// in the fluid and on no cell face of any refinement level ([`finest_centre`] gives such
    /// points).
    pub location: [f64; 3],
    /// MPI ranks `Allrun` runs the solver on (1: serial; at most [`MAX_PROCESSES`]).
    pub processes: u32,
    /// The Courant number the time step adapts to ([`dicts::MAX_COURANT`] by default; a
    /// compressible case's at most [`gas::MAX_COURANT`]).
    pub max_courant: f64,
    /// Compressible flow of an ideal gas (`rhoCentralFoam`) in these conditions, instead of
    /// incompressible flow (`pimpleFoam`).
    pub gas: Option<GasConditions>,
}

impl CaseOptions {
    /// The options with their defaults, meshed as `mesh`, around the point `location` in the
    /// fluid.
    pub fn new(mesh: MeshOptions, location: [f64; 3]) -> Self {
        CaseOptions { mesh, location, processes: 1, max_courant: dicts::MAX_COURANT, gas: None }
    }
}

/// A written case.
#[derive(Debug, Clone, PartialEq)]
pub struct Case {
    /// Case-relative path (with `/`) → text, in path order.
    pub files: BTreeMap<String, String>,
    /// OpenFOAM patch → what it is in the model (`x_min`, or `Body1/F2:top`).
    pub patches: BTreeMap<String, String>,
    /// The background block's cells along x, y and z.
    pub background_cells: [u64; 3],
    /// The solver application.
    pub solver: &'static str,
    /// What the case leaves out or does differently from the lattice Boltzmann step.
    pub notes: Vec<String>,
    /// The fluid's density (kg/m³): OpenFOAM's incompressible pressure is kinematic, p/ρ; a
    /// compressible case's reference density, the ambient pressure over R_s T.
    pub density: f64,
    /// The gas of a compressible case.
    pub gas: Option<Gas>,
    /// The background cell (metres).
    pub background: f64,
    /// Significant digits of the time directories' names.
    pub time_precision: u32,
}

fn unsupported(what: impl Into<String>) -> Error {
    Error::Unsupported { format: FORMAT, what: what.into() }
}

/// The background cell and the finest cell of a case (metres): the spacing times the cell
/// factor, and that cell halved once per level of the finest refinement (obstacle surfaces at
/// the spacing, refinement boxes at their own).
pub fn cell_sizes(domain: &FluidDomain, mesh: &MeshOptions) -> Result<(f64, f64)> {
    if !CELL_FACTORS.contains(&mesh.cell_factor) {
        return Err(unsupported(format!("a cell factor other than {CELL_FACTORS:?}")));
    }
    let background = domain.spacing * f64::from(mesh.cell_factor);
    Ok((background, background / f64::from(1u32 << finest_level(domain, mesh)?)))
}

/// The finest refinement level below the background cell: the obstacle surfaces' (the cell
/// factor's power of two, back to the spacing) or a refinement box's, whichever is deeper.
fn finest_level(domain: &FluidDomain, mesh: &MeshOptions) -> Result<u32> {
    Ok(levels(domain, mesh)?.into_iter().fold(mesh.cell_factor.trailing_zeros(), u32::max))
}

/// The centre of the finest cell `k` (counted along each axis from the domain's lower corner):
/// `min + (k + ½) h_finest`, with `h_finest` the block's extent over its finest cell count, as
/// `blockMesh` divides it (the domain's extent may differ from a whole number of spacings by the
/// tolerance of [`FluidDomain::cells`], which the spacing would carry into the centre).
/// `snappyHexMesh` needs its location inside a cell, never on a face, and every refinement
/// level's faces lie on whole multiples of the finest cell from the domain's lower corner, so such
/// a centre is on none of them.
pub fn finest_centre(domain: &FluidDomain, mesh: &MeshOptions, k: [u64; 3]) -> Result<[f64; 3]> {
    let n = finest_cells(domain, mesh)?;
    Ok([0, 1, 2].map(|a| domain.min[a] + (k[a] as f64 + 0.5) * (domain.max[a] - domain.min[a]) / n[a] as f64))
}

/// The finest cells along each axis of the block: each background cell is `cell_factor = 2^l`
/// lattice cells and halves down to the finest level.
pub fn finest_cells(domain: &FluidDomain, mesh: &MeshOptions) -> Result<[u64; 3]> {
    cell_sizes(domain, mesh)?;
    let lattice = domain.cells()?;
    let per_lattice_cell = 1u64 << (finest_level(domain, mesh)? - mesh.cell_factor.trailing_zeros());
    let mut n = [0u64; 3];
    for (a, slot) in n.iter_mut().enumerate() {
        *slot = lattice[a].checked_mul(per_lattice_cell).ok_or_else(|| Error::Limit("a finest cell count beyond 64 bits".into()))?;
    }
    Ok(n)
}

/// The finest cell holding `p` (counted from the domain's lower corner): a point on a face
/// between two cells is in the upper one. A point not strictly inside the domain box is refused.
pub fn finest_cell_of(domain: &FluidDomain, mesh: &MeshOptions, p: [f64; 3]) -> Result<[u64; 3]> {
    if p.iter().any(|v| !v.is_finite()) || (0..3).any(|a| p[a] <= domain.min[a] || p[a] >= domain.max[a]) {
        return Err(unsupported("a location that is not inside the domain box"));
    }
    let n = finest_cells(domain, mesh)?;
    let mut k = [0u64; 3];
    for (a, slot) in k.iter_mut().enumerate() {
        let t = ((p[a] - domain.min[a]) / (domain.max[a] - domain.min[a]) * n[a] as f64).floor();
        // Inside the box, so 0 ≤ t ≤ n; round-off at the upper face keeps it in the last cell.
        *slot = (t as u64).min(n[a].saturating_sub(1));
    }
    Ok(k)
}

/// Whether the case can be written at all: the step, its fluid domain, its material and its
/// turbulence model are all ones the case expresses. `write_case` makes the same checks; this
/// lets a caller refuse before costly preparation (such as finding the location).
pub fn check(model: &Model, step: Option<&str>, gas: Option<&GasConditions>) -> Result<()> {
    let step = fluid_step(model, step)?;
    let Procedure::Fluid { domain: domain_name, turbulence, .. } = &step.procedure else {
        return Err(unsupported(format!("step `{}` is not a fluid step", step.name)));
    };
    let domain =
        model.fluid_domains().get(domain_name).ok_or_else(|| unsupported(format!("step `{}` names no fluid domain `{domain_name}`", step.name)))?;
    refuse_unsupported(domain, &mut Vec::new())?;
    match gas {
        Some(conditions) => {
            refuse_turbulence_with_gas(turbulence.as_ref())?;
            Gas::of(model, domain, conditions)?;
        }
        None => {
            fluid_properties(model, domain)?;
            dicts::turbulence_properties(turbulence.as_ref(), &mut Vec::new())?;
        }
    }
    Ok(())
}

/// A compressible case is laminar: a turbulence model is refused by name, never dropped.
fn refuse_turbulence_with_gas(turbulence: Option<&fieldcraft_model::Turbulence>) -> Result<()> {
    match turbulence {
        Some(t) => Err(unsupported(format!("the {} turbulence model in a compressible case (laminar only so far)", t.name()))),
        None => Ok(()),
    }
}

/// Each refinement box's level below the background: the box's spacing (the lattice spacing
/// when left out) as the background cell over a power of two, the nearest on a logarithmic
/// scale.
fn levels(domain: &FluidDomain, mesh: &MeshOptions) -> Result<Vec<u32>> {
    if mesh.refine_boxes.len() > MAX_REFINE_BOXES {
        return Err(Error::Limit(format!("{} refinement boxes (at most {MAX_REFINE_BOXES})", mesh.refine_boxes.len())));
    }
    let background = domain.spacing * f64::from(mesh.cell_factor);
    let mut out = Vec::with_capacity(mesh.refine_boxes.len());
    for b in &mesh.refine_boxes {
        if b.min.iter().chain(&b.max).any(|v| !v.is_finite()) || (0..3).any(|i| b.min[i] >= b.max[i]) {
            return Err(unsupported("a refinement box whose min is not below its max"));
        }
        let h = b.spacing.unwrap_or(domain.spacing);
        if !(h.is_finite() && h > 0.0) {
            return Err(unsupported("a refinement box spacing that is not positive"));
        }
        let l = libm::log2(background / h).round();
        if !(0.0..=f64::from(MAX_BOX_LEVEL)).contains(&l) {
            return Err(unsupported(format!("a refinement box spacing that is not the background cell over 1 to 2^{MAX_BOX_LEVEL}")));
        }
        out.push(l as u32);
    }
    Ok(out)
}

/// The fluid step `step` (the model's only fluid step when `None`) as an OpenFOAM case.
pub fn write_case(model: &Model, step: Option<&str>, obstacles: &[ObstacleSurface<'_>], options: &CaseOptions) -> Result<Case> {
    let step = fluid_step(model, step)?;
    let Procedure::Fluid { domain: domain_name, end_time, steps, time_step, frames, turbulence, steady, .. } = &step.procedure else {
        return Err(unsupported(format!("step `{}` is not a fluid step", step.name)));
    };
    let domain =
        model.fluid_domains().get(domain_name).ok_or_else(|| unsupported(format!("step `{}` names no fluid domain `{domain_name}`", step.name)))?;
    let mut notes = Vec::new();
    refuse_unsupported(domain, &mut notes)?;
    let gas = match &options.gas {
        Some(conditions) => {
            refuse_turbulence_with_gas(turbulence.as_ref())?;
            notes.push("compressible: an ideal gas from the material's specific heat, viscosity and conductivity; its density is not used".into());
            let gas = Gas::of(model, domain, conditions)?;
            if gas.gamma() > gas::MAX_IDEAL_GAMMA {
                notes.push(format!(
                    "the gas's ratio of specific heats γ = {:.4} is above 5/3, which an ideal gas does not reach (c_v >= 3R/2): the case is written from the data as given",
                    gas.gamma()
                ));
            }
            Some(gas)
        }
        None => None,
    };
    let (density, viscosity) = match &gas {
        Some(g) => (g.reference_density(), g.viscosity),
        None => fluid_properties(model, domain)?,
    };
    let duration = match (end_time, steps, time_step) {
        (Some(t), _, _) => *t,
        (None, Some(n), Some(dt)) => *n as f64 * dt,
        _ => return Err(unsupported("a step length in lattice steps without a time step (give end_time)")),
    };
    if !(duration.is_finite() && duration > 0.0) {
        return Err(unsupported("a step length that is not a positive time"));
    }
    if time_step.is_some() {
        notes.push("the lattice time step is not used: OpenFOAM's adapts to the Courant limit".into());
    }
    let average_from = steady.as_ref().and_then(|s| s.averaging).map(|a| a.start);
    if steady.is_some() {
        notes.push("the steady-state stop rule is not mapped: the run goes to the end time".into());
    }
    let (_, finest) = cell_sizes(domain, &options.mesh)?;
    let diagonal = (0..3).map(|a| (domain.max[a] - domain.min[a]).powi(2)).sum::<f64>().sqrt();
    if finest <= dicts::MERGE_TOLERANCE * diagonal {
        return Err(unsupported(format!("a finest cell below the mesher's merge tolerance ({} of the domain's diagonal)", dicts::MERGE_TOLERANCE)));
    }
    let max_courant = if gas.is_some() { gas::MAX_COURANT } else { dicts::MAX_COURANT };
    if !(options.max_courant.is_finite() && options.max_courant > 0.0 && options.max_courant <= max_courant) {
        return Err(unsupported(format!("a Courant limit outside (0, {max_courant}]")));
    }
    if !(1..=MAX_PROCESSES).contains(&options.processes) {
        return Err(unsupported(format!("{} processes (1 to {MAX_PROCESSES})", options.processes)));
    }
    if options.location.iter().any(|v| !v.is_finite()) || (0..3).any(|a| options.location[a] <= domain.min[a] || options.location[a] >= domain.max[a])
    {
        return Err(unsupported("a location that is not inside the domain box"));
    }

    let mut files = BTreeMap::new();
    let mut patches = BTreeMap::new();
    for face in BoxFace::ALL {
        patches.insert(face.name().to_string(), face.name().to_string());
    }

    // Obstacles: one STL per body, a solid per named face.
    let total: usize = obstacles.iter().map(|o| o.triangles.len()).sum();
    if total > MAX_TRIANGLES {
        return Err(Error::Limit(format!("{total} obstacle triangles (at most {MAX_TRIANGLES})")));
    }
    let mut surfaces = Vec::new();
    // One name set for the case: the box patches, the refinement boxes' geometry (`refine{i}` in
    // the snappy dictionary), then each body's name and its `<name>_<region>` patches.
    let mut used: BTreeSet<String> = patches.keys().cloned().collect();
    used.extend((0..options.mesh.refine_boxes.len()).map(|i| format!("refine{i}")));
    for o in obstacles {
        if !domain.obstacles.iter().any(|b| b == o.body) {
            continue;
        }
        let (stl, regions) = stl(o)?;
        let name = claim(&sanitize(o.body), &regions, &mut used)?;
        let mut rs = Vec::new();
        for (region, face) in regions {
            let patch = format!("{name}_{region}");
            let condition = surface_condition(domain, o.body, &face);
            patches.insert(patch.clone(), format!("{}/{face}", o.body));
            rs.push(Region { patch, condition });
        }
        files.insert(format!("constant/triSurface/{name}.stl"), stl);
        surfaces.push(Surface { name, regions: rs });
    }
    for body in &domain.obstacles {
        if !obstacles.iter().any(|o| o.body == body) {
            return Err(unsupported(format!("obstacle `{body}` has no surface")));
        }
    }

    let lattice = domain.cells()?;
    let mut cells = [0u64; 3];
    for (axis, c) in cells.iter_mut().enumerate() {
        let n = lattice[axis];
        if n % u64::from(options.mesh.cell_factor) != 0 {
            return Err(unsupported(format!(
                "a cell factor of {} on a box {n} spacings long along {}: the background cells must fit it whole",
                options.mesh.cell_factor,
                ["x", "y", "z"][axis]
            )));
        }
        *c = n / u64::from(options.mesh.cell_factor);
    }
    let level = options.mesh.cell_factor.trailing_zeros();
    let boxes: Vec<(RefineBox, u32)> = options.mesh.refine_boxes.iter().copied().zip(levels(domain, &options.mesh)?).collect();

    files.insert("system/blockMeshDict".into(), block_mesh(domain, cells));
    files.insert("system/surfaceFeatureExtractDict".into(), dicts::feature_extract(&surfaces));
    let background = domain.spacing * f64::from(options.mesh.cell_factor);
    files.insert("system/snappyHexMeshDict".into(), dicts::snappy(&surfaces, level, &boxes, options.location, background, cells));
    let (solver, speed, set_fields) = match (&gas, &options.gas) {
        (Some(g), Some(conditions)) => {
            for (field, text) in gas::fields(domain, g, &surfaces, &mut notes)? {
                files.insert(format!("0/{field}"), text);
            }
            files.insert("constant/thermophysicalProperties".into(), gas::thermophysical(g)?);
            files.insert("constant/turbulenceProperties".into(), dicts::turbulence_properties(None, &mut notes)?);
            files.insert("system/fvSchemes".into(), gas::fv_schemes(dicts::NON_ORTHOGONAL_LIMIT));
            files.insert("system/fvSolution".into(), gas::fv_solution(dicts::solver_tolerance()));
            let boxes = !conditions.initial_pressure.is_empty();
            if boxes {
                files.insert("system/setFieldsDict".into(), gas::set_fields(g, &conditions.initial_pressure)?);
            }
            (gas::SOLVER, gas::courant_speed(domain, g), boxes)
        }
        _ => {
            for (field, text) in fields(domain, density, &surfaces) {
                files.insert(format!("0/{field}"), text);
            }
            files.insert("constant/transportProperties".into(), dicts::transport(viscosity / density));
            files.insert("constant/turbulenceProperties".into(), dicts::turbulence_properties(turbulence.as_ref(), &mut notes)?);
            let les = turbulence.is_some();
            let wall_distance = matches!(turbulence, Some(fieldcraft_model::Turbulence::Smagorinsky { van_driest: Some(_), .. }));
            files.insert("system/fvSchemes".into(), dicts::fv_schemes(les, wall_distance));
            // Orthogonal: no obstacle surface and nothing refined (each box carries its own level).
            let orthogonal = surfaces.is_empty() && boxes.iter().all(|(_, box_level)| *box_level == 0);
            files.insert("system/fvSolution".into(), dicts::fv_solution(orthogonal));
            (dicts::SOLVER, speed_scale(domain, density), false)
        }
    };
    let frames = (*frames).max(1);
    let control = dicts::Control { solver, duration, frames, courant: options.max_courant, finest, speed, average_from };
    files.insert("system/controlDict".into(), dicts::control(&control));
    files.insert("system/decomposeParDict".into(), dicts::decompose(options.processes));
    files.insert("Allrun".into(), allrun(&step.name, solver, options.processes, set_fields));
    files.insert("case.foam".into(), String::new());

    Ok(Case { files, patches, background_cells: cells, solver, notes, density, gas, background, time_precision: dicts::time_precision(frames) })
}

/// The fluid's density and dynamic viscosity: both required, constant and Newtonian.
fn fluid_properties(model: &Model, domain: &FluidDomain) -> Result<(f64, f64)> {
    let name = &domain.material;
    let material = model.materials().get(name).ok_or_else(|| unsupported(format!("the fluid domain names no material `{name}`")))?;
    let density = material.density.ok_or_else(|| unsupported(format!("material `{name}` has no density")))?;
    if !(density.is_finite() && density > 0.0) {
        return Err(unsupported(format!("material `{name}`: a density that is not positive")));
    }
    let viscosity = match material.fluid {
        Some(Fluid::Newtonian { viscosity }) if viscosity.is_finite() && viscosity > 0.0 => viscosity,
        Some(Fluid::Newtonian { .. }) => return Err(unsupported(format!("material `{name}`: a viscosity that is not positive"))),
        Some(other) => return Err(unsupported(format!("material `{name}`: {} viscosity (only Newtonian fluids so far)", other.name()))),
        None => return Err(unsupported(format!("material `{name}` has no fluid behaviour"))),
    };
    if let Some(t) = &material.temperature_dependent
        && let Some((what, _)) = t.tables().find(|(what, _)| matches!(*what, "density" | "viscosity"))
    {
        return Err(unsupported(format!("material `{name}`: a {what} that follows the temperature (the case holds it constant)")));
    }
    Ok((density, viscosity))
}

/// The largest speed the boundaries and the initial state impose (m/s), measured as the Courant
/// number counts it: on a hex cell of side h in a uniform velocity, Co = Δt Σ_faces |φ| / (2V)
/// (doc.openfoam.com v2212, "CourantNo") is Δt (|u| + |v| + |w|) / h, so a velocity counts by its
/// L1 norm. Face and wall velocities and the initial velocity; and the speed `√(2 |Δp| / ρ)` a
/// pressure difference between faces drives (Bernoulli), normal to the face, where the norms agree.
/// The first time step is sized from it.
fn speed_scale(d: &FluidDomain, density: f64) -> f64 {
    let norm = |v: [f64; 3]| v[0].abs() + v[1].abs() + v[2].abs();
    let mut speed = norm(d.initial_velocity);
    let (mut p_lo, mut p_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for face in BoxFace::ALL {
        match *d.faces.get(face) {
            FaceCondition::Velocity { velocity, .. } | FaceCondition::MovingWall { velocity } => speed = speed.max(norm(velocity)),
            FaceCondition::Pressure { pressure, .. } => {
                p_lo = p_lo.min(pressure);
                p_hi = p_hi.max(pressure);
            }
            _ => {}
        }
    }
    for s in &d.surfaces {
        if let SurfaceCondition::MovingWall { velocity } = s.condition {
            speed = speed.max(norm(velocity));
        }
    }
    if p_hi > p_lo {
        speed = speed.max((2.0 * (p_hi - p_lo) / density).sqrt());
    }
    speed
}

fn fluid_step<'m>(model: &'m Model, name: Option<&str>) -> Result<&'m Step> {
    let mut fluid = model.steps().iter().filter(|s| matches!(s.procedure, Procedure::Fluid { .. }));
    match name {
        Some(n) => model.steps().iter().find(|s| s.name == n).ok_or_else(|| unsupported(format!("no step `{n}`"))),
        None => match (fluid.next(), fluid.next()) {
            (Some(s), None) => Ok(s),
            (None, _) => Err(unsupported("a model without a fluid step")),
            (Some(_), Some(_)) => Err(unsupported("several fluid steps without naming one")),
        },
    }
}

fn refuse_unsupported(d: &FluidDomain, notes: &mut Vec<String>) -> Result<()> {
    let refused = [
        (d.thermal.is_some(), "heat in the fluid"),
        (!d.porous.is_empty(), "porous zones"),
        (!d.fans.is_empty(), "fans"),
        (!d.baffles.is_empty(), "thin walls (baffles)"),
        (d.free_surface.is_some(), "a free surface"),
        (d.acceleration.iter().any(|a| *a != 0.0), "a body acceleration"),
    ];
    if let Some((_, what)) = refused.iter().find(|(on, _)| *on) {
        return Err(unsupported(format!("{what} (not mapped yet)")));
    }
    for face in BoxFace::ALL {
        match d.faces.get(face) {
            FaceCondition::Periodic {} => return Err(unsupported(format!("a periodic face ({}; not mapped yet)", face.name()))),
            FaceCondition::Velocity { profile: Profile::Parabolic, .. } => {
                return Err(unsupported(format!("a parabolic inlet profile ({}; not mapped yet)", face.name())));
            }
            FaceCondition::Pressure { outflow, .. } if !outflow.is_uniform() => {
                notes.push(format!("{}: the {outflow:?} outflow is written as an inflow/outflow pressure face", face.name()).to_lowercase());
            }
            _ => {}
        }
    }
    Ok(())
}

fn surface_condition(d: &FluidDomain, body: &str, face: &str) -> SurfaceCondition {
    d.surfaces
        .iter()
        .find(|s| {
            s.face.body == body && s.face.entity.name().is_some_and(|n| n == face || face.strip_prefix(n).is_some_and(|rest| rest.starts_with('#')))
        })
        .map_or(SurfaceCondition::Wall {}, |s| s.condition)
}

/// A name OpenFOAM takes as a word and a patch: letters, digits and `_`, starting with a letter
/// (`+` and `-` of face names such as `F1:+x` become `p` and `m`: `F1_px`).
fn sanitize(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    for c in name.chars() {
        match c {
            '+' => out.push('p'),
            '-' => out.push('m'),
            c if c.is_ascii_alphanumeric() || c == '_' => out.push(c),
            _ => out.push('_'),
        }
    }
    if !out.starts_with(|c: char| c.is_ascii_alphabetic()) {
        out.insert_str(0, "s_");
    }
    out
}

/// The body's name: `base`, else `base_2`, `base_3`, …, the first whose name and every
/// `<name>_<region>` patch are all unused; all of them are then used. A candidate fails only on a
/// name already used, and a used name blocks at most `regions + 1` candidates (as the candidate's
/// own name, or as `<candidate>_<region>` for one region each, which fixes the candidate), so
/// `used × (regions + 1) + 1` candidates always hold one that succeeds and bound the search.
fn claim(base: &str, regions: &[(String, String)], used: &mut BTreeSet<String>) -> Result<String> {
    let candidates = used.len().saturating_mul(regions.len().saturating_add(1)).saturating_add(1);
    for k in 1..=candidates {
        let name = if k == 1 { base.to_string() } else { format!("{base}_{k}") };
        let patches: Vec<String> = regions.iter().map(|(region, _)| format!("{name}_{region}")).collect();
        if !used.contains(&name) && !patches.iter().any(|p| used.contains(p)) {
            used.insert(name.clone());
            used.extend(patches);
            return Ok(name);
        }
    }
    Err(Error::Limit(format!("no unused name for `{base}` among {candidates} candidates")))
}

fn unique(base: &str, used: &mut BTreeMap<String, ()>) -> String {
    let mut name = base.to_string();
    let mut k = 2;
    while used.contains_key(&name) {
        name = format!("{base}_{k}");
        k += 1;
    }
    used.insert(name.clone(), ());
    name
}

/// The body as ASCII STL, one solid per face that has triangles (in face order), and the
/// solids' (region, face name) pairs. An imported surface has no face names: its faces are
/// `face0`, `face1`, … as the lattice solver calls them.
fn stl(o: &ObstacleSurface<'_>) -> Result<(String, Vec<(String, String)>)> {
    if o.tri_face.len() != o.triangles.len() {
        return Err(unsupported(format!("obstacle `{}`: a face index per triangle", o.body)));
    }
    let mut by_face = BTreeMap::<u32, Vec<usize>>::new();
    for (t, f) in o.tri_face.iter().enumerate() {
        by_face.entry(*f).or_default().push(t);
    }
    let mut used = BTreeMap::new();
    let mut text = String::new();
    let mut regions = Vec::new();
    for (f, tris) in by_face {
        let face = o.faces.get(f as usize).cloned().unwrap_or_else(|| format!("face{f}"));
        let region = unique(&sanitize(&face), &mut used);
        let _ = writeln!(text, "solid {region}");
        for t in tris {
            let tri = o.triangles.get(t).ok_or_else(|| unsupported("a triangle index out of range"))?;
            let p = |i: u32| {
                o.positions.get(i as usize).copied().ok_or_else(|| unsupported(format!("obstacle `{}`: a vertex index out of range", o.body)))
            };
            let (a, b, c) = (p(tri[0])?, p(tri[1])?, p(tri[2])?);
            if a.iter().chain(&b).chain(&c).any(|x| !x.is_finite()) {
                return Err(unsupported(format!("obstacle `{}`: a vertex that is not finite", o.body)));
            }
            let n = normal(a, b, c);
            let _ = writeln!(text, "  facet normal {} {} {}\n    outer loop", num(n[0]), num(n[1]), num(n[2]));
            for v in [a, b, c] {
                let _ = writeln!(text, "      vertex {} {} {}", num(v[0]), num(v[1]), num(v[2]));
            }
            let _ = writeln!(text, "    endloop\n  endfacet");
        }
        let _ = writeln!(text, "endsolid {region}");
        regions.push((region, face));
    }
    Ok((text, regions))
}

/// The unit normal of a triangle, or zero for a degenerate one (STL readers recompute normals
/// from the vertices; the normal written is informative).
fn normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if l > 0.0 { [n[0] / l, n[1] / l, n[2] / l] } else { [0.0; 3] }
}

/// Shortest round-trip text of a number (the same on every platform).
fn num(x: f64) -> String {
    if x == 0.0 { "0".into() } else { format!("{x:e}") }
}

fn vector(v: [f64; 3]) -> String {
    format!("({} {} {})", num(v[0]), num(v[1]), num(v[2]))
}

/// A dictionary file: the `FoamFile` header (user guide, "Basic input/output file format") and
/// the body.
fn dictionary(object: &str, body: &str) -> String {
    field_file("dictionary", object, body)
}

fn field_file(class: &str, object: &str, body: &str) -> String {
    format!(
        "FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {class};\n    object {object};\n}}\n// Written by FieldCraft for {TARGET}.\n\n{body}"
    )
}

/// One obstacle's STL surface and its regions.
struct Surface {
    name: String,
    regions: Vec<Region>,
}

/// One region (solid) of an obstacle's STL: its patch and the condition on it.
struct Region {
    patch: String,
    condition: SurfaceCondition,
}

/// Hex vertex order (OpenFOAM's block convention) and each box face's vertices, outward.
const BOX_FACES: [(BoxFace, [usize; 4]); 6] = [
    (BoxFace::XMin, [0, 4, 7, 3]),
    (BoxFace::XMax, [2, 6, 5, 1]),
    (BoxFace::YMin, [1, 5, 4, 0]),
    (BoxFace::YMax, [3, 7, 6, 2]),
    (BoxFace::ZMin, [0, 3, 2, 1]),
    (BoxFace::ZMax, [4, 5, 6, 7]),
];

/// The patch type of a box face (user guide §4.2.2, Table 4.1): a symmetry plane, a no-slip wall
/// (`wall`, whose cell distances wall-distance models measure from), or a generic patch (inlets,
/// outlets and slip faces, which FieldCraft's wall treatments do not measure from).
fn patch_type(c: &FaceCondition) -> &'static str {
    match c {
        FaceCondition::Symmetry {} => "symmetryPlane",
        FaceCondition::Wall {} | FaceCondition::MovingWall { .. } => "wall",
        _ => "patch",
    }
}

fn block_mesh(d: &FluidDomain, n: [u64; 3]) -> String {
    let (lo, hi) = (d.min, d.max);
    let corners = [
        [lo[0], lo[1], lo[2]],
        [hi[0], lo[1], lo[2]],
        [hi[0], hi[1], lo[2]],
        [lo[0], hi[1], lo[2]],
        [lo[0], lo[1], hi[2]],
        [hi[0], lo[1], hi[2]],
        [hi[0], hi[1], hi[2]],
        [lo[0], hi[1], hi[2]],
    ];
    let mut s = String::from("scale 1;\n\nvertices\n(\n");
    for c in corners {
        let _ = writeln!(s, "    {}", vector(c));
    }
    let _ = writeln!(s, ");\n\nblocks\n(\n    hex (0 1 2 3 4 5 6 7) ({} {} {}) simpleGrading (1 1 1)\n);\n\nboundary\n(", n[0], n[1], n[2]);
    for (face, v) in BOX_FACES {
        let _ = writeln!(
            s,
            "    {}\n    {{\n        type {};\n        faces (({} {} {} {}));\n    }}",
            face.name(),
            patch_type(d.faces.get(face)),
            v[0],
            v[1],
            v[2],
            v[3]
        );
    }
    s.push_str(");\n");
    dictionary("blockMeshDict", &s)
}

/// `0/U`, `0/p` and `0/nut`: each box face and each obstacle face (by patch name). Walls are
/// resolved (the wall models are not mapped), so `nut` takes no wall function there: zero
/// gradient puts the wall cell's eddy viscosity on the wall face (UG §A.4 "zeroGradient"), as the
/// lattice relaxes its wall cells with ν + ν_t; with Van Driest's damping that cell's ν_t is
/// already damped.
fn fields(d: &FluidDomain, density: f64, surfaces: &[Surface]) -> [(&'static str, String); 3] {
    let mut u = String::new();
    let mut p = String::new();
    let mut nut = String::new();
    for face in BoxFace::ALL {
        let (fu, fp, fn_) = match *d.faces.get(face) {
            FaceCondition::Velocity { velocity, .. } => {
                (format!("type fixedValue; value uniform {};", vector(velocity)), "type zeroGradient;".to_string(), "type zeroGradient;")
            }
            FaceCondition::Pressure { pressure, .. } => (
                "type pressureInletOutletVelocity; tangentialVelocity uniform (0 0 0); value uniform (0 0 0);".to_string(),
                format!("type fixedValue; value uniform {};", num(pressure / density)),
                "type zeroGradient;",
            ),
            FaceCondition::Slip {} => ("type slip;".into(), "type zeroGradient;".into(), "type zeroGradient;"),
            FaceCondition::Symmetry {} => ("type symmetryPlane;".into(), "type symmetryPlane;".into(), "type symmetryPlane;"),
            FaceCondition::Wall {} => ("type noSlip;".into(), "type zeroGradient;".into(), "type zeroGradient;"),
            FaceCondition::MovingWall { velocity } => {
                (format!("type fixedValue; value uniform {};", vector(velocity)), "type zeroGradient;".into(), "type zeroGradient;")
            }
            // Refused before the fields are written (`refuse_unsupported`).
            FaceCondition::Periodic {} => ("type cyclic;".into(), "type cyclic;".into(), "type cyclic;"),
        };
        let n = face.name();
        let _ = writeln!(u, "    {n} {{ {fu} }}");
        let _ = writeln!(p, "    {n} {{ {fp} }}");
        let _ = writeln!(nut, "    {n} {{ {fn_} }}");
    }
    for surface in surfaces {
        for r in &surface.regions {
            let fu = match r.condition {
                SurfaceCondition::Wall {} => "type noSlip;".to_string(),
                SurfaceCondition::MovingWall { velocity } => format!("type fixedValue; value uniform {};", vector(velocity)),
            };
            let _ = writeln!(u, "    {} {{ {fu} }}", r.patch);
            let _ = writeln!(p, "    {} {{ type zeroGradient; }}", r.patch);
            let _ = writeln!(nut, "    {} {{ type zeroGradient; }}", r.patch);
        }
    }
    let file = |class: &str, object: &str, dim: &str, internal: String, b: &str| {
        field_file(class, object, &format!("dimensions {dim};\n\ninternalField uniform {internal};\n\nboundaryField\n{{\n{b}}}\n"))
    };
    [
        ("U", file("volVectorField", "U", "[0 1 -1 0 0 0 0]", vector(d.initial_velocity), &u)),
        ("p", file("volScalarField", "p", "[0 2 -2 0 0 0 0]", "0".into(), &p)),
        ("nut", file("volScalarField", "nut", "[0 2 -1 0 0 0 0]", "0".into(), &nut)),
    ]
}

/// The script. The function name `writeCellCentres` of `postProcess -func` is a keyword of the
/// format; not in the documentation pages read (U4).
fn allrun(step: &str, solver: &str, processes: u32, set_fields: bool) -> String {
    let step: String = step.chars().filter(|c| !c.is_control()).collect();
    let mut s = format!(
        "#!/bin/sh\n# FieldCraft fluid step `{step}` as an OpenFOAM case: mesh, solve, reconstruct.\n\
# Run in this directory with OpenFOAM's environment. The ranks ({processes}) are those of\n\
# system/decomposeParDict: write the case again to change them.\n\
set -e\ncd \"$(dirname \"$0\")\"\n\
run() {{ log=\"log.$1\"; echo \"$*  > $log\"; \"$@\" > \"$log\" 2>&1; }}\n\
run blockMesh\nrun surfaceFeatureExtract\nrun snappyHexMesh -overwrite\nrun checkMesh\n\
run postProcess -func writeCellCentres -time 0\n"
    );
    if set_fields {
        s.push_str("run setFields\n");
    }
    if processes > 1 {
        let _ = write!(
            s,
            "rm -rf processor*\nrun decomposePar -force\n\
echo \"mpirun -np {processes} {solver} -parallel  > log.{solver}\"\nmpirun -np {processes} {solver} -parallel > log.{solver} 2>&1\n\
run reconstructPar\nrm -rf processor*\n",
        );
    } else {
        let _ = writeln!(s, "run {solver}");
    }
    s
}
