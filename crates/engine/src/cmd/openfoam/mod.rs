//! OpenFOAM: write a fluid step as a case directory (`deck.openfoam`), an export for running the
//! case in OpenFOAM, which runs as its own program. The case is written from OpenFOAM's public
//! user guide (`fieldcraft_decks::openfoam`); FieldCraft's own fluid results are checked against
//! theory and published benchmarks (docs/verification.md), never against OpenFOAM's.

mod location;

use std::path::Path;

use fieldcraft_decks::openfoam::{CaseOptions, MeshOptions, ObstacleSurface, RefineBox, write_case};
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
    params: r#"{"path": "case_dir", "step"?: "Flow", "cell_factor"?: 1, "processes"?: 1, "max_courant"?: 1, "refine_boxes"?: [{"min": [x,y,z], "max": [x,y,z], "spacing"?: 0.1}], "location"?: [x,y,z]} — writes a fluid step (the only one when `step` is left out) as an OpenFOAM case directory, for running it in OpenFOAM: a blockMesh background of spacing × cell_factor (1, 2, 4, 8), snappyHexMesh around the obstacles (one patch per named face, refined back to the spacing; `refine_boxes` to theirs, rounded to the background over a power of two), ASCII fields, pimpleFoam with the step's turbulence model and a time step that adapts to max_courant (at most 1), decomposeParDict for `processes` ranks and an Allrun script (mesh, cell centres, solve, reconstruct). Every setting is written into the case. `location` is a point in the fluid, the side snappyHexMesh keeps: found when left out; when given, it is moved to the centre of the finest cell holding it (on no cell face) and refused inside an obstacle. The case directory must be new or empty. Lengths in the input units; the case is SI. Fails naming what the case cannot express yet (heat, porous zones, fans, baffles, free surface, periodic faces, parabolic inlets, non-Newtonian or temperature-dependent fluids) and any fluid property the material lacks. Nothing of the model changes"#,
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
    // What the case cannot express is refused before the obstacles are tessellated and searched.
    fieldcraft_decks::openfoam::check(&s.model, Some(&step_name))?;
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
    if let Some(c) = opt(p, "max_courant")? {
        options.max_courant = c;
    }
    let surfaces: Vec<ObstacleSurface<'_>> = meshes
        .iter()
        .map(|(body, positions, mesh, faces)| ObstacleSurface { body, positions, triangles: &mesh.triangles, tri_face: &mesh.tri_face, faces })
        .collect();
    let mut case = write_case(&s.model, Some(&step_name), &surfaces, &options)?;
    if given.is_some_and(|l| l != location) {
        case.notes.push("`location` moved to the centre of the finest cell holding it, which lies on no cell face".into());
    }
    let mut bytes = 0;
    for (rel, text) in &case.files {
        let path = root.join(rel);
        if let Some(parent) = path.parent() {
            fieldcraft_io::fs::create_dir_all(parent)?;
        }
        fieldcraft_io::fs::write_output(&path, text.as_bytes())?;
        bytes += text.len();
    }
    Ok(json!({"format": "openfoam", "path": dir, "step": step_name, "solver": case.solver, "files": case.files.keys().collect::<Vec<_>>(),
        "bytes": bytes, "background_cells": case.background_cells, "location": options.location, "patches": case.patches,
        "notes": case.notes, "run": format!("cd {} && sh Allrun", shell_quote(&dir))}))
}

#[cfg(test)]
mod tests;
