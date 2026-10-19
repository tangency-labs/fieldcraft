//! Pyramid transitions between a tet cell and the hex cells across its partitions (#439).
//!
//! The tet cell is meshed after the hex cells: its cut faces that face hex cells are replaced by
//! a pyramid on every hex element face there ([`fieldcraft_mesh::transition`]), the tetrahedra
//! fill the rest and keep the pyramids' sides, and the pyramids join the hex nodes to the tet
//! nodes. The mesh is conforming across the cut, so no tie is needed there.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_mesh::cdt::{CdtOptions, tetrahedralize};
use fieldcraft_mesh::surface::BoundarySurface;
use fieldcraft_mesh::tet::{QualityMesh, QualityOptions, quality_mesh};
use fieldcraft_mesh::transition::{Transition, apex_points, pyramid_surface};
use fieldcraft_model::{ElementKind, GeometryRef, Model};
use serde_json::{Value, json};

use super::super::hex::Input;
use super::super::mesh::{BodyTets, add_tets, analytic_surfaces};
use crate::{EngineError, Result};

/// A tet cell meshed against the pyramids on its hex neighbours' faces.
pub(super) struct Pyramids {
    pub(super) tets: BodyTets,
    transition: Transition,
    /// Each pyramid's base: the hex nodes of its face, oriented into the tet cell.
    bases: Vec<[u32; 4]>,
    /// The interfaces it makes (reply entries).
    pub(super) interfaces: Vec<Value>,
}

/// The base name of a cut face (`partition:P:a` of `partition:P:a#1`): its partition and side.
fn cut_of(name: &str) -> Option<(&str, usize)> {
    let base = name.split_once('#').map_or(name, |(b, _)| b);
    let (partition, side) = base.strip_prefix("partition:")?.rsplit_once(':')?;
    Some((partition, usize::from(side == "b")))
}

/// Mesh tet cell `input` with pyramids on the faces of the hex cells `hex` (already in `model`)
/// across its cut faces. `Ok(Err(why))` when there is no such face or the transition cannot be
/// built (the caller meshes the cell plainly and ties it); limits and a cancel are errors.
pub(super) fn mesh_with_pyramids(
    design: &fieldcraft_geometry::Design,
    model: &Model,
    p: &Value,
    input: &Input,
    hex: &[String],
    cdt: &CdtOptions,
    options: &QualityOptions,
    display: &fieldcraft_model::DisplaySystem,
) -> Result<std::result::Result<Pyramids, String>> {
    // The hex element faces across each of the cell's cut faces.
    let sides = super::sides(model, hex);
    let mut cut: Vec<u32> = Vec::new();
    let mut quads: Vec<[u32; 4]> = Vec::new();
    let mut interfaces = Vec::new();
    let mut seen: BTreeSet<(String, String)> = BTreeSet::new();
    for (i, name) in input.face_names.iter().enumerate() {
        let Some((partition, side)) = cut_of(name) else { continue };
        let Some(across) = sides.get(partition).and_then(|s| s.get(1 - side)) else { continue };
        if across.is_empty() {
            continue;
        }
        cut.push(i as u32);
        for other in across {
            if !seen.insert((partition.to_string(), other.body.clone())) {
                continue;
            }
            for f in &other.faces {
                let nodes = model.face_nodes(*f).unwrap_or_default();
                let Ok(q) = <[u32; 4]>::try_from(nodes.as_slice()) else {
                    return Ok(Err(format!("`{}` has quadratic or non-quadrilateral faces on the cut; pyramid transitions are linear", other.body)));
                };
                quads.push(q);
            }
            interfaces.push((partition.to_string(), other.body.clone(), other.faces.len()));
        }
    }
    if quads.is_empty() {
        return Ok(Err("no hex cell faces it across a partition".into()));
    }
    // The hex faces must cover the cut faces exactly (a cut shared with a tet cell, or a hex
    // cell larger than this one, is tied instead).
    let xyz = |n: u32| model.node(n).map(|x| x.xyz).unwrap_or([f64::NAN; 3]);
    let positions: Vec<[[f64; 3]; 4]> = quads.iter().map(|q| q.map(xyz)).collect();
    let quad_area: f64 = positions.iter().map(|q| 0.5 * norm(cross(sub(q[2], q[0]), sub(q[3], q[1])))).sum();
    let m = &input.mesh;
    let cut_area: f64 = m
        .triangles
        .iter()
        .zip(&m.tri_face)
        .filter(|(_, f)| cut.contains(f))
        .map(|(t, _)| {
            let c = t.map(|v| m.positions.get(v as usize).copied().unwrap_or([f64::NAN; 3]));
            0.5 * norm(cross(sub(c[1], c[0]), sub(c[2], c[0])))
        })
        .sum();
    if !((quad_area - cut_area).abs() <= 1e-3 * cut_area) {
        return Ok(Err(format!(
            "the hex cells' faces do not cover its cut faces (their areas differ by {:.3} %)",
            100.0 * (quad_area - cut_area).abs() / cut_area
        )));
    }
    let tag = u32::try_from(input.face_names.len()).map_err(|_| EngineError::Params("too many faces".into()))?;
    let welding = super::super::mesh::surface_options(p, m)?;
    let built = || -> fieldcraft_mesh::Result<(QualityMesh, Transition, usize, usize, f64)> {
        let surface = BoundarySurface::new(&m.positions, &m.triangles, &m.tri_face, &welding)?;
        let t = pyramid_surface(&surface, &cut, &positions, tag)?;
        // Already welded: nothing more merges.
        let exact = fieldcraft_mesh::surface::SurfaceOptions { weld_tolerance: 0.0, weld_distance: 0.0, ..welding };
        let region = BoundarySurface::new(&t.positions, &t.triangles, &t.tri_face, &exact)?;
        let tets = tetrahedralize(&region, cdt)?;
        let mut analytic = design.body(&input.name).map(analytic_surfaces).unwrap_or_default();
        analytic.resize(tag as usize, None);
        analytic.push(None);
        let o = QualityOptions { seeds: input.seeds.clone(), frozen_faces: vec![tag], ..options.clone() };
        let mesh = quality_mesh(&region, &tets, analytic, &o)?;
        Ok((mesh, t, tets.steiner.len(), tets.stats.cavities, region.volume()))
    };
    let (mesh, transition, steiner, cavities, volume) = match built() {
        Ok(x) => x,
        Err(fieldcraft_mesh::Error::Model(fieldcraft_model::Error::Cancelled)) => return Err(EngineError::Model(fieldcraft_model::Error::Cancelled)),
        Err(e) if e.kind() == Some(fieldcraft_mesh::Kind::Limit) => return Err(e.into()),
        Err(e) => return Ok(Err(crate::units::mesh_error_text(&e, display))),
    };
    if let Err(e) = apex_points(&mesh, &transition) {
        return Ok(Err(crate::units::mesh_error_text(&e, display)));
    }
    let tets = BodyTets::new(input.name.clone(), input.face_names.clone(), mesh, steiner, cavities, volume);
    let interfaces = interfaces
        .into_iter()
        .map(|(partition, other, faces)| json!({"partition": partition, "cells": [other, input.name], "joined": "pyramids", "pyramids": faces}))
        .collect();
    Ok(Ok(Pyramids { tets, transition, bases: quads, interfaces }))
}

/// Add the cell's tetrahedra and pyramids to `model`: the tetrahedra's nodes on the hex faces
/// are the hex nodes, the pyramids join them to the apexes.
pub(super) fn add(model: &mut Model, mut c: Pyramids) -> Result<(usize, usize)> {
    let tag = c.tets.face_names.len() as u32;
    let apexes = apex_points(&c.tets.mesh, &c.transition).map_err(|e| EngineError::Params(e.to_string()))?;
    // The pyramid sides are inside the part: no geometry face.
    c.tets.mesh.boundary.retain(|b| b.face != tag);
    let (first_node, _) = add_tets(model, &c.tets.name, &c.tets.face_names, &c.tets.mesh)?;
    // The tet mesh keeps the hex nodes' coordinates exactly: merge its copies into them.
    let hex: BTreeMap<[u64; 3], u32> = c.bases.iter().flatten().filter_map(|&n| model.node(n).map(|x| (x.xyz.map(f64::to_bits), n))).collect();
    let mut merge = BTreeMap::new();
    for (i, x) in c.tets.mesh.points.iter().enumerate() {
        if let Some(&h) = hex.get(&x.map(f64::to_bits)) {
            merge.insert(first_node + i as u32, h);
        }
    }
    let shared = merge.len();
    model.merge_nodes(&merge)?;
    let (_, first) = crate::cmd::fresh_labels(model, 0, c.bases.len())?;
    let body = model.geometry_tag(&GeometryRef::body(c.tets.name.as_str())?)?;
    for (k, (base, apex)) in c.bases.iter().zip(apexes).enumerate() {
        let label = first + k as u32;
        let mut nodes = base.to_vec();
        nodes.push(first_node + apex);
        model.set_element(label, ElementKind::Pyramid5, nodes)?;
        model.tag_element(label, body)?;
    }
    Ok((c.tets.mesh.elements.len() + c.bases.len(), c.tets.mesh.points.len() - shared))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
