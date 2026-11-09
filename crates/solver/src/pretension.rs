//! Constant-jump cut constraints by exact elimination (Cook et al., 2002, §9.2).
//! q = n·(u_first-u_second); the reference point's Ux is its work-conjugate DOF.
//! Coincident conforming cut nodes are paired deterministically. No penalty stiffness.
use crate::constraints::{Row, Source};
use crate::system::Mesh;
use crate::{Error, Result};
use fieldcraft_fem::DofKind;
use fieldcraft_model::{Model, PretensionControl, PretensionSection, Resolved, Step};
use std::collections::BTreeSet;

pub(crate) fn control_dof(mesh: &Mesh, section: &PretensionSection) -> Result<usize> {
    let point = mesh.point(&section.control).ok_or_else(|| Error::Invalid(format!("pretension control `{}` is missing", section.control)))?;
    mesh.dof(mesh.layout.point_slot(point), DofKind::Ux)
}

pub(crate) fn control_dofs(model: &Model, mesh: &Mesh) -> Result<BTreeSet<usize>> {
    model.pretensions().values().map(|p| control_dof(mesh, p)).collect()
}

pub(crate) fn rows<'a>(model: &'a Model, resolved: &Resolved, mesh: &Mesh) -> Result<Vec<Row<'a>>> {
    let mut rows = Vec::new();
    for (name, section) in model.pretensions() {
        if mesh.theory.is_some() {
            return Err(Error::NotSupported(format!("pretension `{name}` in 2D")));
        }
        // The scalar node must not also act as an ordinary load/coupling node.
        if model.reference_point_users(&section.control).len() != 1 {
            return Err(Error::Invalid(format!("pretension `{name}` needs a dedicated control point")));
        }
        let first = resolved.sets.get(&section.first).ok_or_else(|| Error::Invalid(format!("pretension `{name}` first set")))?;
        let second = resolved.sets.get(&section.second).ok_or_else(|| Error::Invalid(format!("pretension `{name}` second set")))?;
        if first.nodes.is_empty() || first.nodes.len() != second.nodes.len() || first.nodes.len() > 10_000 {
            return Err(Error::Invalid(format!("pretension `{name}` requires equally sized conforming cuts of 1..10000 nodes")));
        }
        let q = control_dof(mesh, section)?;
        let coordinates = |label: u32| -> Result<(usize, [f64; 3])> {
            let node = mesh.node(label).ok_or_else(|| Error::Invalid(format!("pretension `{name}`: cut node {label} has no structural DOFs")))?;
            Ok((node, *mesh.xyz.get(node).ok_or_else(|| Error::Internal("cut coordinate".into()))?))
        };
        let a: Vec<_> = first.nodes.iter().map(|n| coordinates(*n)).collect::<Result<_>>()?;
        let b: Vec<_> = second.nodes.iter().map(|n| coordinates(*n)).collect::<Result<_>>()?;
        let origin = a.first().map(|p| p.1).ok_or_else(|| Error::Internal("empty cut".into()))?;
        let extent = a.iter().map(|(_, x)| (0..3).map(|i| (x[i] - origin[i]).powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max);
        if extent <= 0.0 {
            return Err(Error::Invalid(format!("pretension `{name}` cut has zero extent")));
        }
        const RELATIVE_CUT_TOLERANCE: f64 = 1e-9; // unit-factors: dimensionless relative geometric tolerance, not a conversion
        let tolerance = RELATIVE_CUT_TOLERANCE * extent;
        let mut used = BTreeSet::new();
        for (node, x) in &a {
            let distance = (0..3).map(|i| (x[i] - origin[i]) * section.normal[i]).sum::<f64>();
            if distance.abs() > tolerance {
                return Err(Error::Invalid(format!("pretension `{name}` cut is not planar normal to its axis")));
            }
            let matches: Vec<_> = b.iter().filter(|(_, y)| (0..3).all(|i| (x[i] - y[i]).abs() <= tolerance)).collect();
            let [(other, _)] = matches.as_slice() else {
                return Err(Error::Invalid(format!("pretension `{name}` requires unique coincident cut-node pairs")));
            };
            if node == other || !used.insert(*other) {
                return Err(Error::Invalid(format!("pretension `{name}` cut nodes must be distinct and paired once")));
            }
            for (c, kind) in [DofKind::Ux, DofKind::Uy, DofKind::Uz].into_iter().enumerate() {
                let pivot = mesh.dof(*node, kind)?;
                let mut terms = vec![(pivot, 1.0), (mesh.dof(*other, kind)?, -1.0)];
                if section.normal[c] != 0.0 {
                    terms.push((q, -section.normal[c]));
                }
                rows.push(Row { terms, pivot, coupling: name, source: Source::Pretension });
            }
        }
    }
    Ok(rows)
}

pub(crate) fn prescribe(model: &Model, mesh: &Mesh, step: &Step, values: bool, out: &mut [Option<f64>]) -> Result<()> {
    for (name, section) in model.pretensions() {
        let value = match step.pretension.get(name) {
            Some(PretensionControl::Force { .. }) => None,
            Some(PretensionControl::Adjustment { value }) => Some(if values { *value } else { 0.0 }),
            None => Some(0.0), // undeformed cut or locked perturbation increment
        };
        let q = control_dof(mesh, section)?;
        let slot = out.get_mut(q).ok_or_else(|| Error::Internal("pretension prescribed DOF".into()))?;
        if slot.is_some() {
            return Err(Error::Invalid(format!("pretension `{name}` control has an ordinary boundary condition")));
        }
        *slot = value;
    }
    Ok(())
}

pub(crate) fn loads(model: &Model, mesh: &Mesh, step: &Step, out: &mut [f64]) -> Result<()> {
    for (name, control) in &step.pretension {
        if let PretensionControl::Force { value } = control {
            let section = model.pretensions().get(name).ok_or_else(|| Error::Invalid(format!("missing pretension `{name}`")))?;
            let q = control_dof(mesh, section)?;
            *out.get_mut(q).ok_or_else(|| Error::Internal("pretension force DOF".into()))? += value;
        }
    }
    Ok(())
}

/// Generalized cut force recovered from equilibrium, not copied from the load command.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SectionResult {
    pub name: String,
    pub adjustment: f64,
    pub section_force: f64,
    pub applied_force: f64,
}

pub(crate) fn section_dofs(model: &Model, mesh: &Mesh) -> Result<Vec<(String, usize)>> {
    model.pretensions().iter().map(|(name, p)| Ok((name.clone(), control_dof(mesh, p)?))).collect()
}

pub(crate) fn recover(
    sections: &[(String, usize)],
    step: &Step,
    red: &crate::constraints::Reduction,
    u: &[f64],
    internal: &[f64],
    external: &[f64],
) -> Result<Vec<SectionResult>> {
    if sections.is_empty() {
        return Ok(Vec::new());
    }
    let mut balance: Vec<f64> = internal.iter().zip(external).map(|(a, b)| a - b).collect();
    // Remove only the scalar cut loads from the external vector, then project by
    // virtual work onto the retained free/held adjustment coordinate.
    for (_, q) in sections {
        *balance.get_mut(*q).ok_or_else(|| Error::Internal("cut force balance DOF".into()))? += external.get(*q).copied().unwrap_or(0.);
    }
    let free = red.restrict(&balance);
    let held = red.reactions_of(&balance);
    sections
        .iter()
        .map(|(name, q)| {
            let section_force = match step.pretension.get(name) {
                Some(PretensionControl::Force { .. }) => {
                    let e = red
                        .eq_dof
                        .iter()
                        .position(|d| d == q)
                        .ok_or_else(|| Error::Invalid(format!("pretension `{name}` force control is constrained by another condition")))?;
                    free.get(e).copied().ok_or_else(|| Error::Internal("free cut force".into()))?
                }
                _ => held.get(*q).copied().ok_or_else(|| Error::Internal("held cut force".into()))?,
            };
            Ok(SectionResult {
                name: name.clone(),
                adjustment: u.get(*q).copied().ok_or_else(|| Error::Internal("cut adjustment".into()))?,
                section_force,
                applied_force: external.get(*q).copied().unwrap_or(0.),
            })
        })
        .collect()
}
