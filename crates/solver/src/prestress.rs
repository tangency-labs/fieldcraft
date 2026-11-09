//! Equilibrium stiffness for linear perturbations (Bathe 2014, §6.8.2).
//!
//! Recompute the named equilibrium from current model inputs: output stress averages
//! are not restart state. Small-displacement elastic preload uses K + K_sigma.

use crate::analysis::{Setup, budget, reduction, static_on};
use crate::sparse::{Csr, NodeGraph};
use crate::{Error, Limits, Result};
use fieldcraft_model::Procedure;

/// Owned terminal constraint, valid only for the mesh layout used to build it.
#[derive(Debug, Clone)]
pub(crate) struct FrozenRow {
    pub name: String,
    pub terms: Vec<(usize, f64)>,
    pub pivot: usize,
    pub source: crate::constraints::Source,
}
impl FrozenRow {
    pub fn borrowed(&self) -> crate::constraints::Row<'_> {
        crate::constraints::Row { terms: self.terms.clone(), pivot: self.pivot, coupling: &self.name, source: self.source }
    }
    pub fn from_row(row: crate::constraints::Row<'_>) -> Self {
        Self { name: row.coupling.into(), terms: row.terms, pivot: row.pivot, source: row.source }
    }
}

pub(crate) struct Terminal {
    pub stiffness: Csr,
    pub constraints: Vec<FrozenRow>,
}

pub(crate) fn stiffness(setup: &Setup<'_>, graph: &NodeGraph, k: Csr, limits: &Limits, warnings: &mut Vec<String>) -> Result<Terminal> {
    let Some(base) = setup.model.perturbation_base(setup.step)? else { return Ok(Terminal { stiffness: k, constraints: Vec::new() }) };
    if setup.step.perturbation_follower_pressure && !matches!(base.procedure, Procedure::Static {}) {
        return Err(Error::NotSupported(
            "perturbation_follower_pressure requires a linear static shell base; nonlinear shell equilibrium is unavailable".into(),
        ));
    }
    if matches!(base.procedure, Procedure::NonlinearStatic(_)) {
        let (preload, terminal) = crate::nonlinear::terminal_tangent(setup.model, &base.name, limits)?;
        if preload.nodes != setup.mesh.labels || terminal.stiffness.n != k.n {
            return Err(Error::Invalid(format!("base_state `{}` changed the degree-of-freedom layout", base.name)));
        }
        warnings.extend(preload.warnings);
        warnings.push(format!(
            "linear perturbation about `{}`: converged consistent tangent, including stress and symmetric follower-load stiffness",
            base.name
        ));
        return Ok(terminal);
    }
    let mut curved = false;
    for e in &setup.mesh.structural {
        curved |= match e.kernel()? {
            crate::system::StructuralKernel::Mitc8(s) => !s.is_flat().map_err(|err| Error::element(e.label, err))?,
            crate::system::StructuralKernel::Mitc6(s) => !s.is_flat().map_err(|err| Error::element(e.label, err))?,
            _ => false,
        };
    }
    let shell_modal = curved || setup.step.perturbation_follower_pressure;
    if shell_modal {
        if setup.mesh.structural.iter().any(|e| e.kind == fieldcraft_model::ElementKind::Tri6) {
            return Err(Error::NotSupported("curved MITC6-a prestressed frequency/harmonic or MITC6-a follower-pressure modal tangent: preload-field and actual-pencil residual verification remain open".into()));
        }
        if !matches!(setup.step.procedure, Procedure::Frequency { .. }) {
            return Err(Error::NotSupported("curved MITC8 prestressed harmonic response; only modal perturbation verified".into()));
        }
        if !setup.mesh.elements.is_empty()
            || !setup.mesh.trusses.is_empty()
            || setup.mesh.structural.is_empty()
            || setup.mesh.structural.iter().any(|e| !matches!(e.kind, fieldcraft_model::ElementKind::Quad8 | fieldcraft_model::ElementKind::Tri6))
            || !setup.model.pretensions().is_empty()
        {
            return Err(Error::NotSupported("curved shell modal perturbation requires only MITC8 shells; mixed solids, beams, ring shells and pretension sections are not verified".into()));
        }
        if !setup.model.initial_fields().is_empty()
            || base.loads.iter().any(|name| matches!(setup.model.load(name), Some(fieldcraft_model::Load::Temperature { .. })))
        {
            return Err(Error::NotSupported(
                "curved shell modal perturbation with initial stress/temperature fields or thermal preload: shell eigenstrains are not implemented"
                    .into(),
            ));
        }
        if base.boundary_conditions != setup.step.boundary_conditions {
            return Err(Error::NotSupported("curved shell modal perturbation requires unchanged base-state supports".into()));
        }
    }
    crate::system::refuse_contact(setup.model, "linear static base_state; use a nonlinear equilibrium for contact")?;
    let fixed = crate::pretension::resolve_fixed(setup.model, &base.name, limits)?;
    let base = fixed.step(&base.name).ok_or_else(|| Error::Internal("resolved perturbation base".into()))?;
    let preload = static_on(&Setup { model: fixed.as_ref(), resolved: setup.resolved, mesh: setup.mesh, step: base }, limits)?;
    if setup.step.linear_base {
        warnings.extend(preload.warnings);
        warnings.push(format!("linear perturbation about `{}`: reference-geometry elastic stiffness, no initial-stress stiffness", base.name));
        let constraints = crate::tie::rows(setup.model, setup.resolved, setup.mesh, warnings)?.into_iter().map(FrozenRow::from_row).collect();
        return Ok(Terminal { stiffness: k, constraints });
    }
    let mut kg = crate::buckle::initial_stress_stiffness(setup.mesh, graph, &preload)?;
    if shell_modal {
        let (_, order) = crate::analysis::pattern(setup.mesh);
        let red = reduction(setup, false, &order, budget(limits, setup.mesh, graph, 2), warnings)?;
        if !red.is_selection() {
            return Err(Error::NotSupported("curved shell modal perturbation with multipoint or inclined constraints".into()));
        }
        if setup.step.perturbation_follower_pressure {
            let mut pressure = crate::buckle::pressure_stiffness(setup.model, setup.resolved, setup.mesh, base, graph)?;
            let [free_pressure] = red.reduce([&pressure])?;
            crate::buckle::conservative_pressure(&free_pressure)?;
            // Only after checking the selected free matrix; never mask a nonconservative load.
            pressure.mirror_lower();
            kg = kg.add_scaled(1.0, &pressure)?;
            warnings
                .push("quadratic-shell modal perturbation includes conservative reference-surface pressure stiffness; zero-offset faces only".into());
        }
        warnings.push("curved shell modes use reference geometry, reference mass and linear elastic preload; no finite-preload geometry or material-tangent update".into());
    }
    warnings.extend(preload.warnings);
    warnings.push(format!(
        "linear perturbation about `{}`: reference-geometry elastic stiffness plus initial-stress stiffness; {}",
        base.name,
        if setup.step.perturbation_follower_pressure { "conservative pressure load tangent included" } else { "dead loads" }
    ));
    let constraints = crate::tie::rows(setup.model, setup.resolved, setup.mesh, warnings)?.into_iter().map(FrozenRow::from_row).collect();
    Ok(Terminal { stiffness: k.add_scaled(1.0, &kg)?, constraints })
}
