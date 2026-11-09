//! Equilibrium stiffness for linear perturbations (Bathe 2014, §6.8.2).
//!
//! Recompute the named equilibrium from current model inputs: output stress averages
//! are not restart state. Small-displacement elastic preload uses K + K_sigma.

use crate::analysis::{Setup, static_on};
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
    for e in &setup.mesh.structural {
        let curved = match e.kernel()? {
            crate::system::StructuralKernel::Mitc8(s) => !s.is_flat().map_err(|err| Error::element(e.label, err))?,
            crate::system::StructuralKernel::Mitc6(s) => !s.is_flat().map_err(|err| Error::element(e.label, err))?,
            _ => false,
        };
        if curved {
            return Err(Error::NotSupported("curved MITC8/MITC6 prestressed frequency/harmonic; only linear bifurcation verified".into()));
        }
    }
    crate::system::refuse_contact(setup.model, "linear static base_state; use a nonlinear equilibrium for contact")?;
    let preload = static_on(&Setup { model: setup.model, resolved: setup.resolved, mesh: setup.mesh, step: base }, limits)?;
    let kg = crate::buckle::initial_stress_stiffness(setup.mesh, graph, &preload)?;
    warnings.extend(preload.warnings);
    warnings
        .push(format!("linear perturbation about `{}`: reference-geometry elastic stiffness plus initial-stress stiffness; dead loads", base.name));
    let constraints = crate::tie::rows(setup.model, setup.resolved, setup.mesh, warnings)?.into_iter().map(FrozenRow::from_row).collect();
    Ok(Terminal { stiffness: k.add_scaled(1.0, &kg)?, constraints })
}
