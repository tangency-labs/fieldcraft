//! Model assembly and result recovery for the harmonic matrix kernels.

use std::collections::BTreeMap;

use fieldcraft_fem::element::extrapolation;
use fieldcraft_model::harmonic::{Damping, Method};
use fieldcraft_model::{Model, OutputVariable, Procedure, Resolved, SectionKind, Step};
use serde::{Deserialize, Serialize};

use super::{Matrices, Modal};
use crate::analysis::{self, Setup, Stats};
use crate::complex::Vector;
use crate::sparse::{Csr, NodeGraph, add_block, add_diagonal, add_entry, add_scalar};
use crate::system::{self, Mesh, StructuralKernel};
use crate::{Error, Limits, MassKind, Result};

/// Physical displacement phasors and optional elastic stress phasors at one frequency.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    pub frequency: f64,
    pub real: Vec<[f64; 3]>,
    pub imaginary: Vec<[f64; 3]>,
    pub rotation_real: Vec<[f64; 3]>,
    pub rotation_imaginary: Vec<[f64; 3]>,
    pub points_real: Vec<[f64; 6]>,
    pub points_imaginary: Vec<[f64; 6]>,
    pub stress_real: Option<Vec<[f64; 6]>>,
    pub stress_imaginary: Option<Vec<[f64; 6]>>,
    pub relative_residual: f64,
    /// Algebraic L2 residual in internal SI equation scaling, without display conversion.
    /// Translational-only systems have force units; mixed rotational systems do not.
    pub residual_norm_si: f64,
    pub projected_residual: Option<f64>,
    /// Generalized reduced-system input power (equals physical input power for homogeneous supports).
    pub input_power: f64,
    pub dissipated_power: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HarmonicResult {
    /// Mass discretization used for both direct and modal harmonic operators.
    #[serde(default)]
    pub mass: MassKind,
    /// Min/max diagonal of the reduced physical operators, before modal projection.
    /// Generalized SI coordinates; mass is kg and stiffness N/m for translation-only systems.
    #[serde(default)]
    pub operator_diagonal_bounds: Option<[[f64; 2]; 2]>,
    #[serde(skip)]
    pub fatigue: Option<crate::fatigue::field::FatigueField>,
    pub step: String,
    pub nodes: Vec<u32>,
    pub reference_points: Vec<(String, u32)>,
    pub frames: Vec<Frame>,
    pub natural_frequencies: Vec<f64>,
    #[serde(default)]
    pub residual_vectors: usize,
    pub stats: Stats,
    pub warnings: Vec<String>,
}

fn scaled(values: &[f64], factor: f64) -> Vec<f64> {
    values.iter().map(|v| factor * v).collect()
}

fn at_frequency(error: Error, frequency: f64) -> Error {
    let context = |message| format!("harmonic frequency {frequency} Hz: {message}");
    match error {
        Error::Invalid(message) => Error::Invalid(context(message)),
        Error::NotSupported(message) => Error::NotSupported(context(message)),
        Error::Singular(message) => Error::Singular(context(message)),
        Error::Limit(message) => Error::Limit(context(message)),
        Error::NotConverged(message) => Error::NotConverged(context(message)),
        Error::Internal(message) => Error::Internal(context(message)),
        other => other,
    }
}

fn material_damping(model: &Model, resolved: &Resolved) -> BTreeMap<u32, Damping> {
    let mut out = BTreeMap::new();
    for section in model.sections().values() {
        let damping = model.materials().get(&section.material).and_then(|m| m.damping).unwrap_or_default();
        if damping == Damping::default() {
            continue;
        }
        if let Some(set) = resolved.sets.get(&section.set) {
            for &label in &set.elements {
                out.insert(label, damping);
            }
        }
    }
    out
}

/// Element damping is assembled on the existing graph; no averaging across materials.
fn damping_matrices(mesh: &Mesh, graph: &NodeGraph, materials: &BTreeMap<u32, Damping>, mass: MassKind) -> Result<(Option<Csr>, Option<Csr>)> {
    let viscous = materials.values().any(|d| d.alpha != 0.0 || d.beta != 0.0)
        || !mesh.discrete.dashpots().is_empty()
        || mesh.foundations.iter().any(|f| f.normal_damping != 0.0 || f.tangential_damping != 0.0);
    let loss = materials.values().any(|d| d.loss_factor != 0.0) || mesh.foundations.iter().any(|f| f.loss_factor != 0.0);
    let mut c = viscous.then(|| graph.zero_matrix(&mesh.layout));
    let mut g = loss.then(|| graph.zero_matrix(&mesh.layout));
    for e in &mesh.elements {
        let Some(d) = materials.get(&e.label) else { continue };
        let kernel = e.kernel()?;
        let stiffness = e.stiffness()?.data;
        let kinds = kernel.dofs().node;
        if let Some(c) = &mut c {
            add_block(c, &mesh.layout, &e.nodes, kinds, &scaled(&stiffness, d.beta))?;
            if d.alpha != 0.0 {
                let density = e.density.ok_or_else(|| Error::Invalid(format!("element {} has no density", e.label)))?;
                let (storage, values) = e.mass_matrix(mass, density)?;
                match storage {
                    MassKind::Consistent => add_scalar(c, &mesh.layout, &e.nodes, kinds, &scaled(&values, d.alpha))?,
                    MassKind::Lumped => add_diagonal(c, &mesh.layout, &e.nodes, kinds, &scaled(&values, d.alpha))?,
                }
            }
        }
        if let Some(g) = &mut g {
            add_block(g, &mesh.layout, &e.nodes, kinds, &scaled(&stiffness, d.loss_factor))?;
        }
    }
    for t in &mesh.trusses {
        let Some(d) = materials.get(&t.label) else { continue };
        let kernel = t.truss()?;
        let at = |err| Error::element(t.label, err);
        let kinds = mesh.translation_kinds();
        let stiffness = system::translation_block(&kernel.stiffness(t.young, t.area).map_err(at)?.data, t.nodes.len(), kinds);
        if let Some(c) = &mut c {
            add_block(c, &mesh.layout, &t.nodes, kinds, &scaled(&stiffness, d.beta))?;
            if d.alpha != 0.0 {
                let density = t.density.ok_or_else(|| Error::Invalid(format!("element {} has no density", t.label)))?;
                match mass {
                    MassKind::Consistent => {
                        add_scalar(c, &mesh.layout, &t.nodes, kinds, &scaled(&kernel.mass(density, t.area).map_err(at)?.data, d.alpha))?
                    }
                    MassKind::Lumped => {
                        add_diagonal(c, &mesh.layout, &t.nodes, kinds, &scaled(&kernel.lumped_mass(density, t.area).map_err(at)?, d.alpha))?
                    }
                }
            }
        }
        if let Some(g) = &mut g {
            add_block(g, &mesh.layout, &t.nodes, kinds, &scaled(&stiffness, d.loss_factor))?;
        }
    }
    for e in &mesh.structural {
        let Some(d) = materials.get(&e.label) else { continue };
        let kernel = e.kernel()?;
        let at = |err| Error::element(e.label, err);
        let stiffness = kernel.stiffness().map_err(at)?.data;
        if let Some(c) = &mut c {
            add_block(c, &mesh.layout, &e.nodes, e.kinds(), &scaled(&stiffness, d.beta))?;
            if d.alpha != 0.0 {
                let density = e.density.ok_or_else(|| Error::Invalid(format!("element {} has no density", e.label)))?;
                add_block(c, &mesh.layout, &e.nodes, e.kinds(), &scaled(&kernel.mass(density, mass).map_err(at)?.data, d.alpha))?;
            }
        }
        if let Some(g) = &mut g {
            add_block(g, &mesh.layout, &e.nodes, e.kinds(), &scaled(&stiffness, d.loss_factor))?;
        }
    }
    for f in &mesh.foundations {
        if let Some(c) = &mut c {
            f.add(mesh, c, f.normal_damping, f.tangential_damping)?;
        }
        if let Some(g) = &mut g {
            f.add(mesh, g, f.loss_factor * f.normal, f.loss_factor * f.tangential)?;
        }
    }
    if let Some(c) = &mut c {
        for spring in mesh.discrete.springs.iter().filter(|s| s.dashpot) {
            for &(i, ci) in &spring.terms {
                for &(j, cj) in &spring.terms {
                    add_entry(c, i, j, spring.value * ci * cj)?;
                }
            }
        }
        c.mirror_lower();
    }
    if let Some(g) = &mut g {
        g.mirror_lower();
    }
    Ok((c, g))
}

pub(crate) fn stresses(mesh: &Mesh, u: &[f64], displacement: &[[f64; 3]]) -> Result<Vec<[f64; 6]>> {
    let mut sum = vec![[0.0; 6]; mesh.labels.len()];
    let mut count = vec![0usize; mesh.labels.len()];
    let mut accumulate = |nodes: &[usize], values: &[[f64; 6]]| {
        for (node, stress) in nodes.iter().zip(values) {
            if let (Some(s), Some(n)) = (sum.get_mut(*node), count.get_mut(*node)) {
                for (s, v) in s.iter_mut().zip(stress) {
                    *s += v;
                }
                *n += 1;
            }
        }
    };
    for e in &mesh.elements {
        let kernel = e.kernel()?;
        let at = |err| Error::element(e.label, err);
        let ue: Vec<_> = e.nodes.iter().map(|n| displacement.get(*n).copied().unwrap_or_default()).collect();
        let points = e.stresses_from(&ue, None, None)?;
        let values = points.iter().map(|p| p.stress).collect::<Vec<_>>();
        let nodal = extrapolation(kernel.ty(), kernel.rule()).map_err(at)?.apply(&values);
        accumulate(&e.nodes, &nodal);
    }
    for e in &mesh.trusses {
        let kernel = e.truss()?;
        let at = |err| Error::element(e.label, err);
        let ue: Vec<_> = e.nodes.iter().map(|n| displacement.get(*n).copied().unwrap_or_default()).collect();
        let points = kernel.stresses_with(e.young, &ue, None).map_err(at)?;
        let values = points.iter().map(|p| p.stress).collect::<Vec<_>>();
        accumulate(&e.nodes, &kernel.extrapolate(&values).map_err(at)?);
    }
    for e in &mesh.structural {
        let dofs = mesh.layout.element_dofs(&e.nodes, e.kinds()).ok_or_else(|| Error::Internal("harmonic element dofs".into()))?;
        let ue: Vec<_> = dofs.iter().map(|d| u.get(*d).copied().unwrap_or(0.0)).collect();
        if let StructuralKernel::Axishell(shell) = e.kernel()? {
            accumulate(&e.nodes, &shell.stresses(&ue).map_err(|err| Error::element(e.label, err))?.top);
        }
        if let StructuralKernel::Shell(shell) = e.kernel()? {
            accumulate(&e.nodes, &shell.stresses(&ue).map_err(|err| Error::element(e.label, err))?.top);
        }
    }
    for (s, n) in sum.iter_mut().zip(count) {
        if n != 0 {
            for v in s {
                *v /= n as f64;
            }
        }
    }
    Ok(sum)
}

/// Shared physical load assembly for primal and differentiated harmonic systems.
pub(crate) fn load_vector(model: &Model, resolved: &Resolved, mesh: &Mesh, step: &Step, active: &[bool]) -> Result<Vector> {
    let Procedure::Harmonic(h) = &step.procedure else {
        return Err(Error::Internal("harmonic load step".into()));
    };
    let mut force = Vector::zeros(mesh.dofs());
    for name in &step.loads {
        let mut single = Step::new(step.name.clone(), step.procedure.clone());
        single.loads.insert(name.clone());
        let f = system::loads(model, resolved, mesh, &single, active)?;
        let factor = h.load_factors.get(name).copied().unwrap_or_default();
        let (re, im) = factor.components();
        for ((r, i), value) in force.real.iter_mut().zip(&mut force.imaginary).zip(f) {
            *r += re * value;
            *i += im * value;
        }
    }
    Ok(force)
}

pub(crate) fn damping_operators(model: &Model, resolved: &Resolved, mesh: &Mesh, graph: &NodeGraph) -> Result<(Option<Csr>, Option<Csr>)> {
    damping_operators_with_mass(model, resolved, mesh, graph, MassKind::Consistent)
}

pub(crate) fn damping_operators_with_mass(
    model: &Model,
    resolved: &Resolved,
    mesh: &Mesh,
    graph: &NodeGraph,
    mass: MassKind,
) -> Result<(Option<Csr>, Option<Csr>)> {
    damping_matrices(mesh, graph, &material_damping(model, resolved), mass)
}

/// Solve a harmonic model step with consistent mass and the original boundary
/// reduction. Existing spatial load definitions receive per-step phasor factors.
pub fn solve_harmonic(model: &Model, name: Option<&str>, limits: &Limits) -> Result<HarmonicResult> {
    solve_harmonic_with_mass(model, name, MassKind::Consistent, limits)
}

/// Solve using the same mass choice for inertia, material damping and modal reduction.
pub fn solve_harmonic_with_mass(model: &Model, name: Option<&str>, mass: MassKind, limits: &Limits) -> Result<HarmonicResult> {
    solve_harmonic_impl(model, name, mass, limits, None)
}

/// Reuse one exactly matched modal basis across calls, while retaining per-call assembly,
/// limits, response recovery and physical residual checks. Prestressed models never reuse it.
pub fn solve_harmonic_with_mass_cached(
    model: &Model,
    name: Option<&str>,
    mass: MassKind,
    limits: &Limits,
    cache: &mut super::ModalCache,
) -> Result<HarmonicResult> {
    solve_harmonic_impl(model, name, mass, limits, Some(cache))
}

fn solve_harmonic_impl(
    model: &Model,
    name: Option<&str>,
    mass: MassKind,
    limits: &Limits,
    mut cache: Option<&mut super::ModalCache>,
) -> Result<HarmonicResult> {
    fieldcraft_model::progress::check()?;
    if model.pretensions().values().any(|p| p.face.is_some()) {
        let resolved = model.resolve()?;
        let prepared = crate::pretension::prepare(model, &resolved, limits.max_dofs)?;
        let mut result = solve_harmonic_with_mass(prepared.model.as_ref(), name, mass, limits)?;
        prepared.project_harmonic(&mut result);
        return Ok(result);
    }
    let resolved = model.resolve()?;
    let step = model.select_step(name, |p| matches!(p, Procedure::Harmonic(_)), "harmonic")?;
    let Procedure::Harmonic(h) = &step.procedure else {
        return Err(Error::Internal("harmonic step".into()));
    };
    if step.base_state.is_none() {
        system::refuse_contact(model, "harmonic")?;
    }
    for name in &step.loads {
        if matches!(
            model.load(name),
            Some(
                fieldcraft_model::Load::Temperature { .. }
                    | fieldcraft_model::Load::SurfaceFlux { .. }
                    | fieldcraft_model::Load::Convection { .. }
                    | fieldcraft_model::Load::HeatSource { .. }
            )
        ) {
            return Err(Error::NotSupported(format!("harmonic thermal load `{name}`")));
        }
    }
    // A laminate has no single section material. Until ply-weighted damping is
    // assembled, refuse it explicitly instead of silently dropping ply damping.
    for (name, section) in model.sections() {
        if let SectionKind::CompositeShell { layup } = &section.kind {
            for ply in &layup.plies {
                if model.materials().get(&ply.material).and_then(|m| m.damping).is_some_and(|d| d != Damping::default()) {
                    return Err(Error::NotSupported(format!("section `{name}`: harmonic per-ply material damping; use analysis-step damping")));
                }
            }
        }
    }
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    if step.output.contains(&OutputVariable::Stress) && mesh.structural.iter().any(|e| matches!(e.section, system::Section::Beam { .. })) {
        return Err(Error::NotSupported("harmonic beam section stress; select displacement output for beam response".into()));
    }
    let frequencies = h.frequencies.points()?;
    let per_frame = mesh
        .labels
        .len()
        .checked_mul(if step.output.contains(&OutputVariable::Stress) { 24 } else { 12 })
        .and_then(|n| mesh.points.len().checked_mul(12).and_then(|p| n.checked_add(p)));
    if per_frame.and_then(|n| n.checked_mul(frequencies.len())).is_none_or(|n| n > limits.max_factor_entries) {
        return Err(Error::Limit("harmonic stored frames exceed the entry budget; request fewer frequencies or omit stress".into()));
    }
    let mut warnings: Vec<String> = resolved.warnings.iter().map(ToString::to_string).chain(mesh.warnings.iter().cloned()).collect();
    warnings.extend(system::inelastic_warnings(model));
    warnings.extend(fieldcraft_model::initial::left_out(model, step));
    let (graph, order) = analysis::pattern(&mesh);
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let assembly_timing = crate::timing::Span::new("harmonic_assembly_reduction");
    let assembled = analysis::assemble(&mesh, &graph, Some(mass))?;
    let terminal = crate::prestress::stiffness(&setup, &graph, assembled.k, limits, &mut warnings)?;
    let red = analysis::reduction_with_frozen(
        &setup,
        true,
        &order,
        analysis::budget(limits, &mesh, &graph, 4),
        &mut warnings,
        step.base_state.as_ref().map(|_| terminal.constraints.as_slice()),
    )?;
    if red.neq() == 0 {
        return Err(Error::Invalid("every harmonic degree of freedom is prescribed".into()));
    }
    let mut k = terminal.stiffness;
    let mut m = assembled.m.ok_or_else(|| Error::Internal("harmonic mass".into()))?;
    k.mirror_lower();
    m.mirror_lower();
    let (c, g) = damping_matrices(&mesh, &graph, &material_damping(model, &resolved), mass)?;
    let [mut kr, mut mr] = red.reduce([&k, &m])?;
    kr.mirror_lower();
    mr.mirror_lower();
    let reduce_optional = |a: &Option<Csr>| -> Result<Option<Csr>> {
        a.as_ref()
            .map(|a| {
                let [mut reduced] = red.reduce([a])?;
                reduced.mirror_lower();
                Ok(reduced)
            })
            .transpose()
    };
    let cr = reduce_optional(&c)?;
    let gr = reduce_optional(&g)?;
    let force = load_vector(model, &resolved, &mesh, step, &red.active)?;
    let full = (red.largest_prescribed() != 0.0).then_some((k, m, c, g));
    let bounds = |a: &Csr| a.diagonal().iter().fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], &v| [lo.min(v), hi.max(v)]);
    let operator_diagonal_bounds = Some([bounds(&kr), bounds(&mr)]);
    let matrices = Matrices { stiffness: &kr, mass: &mr, viscous: cr.as_ref(), loss: gr.as_ref() };
    let reduced_force = Vector { real: red.restrict(&force.real), imaginary: red.restrict(&force.imaginary) };
    drop(assembly_timing);
    let basis_timing = crate::timing::Span::new("harmonic_basis");
    let modal = match h.method {
        Method::Direct => None,
        Method::Modal { modes, residual_vectors, .. } => {
            if residual_vectors > 0 && full.is_some() {
                return Err(Error::NotSupported("harmonic residual vectors with nonzero prescribed motion".into()));
            }
            let reusable = model.pretensions().is_empty() && step.base_state.is_none() && full.is_none();
            Some(if let Some(cache) = cache.as_mut().filter(|_| reusable) {
                let context = super::cache::BasisContext {
                    mass,
                    sections: model.sections().iter().map(|(n, s)| (n.clone(), s.formulation)).collect(),
                    nodes: mesh.labels.clone(),
                    order: order.clone(),
                    constraints: red.cache_identity(),
                    scalar_projection: std::env::var("FIELDCRAFT_PROJECTION").is_ok_and(|v| v == "scalar"),
                    cached_dense: std::env::var("FIELDCRAFT_DENSE_UPDATE").is_ok_and(|v| v == "cached"),
                };
                cache.get(&matrices, modes as usize, &reduced_force, residual_vectors as usize, limits, &context)?
            } else {
                Modal::with_residual_vectors(&matrices, modes as usize, &reduced_force, residual_vectors as usize, limits)?
            })
        }
    };
    drop(basis_timing);
    let mut frames = Vec::with_capacity(frequencies.len());
    let mut stats = Stats { dofs: mesh.dofs(), equations: red.neq(), factor_entries: 0, lanczos_vectors: 0 };
    if let Some(m) = &modal {
        stats.factor_entries = m.factor_entries;
        stats.lanczos_vectors = m.lanczos_vectors;
    }
    let response_timing = crate::timing::Span::new("harmonic_response_recovery");
    let batch = modal.is_some() && full.is_none() && model.pretensions().is_empty() && step.base_state.is_none() && frequencies.len() > 1;
    let prepared = if let Some(modal) = modal.as_ref().filter(|_| batch) {
        let _timing = crate::timing::Span::new("harmonic_modal_coordinates");
        let Method::Modal { damping_ratios, .. } = &h.method else {
            return Err(Error::Internal("harmonic modal method".into()));
        };
        let ratios = if damping_ratios.len() > 1 { damping_ratios.get(..modal.eigenvalues().len()).unwrap_or_default() } else { damping_ratios };
        let mut projected = Vector::zeros(modal.basis_vectors().len());
        for (j, phi) in modal.basis_vectors().iter().enumerate() {
            projected.real[j] = super::dot(phi, &reduced_force.real);
            projected.imaginary[j] = super::dot(phi, &reduced_force.imaginary);
        }
        frequencies
            .iter()
            .map(|&frequency| {
                fieldcraft_model::progress::check()?;
                modal
                    .prepare_coordinates(2.0 * std::f64::consts::PI * frequency, h.damping, ratios, &reduced_force, limits, Some(&projected))
                    .map_err(|e| at_frequency(e, frequency))
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        Vec::new()
    };
    let mut recovered = std::collections::VecDeque::new();
    // At most four frequency columns and 256 MiB of temporary Phi Q buffers.
    let panel = (256 * 1024 * 1024usize / red.neq().saturating_mul(64).max(1)).clamp(1, 4);
    let mut stress_frames = if step.output.contains(&OutputVariable::Stress) && !prepared.is_empty() {
        let _timing = crate::timing::Span::new("harmonic_modal_stress");
        let modal = modal.as_ref().ok_or_else(|| Error::Internal("harmonic stress basis".into()))?;
        let mut fields = vec![(vec![[0.0; 6]; mesh.labels.len()], vec![[0.0; 6]; mesh.labels.len()]); frequencies.len()];
        // Stream one modal stress field at a time; do not retain a modes*nodes tensor.
        for (j, phi) in modal.basis_vectors().iter().enumerate() {
            fieldcraft_model::progress::check()?;
            let u = red.expand(phi, false);
            let (translation, _, _) = analysis::split(&u, &mesh);
            let stress = stresses(&mesh, &u, &translation)?;
            fieldcraft_parallel::for_each_mut(&mut fields, |f, (real, imaginary)| {
                for ((r, i), values) in real.iter_mut().zip(imaginary).zip(&stress) {
                    for k in 0..6 {
                        r[k] += values[k] * prepared[f].q.real[j];
                        i[k] += values[k] * prepared[f].q.imaginary[j];
                    }
                }
            });
        }
        fields.into_iter().collect::<std::collections::VecDeque<_>>()
    } else {
        std::collections::VecDeque::new()
    };
    for (frequency_index, frequency) in frequencies.into_iter().enumerate() {
        fieldcraft_model::progress::check()?;
        let omega = 2.0 * std::f64::consts::PI * frequency;
        let f = if let Some((k, m, c, g)) = &full {
            let (a, b) = Matrices { stiffness: k, mass: m, viscous: c.as_ref(), loss: g.as_ref() }.dynamic(omega, h.damping, limits)?;
            Vector { real: red.rhs(&force.real, &a), imaginary: red.rhs(&force.imaginary, &b) }
        } else {
            reduced_force.clone()
        };
        let (q, relative_residual, residual_norm, projected_residual, modal_power, loss_operator) = if let Some(modal) = &modal {
            let Method::Modal { damping_ratios, .. } = &h.method else {
                return Err(Error::Internal("harmonic modal method".into()));
            };
            let ratios = if damping_ratios.len() > 1 { damping_ratios.get(..modal.eigenvalues().len()).unwrap_or_default() } else { damping_ratios };
            let (out, loss) = if prepared.is_empty() {
                modal.solve_with_operator(omega, h.damping, ratios, &f, limits)
            } else {
                if recovered.is_empty() {
                    recovered = modal.recover_batch(&prepared[frequency_index..(frequency_index + panel).min(prepared.len())]);
                }
                let (u, d) = recovered.pop_front().ok_or_else(|| Error::Internal("harmonic response panel".into()))?;
                modal.finish_response(h.damping, &f, limits, &prepared[frequency_index], u, d)
            }
            .map_err(|e| at_frequency(e, frequency))?;
            let power = modal
                .eigenvalues()
                .iter()
                .enumerate()
                .map(|(j, lambda)| {
                    let ratio = ratios.get(if ratios.len() == 1 { 0 } else { j }).copied().unwrap_or(0.0);
                    let (real, imaginary) = (out.coordinates.real[j], out.coordinates.imaginary[j]);
                    let square = real * real + imaginary * imaginary;
                    omega * omega * ratio * lambda.max(0.0).sqrt() * square
                })
                .sum::<f64>();
            (out.displacement, out.full_relative_residual, out.full_residual_norm, Some(out.projected_residual.relative), power, Some(loss))
        } else {
            let out = matrices.direct(omega, h.damping, &f, limits).map_err(|e| at_frequency(e, frequency))?;
            stats.factor_entries = stats.factor_entries.max(out.factor_entries);
            (out.displacement, out.residual.relative, out.residual.norm, None, 0.0, None)
        };
        let input_power = 0.5
            * omega
            * f.real.iter().zip(&f.imaginary).zip(q.real.iter().zip(&q.imaginary)).map(|((fr, fi), (ur, ui))| fi * ur - fr * ui).sum::<f64>();
        let b = match loss_operator {
            Some(b) => b,
            None => matrices.dynamic(omega, h.damping, limits)?.1,
        };
        let mut br = vec![0.0; b.n];
        let mut bi = vec![0.0; b.n];
        b.mul(&q.real, &mut br);
        b.mul(&q.imaginary, &mut bi);
        let dissipated_power = 0.5 * omega * q.real.iter().zip(br).chain(q.imaginary.iter().zip(bi)).map(|(u, b)| u * b).sum::<f64>() + modal_power;
        let ur = red.expand(&q.real, true);
        let ui = red.expand(&q.imaginary, false);
        let (real, rotation_real, points_real) = analysis::split(&ur, &mesh);
        let (imaginary, rotation_imaginary, points_imaginary) = analysis::split(&ui, &mesh);
        let (stress_real, stress_imaginary) = if let Some((real, imaginary)) = stress_frames.pop_front() {
            (Some(real), Some(imaginary))
        } else {
            (
                step.output.contains(&OutputVariable::Stress).then(|| stresses(&mesh, &ur, &real)).transpose()?,
                step.output.contains(&OutputVariable::Stress).then(|| stresses(&mesh, &ui, &imaginary)).transpose()?,
            )
        };
        frames.push(Frame {
            frequency,
            real,
            imaginary,
            rotation_real,
            rotation_imaginary,
            points_real,
            points_imaginary,
            stress_real,
            stress_imaginary,
            relative_residual,
            residual_norm_si: residual_norm,
            projected_residual,
            input_power,
            dissipated_power,
        });
    }
    drop(response_timing);
    let natural_frequencies =
        modal.as_ref().map(|m| m.eigenvalues().iter().map(|l| l.max(0.0).sqrt() / (2.0 * std::f64::consts::PI)).collect()).unwrap_or_default();
    Ok(HarmonicResult {
        mass,
        operator_diagonal_bounds,
        fatigue: None,
        step: step.name.clone(),
        nodes: mesh.labels.clone(),
        reference_points: mesh.points.iter().map(|p| (p.name.clone(), p.label)).collect(),
        frames,
        natural_frequencies,
        residual_vectors: modal.as_ref().map_or(0, Modal::residual_vectors),
        stats,
        warnings,
    })
}

mod matrices;
pub use matrices::{MatrixExport, export_matrices, export_matrices_with_mass};

mod element_stress;
pub use element_stress::{ElementStressRecovery, StressLocation};
