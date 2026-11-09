//! Finite-element adapter shared by direct and modal time histories and spectra.
use super::Integrator;
use crate::analysis::{self, MassKind, Setup};
use crate::harmonic::{Matrices, Modal};
use crate::sparse::Csr;
use crate::system::{self, Mesh};
use crate::{Error, Limits, Result};
use fieldcraft_model::dynamics::{Combination, Integration, Method};
use fieldcraft_model::{Model, OutputVariable, Procedure, Step};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    /// None for a response-spectrum peak estimate: there is no physical time.
    pub time: Option<f64>,
    pub displacement: Vec<[f64; 3]>,
    pub rotation: Vec<[f64; 3]>,
    pub velocity: Vec<[f64; 3]>,
    /// Relative acceleration for transient, pseudo-acceleration peaks for spectrum.
    pub acceleration: Vec<[f64; 3]>,
    pub stress: Option<Vec<[f64; 6]>>,
    pub kinetic_energy: Option<f64>,
    pub strain_energy: Option<f64>,
}
/// Cumulative work in the integrated coordinates (SI joules).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Energy {
    pub damping_dissipation: f64,
    pub external_work: f64,
    /// E(t)-E(0)+damping_dissipation-external_work; includes algorithmic loss.
    pub energy_balance: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Response {
    /// Assembly choice used by this solve; preserved with the result.
    pub mass: MassKind,
    pub step: String,
    pub nodes: Vec<u32>,
    pub frames: Vec<Frame>,
    pub spectrum: bool,
    pub natural_frequencies: Vec<f64>,
    pub warnings: Vec<String>,
    /// One entry per stored frame; absent in legacy files and spectra.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub energy: Option<Vec<Energy>>,
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn project(a: &Csr, basis: &[Vec<f64>]) -> Csr {
    let n = basis.len();
    let mut val = vec![0.0; n * n];
    let mut av = vec![0.0; a.n];
    for j in 0..n {
        a.mul(&basis[j], &mut av);
        for i in j..n {
            let v = dot(&basis[i], &av);
            val[i * n + j] = v;
            val[j * n + i] = v;
        }
    }
    Csr { n, ptr: (0..=n).map(|i| i * n).collect(), idx: (0..n).flat_map(|_| 0..n).collect(), val }
}
fn project_force(f: &[f64], basis: &[Vec<f64>]) -> Vec<f64> {
    basis.iter().map(|v| dot(v, f)).collect()
}
fn expand(q: &[f64], basis: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    for (value, mode) in q.iter().zip(basis) {
        for (x, p) in out.iter_mut().zip(mode) {
            *x += value * p;
        }
    }
    out
}
pub(super) fn initial(values: &BTreeMap<u32, [f64; 3]>, mesh: &Mesh, red: &crate::constraints::Reduction) -> Result<Vec<f64>> {
    let mut full = vec![0.0; mesh.dofs()];
    for (label, v) in values {
        let slot = mesh.node(*label).ok_or_else(|| Error::Invalid(format!("initial node {label} not present")))?;
        for (c, x) in v.iter().enumerate() {
            let kind = fieldcraft_fem::DofKind::translation(c).ok_or_else(|| Error::Internal("translation direction".into()))?;
            match mesh.layout.dof(slot, kind).and_then(|d| full.get_mut(d)) {
                Some(target) => *target = *x,
                None if *x == 0.0 => {}
                _ => return Err(Error::Invalid(format!("initial node {label} has no requested degree of freedom"))),
            }
        }
    }
    let q: Vec<_> = red.eq_dof.iter().map(|d| full.get(*d).copied().unwrap_or(0.0)).collect();
    let reconstructed = red.expand(&q, false);
    let scale = full.iter().map(|x| x.abs()).fold(0.0, f64::max);
    if full.iter().zip(reconstructed).any(|(a, b)| (a - b).abs() > 1e-10 * scale) {
        return Err(Error::Invalid("initial state violates the active constraints".into()));
    }
    Ok(q)
}
fn translation(mesh: &Mesh, direction: [f64; 3]) -> Vec<f64> {
    let mut out = vec![0.0; mesh.dofs()];
    for slot in 0..mesh.labels.len() {
        for (c, value) in direction.iter().enumerate() {
            if let Some(d) = fieldcraft_fem::DofKind::translation(c).and_then(|k| mesh.layout.dof(slot, k))
                && let Some(x) = out.get_mut(d)
            {
                *x = *value;
            }
        }
    }
    out
}
fn energy(a: &Csr, q: &[f64]) -> f64 {
    let mut v = vec![0.0; a.n];
    a.mul(q, &mut v);
    0.5 * dot(q, &v)
}

pub fn solve(model: &Model, name: Option<&str>, limits: &Limits, spectrum: bool) -> Result<Response> {
    solve_with_mass(model, name, MassKind::Consistent, limits, spectrum)
}

/// Solve with one mass choice for assembly, damping, modes, initial balance and energy.
pub fn solve_with_mass(model: &Model, name: Option<&str>, mass: MassKind, limits: &Limits, spectrum: bool) -> Result<Response> {
    let resolved = model.resolve()?;
    let wanted: fn(&Procedure) -> bool =
        if spectrum { |p| matches!(p, Procedure::ResponseSpectrum(_)) } else { |p| matches!(p, Procedure::Transient(_)) };
    let step = model.select_step(name, wanted, if spectrum { "response_spectrum" } else { "transient" })?;
    step.procedure.validate()?;
    system::refuse_contact(model, "linear dynamics")?;
    // Reject data with no implemented linear time-domain meaning rather than drop it.
    if !model.initial_fields().is_empty() {
        return Err(Error::NotSupported("linear dynamics initial fields; use transient nodal initial states".into()));
    }
    for (name, m) in model.materials() {
        if m.damping.is_some_and(|d| d.loss_factor != 0.0) {
            return Err(Error::NotSupported(format!("material `{name}` structural loss in time/spectrum response")));
        }
    }
    for (name, f) in model.foundations() {
        if f.loss_factor != 0.0 {
            return Err(Error::NotSupported(format!("foundation `{name}` structural loss in time/spectrum response")));
        }
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
            return Err(Error::NotSupported(format!("linear dynamics thermal load `{name}`")));
        }
    }
    for (name, s) in model.sections() {
        if let fieldcraft_model::SectionKind::CompositeShell { layup } = &s.kind
            && layup.plies.iter().any(|p| model.materials().get(&p.material).and_then(|m| m.damping).is_some_and(|d| d != Default::default()))
        {
            return Err(Error::NotSupported(format!("section `{name}` per-ply damping in linear dynamics")));
        }
    }
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    if model.sections().values().any(|s| s.formulation.axisymmetric_shell) {
        return Err(Error::NotSupported("axisymmetric conical-shell transient/response spectrum".into()));
    }
    if model.sections().values().any(|s| s.formulation.solid_shell) {
        return Err(Error::NotSupported(
            "solid-shell transient and response-spectrum procedures are not validated; use frequency or harmonic steps".into(),
        ));
    }
    if mesh.structural.iter().any(|e| matches!(e.kind, fieldcraft_model::ElementKind::Quad8 | fieldcraft_model::ElementKind::Tri6)) {
        return Err(Error::NotSupported(
            "MITC8/MITC6-a transient and response-spectrum procedures are not validated; use frequency or harmonic steps".into(),
        ));
    }
    if matches!(&step.procedure, Procedure::Transient(t) if t.adaptive.is_some())
        && (0..mesh.dofs())
            .any(|d| !matches!(mesh.kind(d), Some(fieldcraft_fem::DofKind::Ux | fieldcraft_fem::DofKind::Uy | fieldcraft_fem::DofKind::Uz)))
    {
        return Err(Error::NotSupported(
            "adaptive transient requires translational degrees of freedom; rotation/pressure/thermal coordinates are unsupported".into(),
        ));
    }
    if !mesh.points.is_empty() {
        return Err(Error::NotSupported("linear dynamics reference-point output".into()));
    }
    if step.output.iter().any(|o| !matches!(o, OutputVariable::Displacement | OutputVariable::Stress)) {
        return Err(Error::NotSupported("linear dynamics output other than displacement/stress".into()));
    }
    let stress = step.output.contains(&OutputVariable::Stress);
    if stress && mesh.structural.iter().any(|e| matches!(e.section, system::Section::Beam { .. })) {
        return Err(Error::NotSupported("linear dynamics beam section stress".into()));
    }
    let frames = match &step.procedure {
        Procedure::Transient(t) => t.increments / t.output_every + 2,
        _ => 1,
    };
    if mesh
        .labels
        .len()
        .checked_mul(frames as usize)
        .and_then(|n| n.checked_mul(if stress { 18 } else { 12 }))
        .is_none_or(|n| n > limits.max_factor_entries.min(50_000_000))
    {
        return Err(Error::Limit("dynamic output budget exceeded; increase output_every or omit stress".into()));
    }
    let mut warnings: Vec<String> = resolved.warnings.iter().map(ToString::to_string).chain(mesh.warnings.iter().cloned()).collect();
    warnings.extend(system::inelastic_warnings(model));
    let (graph, order) = analysis::pattern(&mesh);
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let red = analysis::reduction(&setup, true, &order, analysis::budget(limits, &mesh, &graph, 4), &mut warnings)?;
    if red.largest_prescribed() != 0.0 {
        return Err(Error::NotSupported("linear dynamics nonzero prescribed motion; use base_acceleration".into()));
    }
    let assembled = analysis::assemble(&mesh, &graph, Some(mass))?;
    let mut k = assembled.k;
    let mut m = assembled.m.ok_or_else(|| Error::Internal("dynamic mass".into()))?;
    k.mirror_lower();
    m.mirror_lower();
    let (viscous, loss) = crate::harmonic::damping_operators_with_mass(model, &resolved, &mesh, &graph, mass)?;
    if loss.as_ref().is_some_and(|g| g.val.iter().any(|x| *x != 0.0)) {
        return Err(Error::NotSupported("structural loss in linear dynamics".into()));
    }
    let mut c = viscous.unwrap_or_else(|| Csr { val: vec![0.0; k.val.len()], ..k.clone() });
    if let Procedure::Transient(t) = &step.procedure {
        c = c.add_scaled(t.damping.alpha, &m)?.add_scaled(t.damping.beta, &k)?;
    }
    let [mut kr, mut mr, mut cr] = red.reduce([&k, &m, &c])?;
    kr.mirror_lower();
    mr.mirror_lower();
    cr.mirror_lower();
    let mut output = Response {
        mass,
        step: step.name.clone(),
        nodes: mesh.labels.clone(),
        frames: Vec::new(),
        spectrum,
        natural_frequencies: Vec::new(),
        warnings,
        energy: None,
    };
    if let Procedure::ResponseSpectrum(spec) = &step.procedure {
        if !step.loads.is_empty() || cr.val.iter().any(|c| *c != 0.0) {
            return Err(Error::NotSupported(
                "response spectrum requires no spatial loads or assembled damping; damping_ratio defines the input spectrum".into(),
            ));
        }
        let modal = Modal::new(&Matrices { stiffness: &kr, mass: &mr, viscous: None, loss: None }, spec.modes as usize, limits)?;
        let modes = modal.shapes();
        if modes
            .len()
            .checked_mul(mesh.labels.len())
            .and_then(|n| n.checked_mul(if stress { 15 } else { 9 }))
            .is_none_or(|n| n > limits.max_factor_entries)
        {
            return Err(Error::Limit("spectrum modal recovery exceeds entry budget".into()));
        }
        for &lambda in modal.eigenvalues() {
            if lambda <= 0.0 {
                return Err(Error::Invalid("spectrum rigid-body mode has no spectral displacement".into()));
            }
            output.natural_frequencies.push(libm::sqrt(lambda) / (2.0 * std::f64::consts::PI));
        }
        let damping: Vec<_> = (0..modes.len()).map(|i| spec.mode_damping(i)).collect();
        let ordinates: Vec<_> =
            output.natural_frequencies.iter().zip(&damping).map(|(f, z)| spec.ordinate(*f, *z)).collect::<std::result::Result<_, _>>()?;
        let mut direction_u = Vec::new();
        let mut direction_r = Vec::new();
        let mut direction_a = Vec::new();
        let mut direction_s = Vec::new();
        for excitation in spec.excitations() {
            fieldcraft_model::progress::check()?;
            let mut mrbase = vec![0.0; m.n];
            m.mul(&translation(&mesh, excitation.direction), &mut mrbase);
            let participation = project_force(&red.restrict(&mrbase), modes);
            let mut uf = Vec::new();
            let mut rf = Vec::new();
            let mut af = Vec::new();
            let mut sf = Vec::new();
            for (i, (&lambda, shape)) in modal.eigenvalues().iter().zip(modes).enumerate() {
                let q = participation.get(i).copied().unwrap_or(0.0) * ordinates.get(i).copied().unwrap_or(0.0) * excitation.scale / lambda;
                let full = red.expand(&shape.iter().map(|x| x * q).collect::<Vec<_>>(), false);
                let (u, r, _) = analysis::split(&full, &mesh);
                af.push(u.iter().map(|v| v.map(|x| x * lambda)).collect::<Vec<_>>());
                if stress {
                    sf.push(crate::harmonic::stresses(&mesh, &full, &u)?);
                }
                uf.push(u);
                rf.push(r);
            }
            let combine = |values: &[f64]| match spec.combination {
                Combination::Srss => crate::spectrum::srss(values),
                Combination::Cqc => crate::spectrum::cqc(values, &output.natural_frequencies, &damping),
            };
            direction_u.push(combine_fields(&uf, mesh.labels.len(), &combine)?);
            direction_r.push(combine_fields(&rf, mesh.labels.len(), &combine)?);
            direction_a.push(combine_fields(&af, mesh.labels.len(), &combine)?);
            if stress {
                direction_s.push(combine_fields(&sf, mesh.labels.len(), &combine)?);
            }
        }
        let combine_direction = |peaks: &[f64]| match spec.direction_combination {
            fieldcraft_model::dynamics::DirectionCombination::Srss => crate::spectrum::srss(peaks),
            fieldcraft_model::dynamics::DirectionCombination::Percent1004040 => {
                let mut xyz = [0.0; 3];
                for (x, p) in xyz.iter_mut().zip(peaks) {
                    *x = *p;
                }
                crate::spectrum::directional_100_40_40(xyz)
            }
        };
        output.frames.push(Frame {
            time: None,
            displacement: combine_fields(&direction_u, mesh.labels.len(), &combine_direction)?,
            rotation: combine_fields(&direction_r, mesh.labels.len(), &combine_direction)?,
            velocity: Vec::new(),
            acceleration: combine_fields(&direction_a, mesh.labels.len(), &combine_direction)?,
            stress: if stress { Some(combine_fields(&direction_s, mesh.labels.len(), &combine_direction)?) } else { None },
            kinetic_energy: None,
            strain_energy: None,
        });
        output.warnings.push(
            "Spectrum component peaks are not simultaneous fields; acceleration is pseudo-acceleration, not true absolute acceleration.".into(),
        );
    } else if let Procedure::Transient(t) = &step.procedure {
        for name in t.load_amplitudes.keys() {
            if !step.loads.contains(name) {
                return Err(Error::Invalid(format!("amplitude for inactive load `{name}`")));
            }
        }
        if step.loads.len().checked_mul(kr.n).is_none_or(|n| n > limits.max_factor_entries) {
            return Err(Error::Limit("dynamic load vectors exceed the entry budget".into()));
        }
        let mut load_vectors = Vec::new();
        for name in &step.loads {
            let mut one = Step::new(step.name.clone(), step.procedure.clone());
            one.loads.insert(name.clone());
            let f = red.restrict(&system::loads(model, &resolved, &mesh, &one, &red.active)?);
            load_vectors.push((f, t.load_amplitudes.get(name).map(|a| a.evaluator()).transpose()?));
        }
        let method = match t.integration {
            Integration::AverageAcceleration => super::Integration::AverageAcceleration,
            Integration::LinearAcceleration => super::Integration::LinearAcceleration,
            Integration::Hht { alpha } => super::Integration::Hht { alpha },
        };
        let base = if let Some(b) = &t.base_acceleration {
            let mut f = vec![0.0; m.n];
            m.mul(&translation(&mesh, b.direction), &mut f);
            Some((red.restrict(&f).iter().map(|x| -b.magnitude * x).collect::<Vec<_>>(), b.amplitude.evaluator()?))
        } else {
            None
        };
        // Enforce prescribed base positions with the same Newmark kinematics as
        // the free coordinates; corners in tables are resolved at the fixed dt.
        let mut motion = if let Some(b) = &t.base_displacement {
            let mut f = vec![0.0; m.n];
            m.mul(&translation(&mesh, b.direction), &mut f);
            let curve = b.amplitude.evaluator()?;
            let position = b.magnitude * curve.value(0.0)?;
            if !position.is_finite() {
                return Err(Error::Invalid("base displacement overflow".into()));
            }
            Some((red.restrict(&f), curve, b.magnitude, position, b.initial_velocity, b.initial_acceleration))
        } else {
            None
        };
        let (beta, gamma, _) = method.parameters()?;
        let mut force = |time: f64| -> Result<Vec<f64>> {
            let mut out = vec![0.0; kr.n];
            for (f, a) in &load_vectors {
                let value = match a {
                    Some(a) => a.value(time)?,
                    None => 1.0,
                };
                for (o, x) in out.iter_mut().zip(f) {
                    *o += value * x;
                }
            }
            if let Some((f, a)) = &base {
                let value = a.value(time)?;
                for (o, x) in out.iter_mut().zip(f) {
                    *o += value * x;
                }
            }
            if let Some((mr, curve, magnitude, position, velocity, acceleration)) = &mut motion {
                if time > 0.0 {
                    let target = *magnitude * curve.value(time)?;
                    let next_a = (target - *position - t.time_step * *velocity - t.time_step * t.time_step * (0.5 - beta) * *acceleration)
                        / (beta * t.time_step * t.time_step);
                    let next_v = *velocity + t.time_step * ((1.0 - gamma) * *acceleration + gamma * next_a);
                    if [target, next_v, next_a].iter().any(|x| !x.is_finite()) {
                        return Err(Error::Invalid("prescribed base motion overflow".into()));
                    }
                    *position = target;
                    *velocity = next_v;
                    *acceleration = next_a;
                }
                for (o, x) in out.iter_mut().zip(mr) {
                    *o -= *acceleration * *x;
                }
            }
            Ok(out)
        };
        let u = initial(&t.initial_displacement, &mesh, &red)?;
        let v = initial(&t.initial_velocity, &mesh, &red)?;

        let modal = if let Method::Modal { modes, .. } = t.method {
            Some(Modal::new(&Matrices { stiffness: &kr, mass: &mr, viscous: Some(&cr), loss: None }, modes as usize, limits)?)
        } else {
            None
        };
        let (ki, mi, ci, ui, vi) = if let Some(modal) = &modal {
            let basis = modal.shapes();
            if basis.len().checked_mul(basis.len()).and_then(|n| n.checked_mul(3)).is_none_or(|n| n > limits.max_factor_entries) {
                return Err(Error::Limit("projected dynamic operators exceed the entry budget".into()));
            }
            let mut ci = project(&cr, basis);
            if let Method::Modal { damping_ratios, .. } = &t.method {
                for (i, lambda) in modal.eigenvalues().iter().enumerate() {
                    let z = damping_ratios.get(if damping_ratios.len() == 1 { 0 } else { i }).copied().unwrap_or(0.0);
                    if let Some(c) = ci.val.get_mut(i * ci.n + i) {
                        *c += 2.0 * z * libm::sqrt(lambda.max(0.0));
                    }
                }
            }
            output.natural_frequencies = modal.eigenvalues().iter().map(|l| libm::sqrt(l.max(0.0)) / (2.0 * std::f64::consts::PI)).collect();
            let mut mu = vec![0.0; mr.n];
            let mut mv = vec![0.0; mr.n];
            mr.mul(&u, &mut mu);
            mr.mul(&v, &mut mv);
            output.warnings.push("Modal truncation projects initial states and loads; increase retained modes to check convergence.".into());
            (project(&kr, basis), project(&mr, basis), ci, project_force(&mu, basis), project_force(&mv, basis))
        } else {
            (kr.clone(), mr.clone(), cr.clone(), u, v)
        };
        if matches!(t.integration, Integration::LinearAcceleration) {
            output
                .warnings
                .push("Linear acceleration is conditionally stable: dt * highest retained/system angular frequency must not exceed sqrt(12).".into());
        }
        let project_f = |f: Vec<f64>| match &modal {
            Some(m) => project_force(&f, m.shapes()),
            None => f,
        };
        let mut f = project_f(force(0.0)?);
        let mut state = Integrator::initial_state(&ki, &mi, &ci, ui, vi, &f, limits)?;
        let integrator = if t.adaptive.is_none() { Some(Integrator::new(&ki, &mi, &ci, t.time_step, method, limits)?) } else { None };
        let mut controller = super::adaptive::Controller::new(t.time_step);
        let mut knots = Vec::new();
        if t.adaptive.is_some() {
            for curve in t.load_amplitudes.values().chain(t.base_acceleration.iter().map(|b| &b.amplitude)) {
                match curve {
                    fieldcraft_model::amplitude::Amplitude::Tabular { points } | fieldcraft_model::amplitude::Amplitude::SmoothStep { points } => {
                        knots.extend(points.iter().map(|p| p[0]))
                    }
                    fieldcraft_model::amplitude::Amplitude::Periodic { origin, .. } => knots.push(*origin),
                }
            }
            knots.sort_by(f64::total_cmp);
            knots.dedup();
        }
        let initial_energy = energy(&mi, &state.velocity) + energy(&ki, &state.displacement);
        let mut damping_dissipation = 0.0;
        let mut external_work = 0.0;
        let mut energies = Vec::new();
        for i in 0..=t.increments {
            fieldcraft_model::progress::check()?;
            let time = i as f64 * t.time_step;
            if i > 0
                && let Some(settings) = &t.adaptive
            {
                let (d, w) = controller.advance_to(
                    &mut state,
                    &mut f,
                    (i - 1) as f64 * t.time_step,
                    time,
                    &ki,
                    &mi,
                    &ci,
                    settings,
                    &knots,
                    &mut |time| Ok(project_f(force(time)?)),
                    limits,
                )?;
                damping_dissipation += d;
                external_work += w;
            } else if i > 0 {
                let next = project_f(force(time)?);
                let previous_u = state.displacement.clone();
                let previous_v = state.velocity.clone();
                integrator.as_ref().ok_or_else(|| Error::Internal("missing fixed-step integrator".into()))?.advance(&mut state, &f, &next)?;
                // Midpoint power / trapezoidal force work. Average-acceleration
                // Newmark satisfies delta-u = dt*v_mid, so its discrete balance
                // closes to roundoff for linear time-invariant operators.
                let midpoint_v: Vec<_> = previous_v.iter().zip(&state.velocity).map(|(a, b)| 0.5 * (a + b)).collect();
                damping_dissipation += 2.0 * t.time_step * energy(&ci, &midpoint_v);
                external_work +=
                    f.iter().zip(&next).zip(previous_u.iter().zip(&state.displacement)).map(|((a, b), (u, v))| 0.5 * (a + b) * (v - u)).sum::<f64>();
                f = next;
            }
            if i % t.output_every != 0 && i != t.increments {
                continue;
            }
            let expand_state = |q: &[f64]| match &modal {
                Some(m) => expand(q, m.shapes(), kr.n),
                None => q.to_vec(),
            };
            let u = expand_state(&state.displacement);
            let v = expand_state(&state.velocity);
            let a = expand_state(&state.acceleration);
            let full = red.expand(&u, false);
            let (displacement, rotation, _) = analysis::split(&full, &mesh);
            let velocity = analysis::split(&red.expand(&v, false), &mesh).0;
            let acceleration = analysis::split(&red.expand(&a, false), &mesh).0;
            let stress = if stress { Some(crate::harmonic::stresses(&mesh, &full, &displacement)?) } else { None };
            energies.push(Energy {
                damping_dissipation,
                external_work,
                energy_balance: energy(&mi, &state.velocity) + energy(&ki, &state.displacement) - initial_energy + damping_dissipation
                    - external_work,
            });
            output.frames.push(Frame {
                time: Some(time),
                displacement,
                rotation,
                velocity,
                acceleration,
                stress,
                kinetic_energy: Some(energy(&mr, &v)),
                strain_energy: Some(energy(&kr, &u)),
            });
        }
        output.energy = Some(energies);
        if t.adaptive.is_some() {
            output.warnings.push(format!("Adaptive local displacement control: accepted={}, rejected={}, min_step_s={}, max_step_s={}; estimates are not global error bounds; saved frames retain the requested grid.",
                controller.accepted, controller.rejected, controller.minimum_used, controller.maximum_used));
        }
    }
    output.validate().map_err(Error::Invalid)?;
    Ok(output)
}
fn combine_fields<const N: usize>(fields: &[Vec<[f64; N]>], nodes: usize, combine: &impl Fn(&[f64]) -> Result<f64>) -> Result<Vec<[f64; N]>> {
    let mut out = vec![[0.0; N]; nodes];
    for (i, v) in out.iter_mut().enumerate() {
        fieldcraft_model::progress::check()?;
        for (c, x) in v.iter_mut().enumerate() {
            let values: Vec<_> = fields.iter().filter_map(|f| f.get(i).map(|v| v[c])).collect();
            *x = combine(&values)?;
        }
    }
    Ok(out)
}
impl Response {
    pub fn validate(&self) -> std::result::Result<(), String> {
        if let Some(e) = &self.energy
            && (self.spectrum
                || e.len() != self.frames.len()
                || e.iter().any(|v| !v.external_work.is_finite() || !v.damping_dissipation.is_finite() || !v.energy_balance.is_finite()))
        {
            return Err("invalid dynamic energy ledger".into());
        }

        if self.nodes.is_empty()
            || self.nodes.first() == Some(&0)
            || self.nodes.windows(2).any(|w| w[0] >= w[1])
            || self.frames.is_empty()
            || self.natural_frequencies.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err("invalid dynamics nodes, frames or natural frequencies".into());
        }
        if self.spectrum && self.frames.len() != 1 {
            return Err("spectrum requires one peak frame".into());
        }
        if !self.spectrum && (self.frames.first().and_then(|f| f.time) != Some(0.0) || self.frames.windows(2).any(|w| w[0].time >= w[1].time)) {
            return Err("transient time axis must ascend from zero".into());
        }
        for f in &self.frames {
            if f.time.is_none() != self.spectrum || f.time.is_some_and(|t| !t.is_finite() || t < 0.0) {
                return Err("invalid dynamic frame time".into());
            }
            for v in [&f.displacement, &f.rotation, &f.acceleration] {
                if v.len() != self.nodes.len() || v.iter().flatten().any(|x| !x.is_finite()) {
                    return Err("dynamic vector length or finite-value mismatch".into());
                }
            }
            if (!self.spectrum && f.velocity.len() != self.nodes.len())
                || (self.spectrum && !f.velocity.is_empty())
                || f.velocity.iter().flatten().any(|x| !x.is_finite())
            {
                return Err("invalid dynamics velocity".into());
            }
            if f.stress.as_ref().is_some_and(|s| s.len() != self.nodes.len() || s.iter().flatten().any(|x| !x.is_finite())) {
                return Err("invalid dynamic stress".into());
            }
            for e in [f.kinetic_energy, f.strain_energy] {
                if e.is_none() != self.spectrum || e.is_some_and(|v| !v.is_finite()) {
                    return Err("invalid dynamic energy".into());
                }
            }
        }
        Ok(())
    }
}
