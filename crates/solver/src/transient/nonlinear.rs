//! Geometric nonlinear Newmark/HHT in fixed Cartesian reference coordinates.
//! Bathe, Finite Element Procedures, nonlinear direct integration; Bathe & Cimento
//! (1980), doi:10.1016/0045-7825(80)90051-1, full Newton iteration.
//! Hughes (1987/2000), ch. 9; Hilber, Hughes & Taylor (1977),
//! doi:10.1002/eqe.4290050306. See docs/manual/structural/nonlinear-transient.md
//! for verified sources and limits. No native result serialization is defined.
use super::{Integration, State, factor, vector};
use crate::analysis::{self, Setup};
use crate::sparse::{Csr, add_block};
use crate::system::{self, Mesh};
use crate::{Error, Limits, MassKind, Result};
use fieldcraft_model::dynamics::{Integration as Choice, Method};
use fieldcraft_model::{ElementKind, Load, Model, Procedure};
use serde::Serialize;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn mul(a: &Csr, v: &[f64]) -> Vec<f64> {
    let mut out = vec![0.; a.n];
    a.mul(v, &mut out);
    out
}
fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0., f64::max)
}

/// Read-only evaluation: a rejected iteration never commits constitutive history.
pub struct Evaluation {
    pub force: Vec<f64>,
    pub tangent: Csr,
    pub strain_energy: f64,
}
/// Atomic time advancement; every Newton correction rebuilds the effective tangent.
#[derive(Clone, Copy)]
pub struct Newton<'a> {
    pub mass: &'a Csr,
    pub damping: &'a Csr,
    pub dt: f64,
    pub method: Integration,
    pub max_iterations: usize,
    pub tolerance: f64,
    pub limits: &'a Limits,
}
impl Newton<'_> {
    pub fn advance(
        &self,
        state: &mut State,
        old_force: &[f64],
        next_force: &[f64],
        evaluate: impl Fn(&[f64]) -> Result<Evaluation>,
    ) -> Result<(usize, f64)> {
        let n = self.mass.n;
        crate::harmonic::Matrices { stiffness: self.mass, mass: self.mass, viscous: Some(self.damping), loss: None }.validate(self.limits)?;
        for v in [&state.displacement[..], &state.velocity, &state.acceleration, old_force, next_force] {
            vector(v, n)?;
        }
        if !self.dt.is_finite()
            || self.dt <= 0.
            || !(1..=100).contains(&self.max_iterations)
            || !self.tolerance.is_finite()
            || self.tolerance <= 0.
            || self.tolerance > 1e-2
        {
            return Err(Error::Invalid("nonlinear time step / Newton controls".into()));
        }
        let (beta, gamma, alpha) = self.method.parameters()?;
        let bdt = beta * self.dt * self.dt;
        if !bdt.is_finite() || bdt <= 0. {
            return Err(Error::Invalid("nonlinear time step under/overflow".into()));
        }
        let up: Vec<_> =
            (0..n).map(|i| state.displacement[i] + self.dt * state.velocity[i] + self.dt * self.dt * (0.5 - beta) * state.acceleration[i]).collect();
        let vp: Vec<_> = (0..n).map(|i| state.velocity[i] + self.dt * (1. - gamma) * state.acceleration[i]).collect();
        let old = evaluate(&state.displacement)?;
        vector(&old.force, n)?;
        let old_cv = mul(self.damping, &state.velocity);
        let mut u = up.clone();
        for iteration in 0..=self.max_iterations {
            fieldcraft_model::progress::check()?;
            vector(&u, n)?;
            let e = evaluate(&u)?;
            vector(&e.force, n)?;
            let a: Vec<_> = u.iter().zip(&up).map(|(u, p)| (u - p) / bdt).collect();
            let v: Vec<_> = vp.iter().zip(&a).map(|(v, a)| v + gamma * self.dt * a).collect();
            vector(&a, n)?;
            vector(&v, n)?;
            let ma = mul(self.mass, &a);
            let cv = mul(self.damping, &v);
            let r: Vec<_> = (0..n)
                .map(|i| ma[i] + (1. + alpha) * (e.force[i] + cv[i] - next_force[i]) - alpha * (old.force[i] + old_cv[i] - old_force[i]))
                .collect();
            vector(&r, n)?;
            let scale = norm(&ma).max(norm(&e.force)).max(norm(&cv)).max(norm(next_force)).max(norm(&old.force)).max(norm(old_force));
            let relative = if scale > 0. { norm(&r) / scale } else { norm(&r) };
            if relative <= self.tolerance {
                *state = State { displacement: u, velocity: v, acceleration: a };
                return Ok((iteration, relative));
            }
            if iteration == self.max_iterations {
                break;
            }
            // Acceleration-scaled effective tangent: M+(1+a)gamma dt C+(1+a)beta dt² Kt.
            let effective = self.mass.add_scaled((1. + alpha) * gamma * self.dt, self.damping)?.add_scaled((1. + alpha) * bdt, &e.tangent)?;
            crate::complex::validate(&effective, self.limits)?;
            let upper = crate::cholesky::UpperColumns::from_symmetric(&effective);
            let fact = if let Some(bytes) = self.limits.max_factor_memory_bytes {
                crate::ldlt::Symmetric::with_memory_limit(&upper, self.limits.max_factor_entries, crate::ldlt::Kind::Auto, bytes).map(|(f, _)| f)
            } else {
                crate::ldlt::Symmetric::new(&upper, self.limits.max_factor_entries, crate::ldlt::Kind::Auto)
            }
            .map_err(|e| match e {
                crate::cholesky::FactorError::Cancelled => fieldcraft_model::Error::Cancelled.into(),
                crate::cholesky::FactorError::TooLarge { entries, limit } => Error::factor_limit(entries, limit),
                other => Error::Invalid(format!("nonlinear effective tangent factor: {other:?}")),
            })?;
            let mut da: Vec<_> = r.iter().map(|x| -x).collect();
            fact.solve(&mut da)?;
            for (u, da) in u.iter_mut().zip(da) {
                *u += bdt * da;
            }
        }
        Err(Error::NotConverged(format!("nonlinear transient Newton exceeded {} corrections; committed state unchanged", self.max_iterations)))
    }
}

/// St Venant–Kirchhoff rod2, total Lagrangian: E_G=(d.d-L²)/(2L²).
/// f = EA E_G d/L; df/dd=EA/L [d⊗d/L²+E_G I].
fn rod(x: &[[f64; 3]], u: &[f64], ea: f64) -> Result<(Vec<f64>, Vec<f64>, f64)> {
    if x.len() != 2 || u.len() != 6 || !ea.is_finite() || ea <= 0. {
        return Err(Error::Invalid("finite rod2 data".into()));
    }
    vector(u, 6)?;
    let dx: [f64; 3] = std::array::from_fn(|i| x[1][i] - x[0][i]);
    let l2 = dot(&dx, &dx);
    let l = libm::sqrt(l2);
    let d: [f64; 3] = std::array::from_fn(|i| dx[i] + u[i + 3] - u[i]);
    let d2 = dot(&d, &d);
    if !l2.is_finite() || l2 <= 0. || !d2.is_finite() || d2 <= 1e-20 * l2 {
        return Err(Error::Invalid("collapsed finite rod2".into()));
    }
    let strain = (d2 - l2) / (2. * l2);
    let f: Vec<_> = (0..6).map(|i| if i < 3 { -1. } else { 1. } * ea / l * strain * d[i % 3]).collect();
    let k: Vec<_> = (0..36)
        .map(|i| {
            let r = i / 6;
            let c = i % 6;
            let sign = if (r < 3) == (c < 3) { 1. } else { -1. };
            sign * ea / l * (d[r % 3] * d[c % 3] / l2 + if r % 3 == c % 3 { strain } else { 0. })
        })
        .collect();
    let energy = 0.5 * ea * l * strain * strain;
    vector(&f, 6)?;
    vector(&k, 36)?;
    if !energy.is_finite() {
        return Err(Error::Invalid("rod energy overflow".into()));
    }
    Ok((f, k, energy))
}

/// Ephemeral SI diagnostics only: deliberately not a linear Response or a native Result.
#[derive(Serialize)]
pub struct Frame {
    pub time: f64,
    pub displacement: Vec<[f64; 3]>,
    pub velocity: Vec<[f64; 3]>,
    pub acceleration: Vec<[f64; 3]>,
    pub kinetic_energy: f64,
    pub strain_energy: f64,
    pub external_work: f64,
    pub physical_dissipation: f64,
    pub energy_balance: f64,
    pub newton_corrections: usize,
    pub relative_residual: f64,
}
#[derive(Serialize)]
pub struct History {
    pub step: String,
    pub mass: MassKind,
    pub nodes: Vec<u32>,
    pub frames: Vec<Frame>,
    pub warnings: Vec<String>,
}

/// Experimental elastic rod2 path, dead nodal forces, fixed supports, constant M/C.
/// Unsupported model data are refused rather than silently linearized.
pub fn solve(model: &Model, name: Option<&str>, mass: MassKind, max_iterations: usize, tolerance: f64, limits: &Limits) -> Result<History> {
    let step = model.select_step(name, |p| matches!(p, Procedure::Transient(_)), "transient")?;
    let Procedure::Transient(t) = &step.procedure else {
        return Err(Error::Internal("transient step".into()));
    };
    t.validate()?;
    if !matches!(t.method, Method::Direct)
        || t.adaptive.is_some()
        || t.base_acceleration.is_some()
        || t.base_displacement.is_some()
        || t.damping.beta != 0.
        || matches!(t.integration, Choice::LinearAcceleration)
    {
        return Err(Error::NotSupported("nonlinear transient: modal, adaptive, base motion, stiffness damping or linear_acceleration".into()));
    }
    if step.output.iter().any(|v| !matches!(v, fieldcraft_model::OutputVariable::Displacement)) {
        return Err(Error::NotSupported("nonlinear transient stress/other requested output; only nodal histories and energies supported".into()));
    }
    system::refuse_contact(model, "nonlinear transient")?;
    if !model.foundations().is_empty()
        || !model.initial_fields().is_empty()
        || !model.couplings().is_empty()
        || !model.ties().is_empty()
        || !model.pretensions().is_empty()
    {
        return Err(Error::NotSupported("nonlinear transient foundations, initial fields, ties, pretension or linear couplings".into()));
    }
    for (name, m) in model.materials() {
        if m.plastic.is_some() || m.hyperelastic.is_some() || m.damping.is_some_and(|d| d != Default::default()) {
            return Err(Error::NotSupported(format!("nonlinear transient material `{name}` plasticity/hyperelasticity/material damping")));
        }
    }
    for name in &step.loads {
        if !matches!(model.load(name), Some(Load::Force { .. })) {
            return Err(Error::NotSupported(format!("nonlinear transient load `{name}`: only dead nodal force supported")));
        }
    }
    let resolved = model.resolve()?;
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    if mesh.trusses.is_empty()
        || !mesh.elements.is_empty()
        || !mesh.structural.is_empty()
        || !mesh.points.is_empty()
        || mesh.trusses.iter().any(|e| e.kind != ElementKind::Rod2)
        || !mesh.discrete.springs.is_empty()
    {
        return Err(Error::NotSupported(
            "nonlinear transient requires rod2, optional point masses, no solids/beams/shells/reference points/springs/dashpots".into(),
        ));
    }
    let mut warnings: Vec<_> = resolved.warnings.iter().map(ToString::to_string).chain(mesh.warnings.iter().cloned()).collect();
    let (graph, order) = analysis::pattern(&mesh);
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let red = analysis::reduction(&setup, true, &order, analysis::budget(limits, &mesh, &graph, 6), &mut warnings)?;
    if red.largest_prescribed() != 0. || red.neq() == 0 {
        return Err(Error::NotSupported("nonlinear transient requires free equations and homogeneous fixed supports".into()));
    }
    let frames = t.increments / t.output_every + 2;
    if mesh.labels.len().checked_mul(frames as usize).and_then(|n| n.checked_mul(9)).is_none_or(|n| n > limits.max_factor_entries.min(5_000_000)) {
        return Err(Error::Limit("nonlinear diagnostic frame budget".into()));
    }
    let assembled = analysis::assemble(&mesh, &graph, Some(mass))?;
    let mut m = assembled.m.ok_or_else(|| Error::Internal("nonlinear mass".into()))?;
    m.mirror_lower();
    let [mut mr] = red.reduce([&m])?;
    mr.mirror_lower();
    let mut cr = mr.clone();
    for c in &mut cr.val {
        *c *= t.damping.alpha;
    }
    let evaluate = |q: &[f64]| -> Result<Evaluation> {
        let full = red.expand(q, false);
        let mut force = vec![0.; mesh.dofs()];
        let mut k = graph.zero_matrix(&mesh.layout);
        let mut energy = 0.;
        for e in &mesh.trusses {
            let dofs = mesh.layout.element_dofs(&e.nodes, fieldcraft_fem::DofSet::TRANSLATIONS).ok_or_else(|| Error::Internal("rod dofs".into()))?;
            let ue: Vec<_> = dofs.iter().map(|d| full.get(*d).copied().unwrap_or(0.)).collect();
            let (f, ke, se) = rod(&e.xyz, &ue, e.young * e.area)?;
            energy += se;
            for (d, f) in dofs.iter().zip(f) {
                if let Some(v) = force.get_mut(*d) {
                    *v += f;
                }
            }
            add_block(&mut k, &mesh.layout, &e.nodes, fieldcraft_fem::DofSet::TRANSLATIONS, &ke)?;
        }
        if !energy.is_finite() {
            return Err(Error::Invalid("nonlinear strain energy overflow".into()));
        }
        let [mut tangent] = red.reduce([&k])?;
        tangent.mirror_lower();
        Ok(Evaluation { force: red.restrict(&force), tangent, strain_energy: energy })
    };
    let mut loads = Vec::new();
    for name in &step.loads {
        let mut one = step.clone();
        one.loads.clear();
        one.loads.insert(name.clone());
        loads.push((
            red.restrict(&system::loads(model, &resolved, &mesh, &one, &red.active)?),
            t.load_amplitudes.get(name).map(|a| a.evaluator()).transpose()?,
        ));
    }
    let force = |time: f64| -> Result<Vec<f64>> {
        let mut out = vec![0.; mr.n];
        for (f, curve) in &loads {
            let v = match curve {
                Some(a) => a.value(time)?,
                None => 1.,
            };
            for (o, f) in out.iter_mut().zip(f) {
                *o += v * f;
            }
        }
        vector(&out, mr.n)?;
        Ok(out)
    };
    let u = super::model::initial(&t.initial_displacement, &mesh, &red)?;
    let v = super::model::initial(&t.initial_velocity, &mesh, &red)?;
    let mut f = force(0.)?;
    let init = evaluate(&u)?;
    let cv = mul(&cr, &v);
    let mut a: Vec<_> = f.iter().zip(&init.force).zip(cv).map(|((f, i), c)| f - i - c).collect();
    factor(&mr, limits)?.solve(&mut a)?;
    vector(&a, mr.n)?;
    let initial_residual =
        norm(&mul(&mr, &a).iter().zip(&init.force).zip(mul(&cr, &v)).zip(&f).map(|(((m, i), c), f)| m + i + c - f).collect::<Vec<_>>());
    let scale = norm(&f).max(norm(&init.force)).max(norm(&mul(&cr, &v))).max(norm(&mul(&mr, &a)));
    let initial_residual = if scale > 0. { initial_residual / scale } else { initial_residual };
    let mut state = State { displacement: u, velocity: v, acceleration: a };
    let energy0 = 0.5 * dot(&state.velocity, &mul(&mr, &state.velocity)) + init.strain_energy;
    let (mut work, mut diss) = (0., 0.);
    let method = match t.integration {
        Choice::Hht { alpha } => Integration::Hht { alpha },
        _ => Integration::AverageAcceleration,
    };
    let newton = Newton { mass: &mr, damping: &cr, dt: t.time_step, method, max_iterations, tolerance, limits };
    let field = |q: &[f64]| -> Result<Vec<[f64; 3]>> {
        let full = red.expand(q, false);
        mesh.labels
            .iter()
            .enumerate()
            .map(|(n, _)| {
                let ds = system::translations(&mesh, n)?;
                Ok(std::array::from_fn(|i| ds.get(i).and_then(|d| full.get(*d)).copied().unwrap_or(0.)))
            })
            .collect()
    };
    let mut total_iterations = 0usize;
    let mut out = History { step: step.name.clone(), mass, nodes: mesh.labels.clone(), frames: Vec::new(), warnings };
    out.warnings.push("Experimental geometric rod2 dynamics: SI diagnostics only, not persisted. Energy balance includes time-integration and work-quadrature error; it is not physical damping.".into());
    for i in 0..=t.increments {
        let (mut iterations, mut residual) = (0, initial_residual);
        if i > 0 {
            if total_iterations >= limits.max_newton_iterations {
                return Err(Error::Limit("nonlinear total Newton budget".into()));
            }
            let nf = force(i as f64 * t.time_step)?;
            let old = state.clone();
            (iterations, residual) = Newton { max_iterations: newton.max_iterations.min(limits.max_newton_iterations - total_iterations), ..newton }
                .advance(&mut state, &f, &nf, evaluate)
                .map_err(|e| match e {
                    Error::NotConverged(m) => {
                        Error::NotConverged(format!("time step {i}, last committed time {}: {m}", (i - 1) as f64 * t.time_step))
                    }
                    other => other,
                })?;
            total_iterations += iterations;
            work +=
                state.displacement.iter().zip(&old.displacement).zip(f.iter().zip(&nf)).map(|((u, o), (f, n))| 0.5 * (u - o) * (f + n)).sum::<f64>();
            diss += 0.5 * t.time_step * (dot(&old.velocity, &mul(&cr, &old.velocity)) + dot(&state.velocity, &mul(&cr, &state.velocity)));
            f = nf;
            if !work.is_finite() || !diss.is_finite() {
                return Err(Error::Invalid("nonlinear work overflow".into()));
            }
        }
        if i % t.output_every == 0 || i == t.increments {
            let kinetic = 0.5 * dot(&state.velocity, &mul(&mr, &state.velocity));
            let strain = evaluate(&state.displacement)?.strain_energy;
            out.frames.push(Frame {
                time: i as f64 * t.time_step,
                displacement: field(&state.displacement)?,
                velocity: field(&state.velocity)?,
                acceleration: field(&state.acceleration)?,
                kinetic_energy: kinetic,
                strain_energy: strain,
                external_work: work,
                physical_dissipation: diss,
                energy_balance: kinetic + strain - energy0 + diss - work,
                newton_corrections: iterations,
                relative_residual: residual,
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn finite_rod_force_tangent_and_energy_are_consistent_under_rotation() {
        let x = [[0., 0., 0.], [1., 0., 0.]];
        let u = [0., 0., 0., -0.8, 1.1, 0.3];
        let (f, k, _) = rod(&x, &u, 17.).unwrap();
        let eps = 1e-6;
        for j in 0..6 {
            let mut p = u;
            let mut m = u;
            p[j] += eps;
            m[j] -= eps;
            let (fp, _, ep) = rod(&x, &p, 17.).unwrap();
            let (fm, _, em) = rod(&x, &m, 17.).unwrap();
            assert!(((ep - em) / (2. * eps) - f[j]).abs() < 1e-8);
            for i in 0..6 {
                assert!(((fp[i] - fm[i]) / (2. * eps) - k[i * 6 + j]).abs() < 1e-8);
            }
        }
        let (f, _, energy) = rod(&x, &[0., 0., 0., -1., 1., 0.], 17.).unwrap();
        assert_eq!(energy, 0.);
        assert!(f.iter().all(|v| *v == 0.));
    }
    #[test]
    fn failure_restores_all_kinematic_state_and_hht_zero_matches_newmark() {
        let matrix = |v| Csr { n: 1, ptr: vec![0, 1], idx: vec![0], val: vec![v] };
        let m = matrix(1.);
        let c = matrix(0.);
        let limits = Limits::default();
        let eval = |u: &[f64]| {
            Ok(Evaluation {
                force: vec![u[0] + u[0].powi(3)],
                tangent: matrix(1. + 3. * u[0] * u[0]),
                strain_energy: 0.5 * u[0] * u[0] + 0.25 * u[0].powi(4),
            })
        };
        let initial = State { displacement: vec![1.], velocity: vec![0.], acceleration: vec![-2.] };
        let mut state = initial.clone();
        let mut solver =
            Newton { mass: &m, damping: &c, dt: 1., method: Integration::AverageAcceleration, max_iterations: 1, tolerance: 1e-10, limits: &limits };
        assert!(solver.advance(&mut state, &[0.], &[4.], eval).is_err());
        assert_eq!(state.displacement, initial.displacement);
        assert_eq!(state.velocity, initial.velocity);
        assert_eq!(state.acceleration, initial.acceleration);
        solver.max_iterations = 20;
        solver.dt = 0.1;
        solver.advance(&mut state, &[0.], &[0.], eval).unwrap();
        let mut other = initial;
        solver.method = Integration::Hht { alpha: 0. };
        solver.advance(&mut other, &[0.], &[0.], eval).unwrap();
        assert_eq!(state.displacement, other.displacement);
        assert_eq!(state.velocity, other.velocity);
    }
}
