//! The time loop: central differences with variable increments, kinematic constraints, the
//! energy account and output.
//!
//! One increment (Belytschko, Liu, Moran, Elkhodary 2014, Box 6.1, with Δt^n = ½(Δt^{n−½} +
//! Δt^{n+½})):
//! 1. u^n = u^{n−1} + Δt^{n−½} v^{n−½}; t^n = t^{n−1} + Δt^{n−½};
//! 2. internal forces f^n from u^n (stresses, hourglass stiffness) and v^{n−½} (hourglass and
//!    bulk viscosity), and the next increment Δt^{n+½} = scale · min_e Δt_e;
//! 3. the contact forces c^n at u^n ([`crate::contact`], with the stable increment Δt^{n+½});
//!    v^{n+½} = v^{n−½} + Δt^n M⁻¹(f_ext + c^n − f^n), then the constraints: rigid walls set the
//!    normal velocity of a node that would cross its wall so it lands on it (a kinematic contact
//!    condition: T. J. R. Hughes, R. L. Taylor, J. L. Sackman, A. Curnier, W. Kanoknukulchai,
//!    *A finite element method for a class of contact-impact problems*, CMAME 8(3):249–276,
//!    1976), with Coulomb friction as a tangential impulse at most μ times the normal one (the
//!    tangential velocity stops, or loses μ Δv_n when it slides), prescribed degrees of
//!    freedom take their velocity, and a symmetry plane not normal to an axis removes its nodes'
//!    velocity along its normal (from their free components: the lumped mass is the same on a
//!    node's three translations, so this is the mass-weighted projection onto the constraint). The change a constraint makes is the impulse of its reaction,
//!    r^n = M Δv / Δt^n;
//! 4. energies: the works of the internal, contact and external forces by the trapezoidal rule,
//!    W += ½ Δu · (f^{n−1} + f^n) (exactly the strain energy for a linear elastic material),
//!    the constraints' work Δt^n r^n · ½(v^{n−½} + v^{n+½}) (exactly the kinetic energy their
//!    impulses add or remove), and the kinetic energy of v^n = v^{n−½} + ½Δt^{n−½} a^n.
//!
//! The run starts from u = 0 and the initial velocities with the half increment
//! v^{½} = v^0 + ½Δt^{½} a^0.

/// Increments between progress reports (each takes a lock).
const PROGRESS_EVERY: u64 = 64;

use fieldcraft_model::{ExplicitDynamic, Model, OutputVariable, Procedure, Step};

use crate::contact::{Evaluation, Interface};
use crate::elements::{Motion, State, Workspace, internal_forces};
use crate::result::{ContactSummary, ControllingElement, Energies, ErodedElement, ExplicitResult, Frame, History};
use crate::system::{NODE_DOFS, System};
use crate::{Error, Limits, Result};

/// How many elements the result names as controlling the time increment.
pub const CONTROLLING: usize = 10;

/// Relative energy balance error above which a run warns.
pub const BALANCE_TOLERANCE: f64 = 0.01;

/// Hourglass energy as a fraction of the largest internal energy of the run above which it warns
/// (the internal energy may fall back to zero, as in a rebound).
pub const HOURGLASS_TOLERANCE: f64 = 0.1;

/// The step named `name`, or the model's only explicit dynamic step.
fn select_step<'m>(model: &'m Model, name: Option<&str>) -> Result<(&'m Step, &'m ExplicitDynamic)> {
    let step = model.select_step(name, Procedure::is_dynamic, "explicit dynamic")?;
    match &step.procedure {
        Procedure::ExplicitDynamic(e) => Ok((step, e)),
        _ => Err(Error::Internal("an explicit dynamic step".into())),
    }
}

/// Times at a fixed interval: whether `t` has reached the next one.
struct Schedule {
    interval: f64,
    next: f64,
}

impl Schedule {
    fn due(&mut self, t: f64) -> bool {
        if t < self.next * (1.0 - 1e-12) {
            return false;
        }
        self.next = ((t / self.interval).floor() + 1.0) * self.interval;
        true
    }
}

/// What the run writes.
struct Output {
    displacement: bool,
    stress: bool,
    strain: bool,
}

/// Velocity constraints over an increment of length `dt` (rigid walls when dt > 0, then
/// prescribed degrees of freedom, then inclined symmetry planes); returns the kinetic energy the walls' friction took (`mass`
/// per degree of freedom).
fn constrain(sys: &System, u: &[f64], v: &mut [f64], dt: f64, mass: &[f64]) -> f64 {
    let mut dissipated = 0.0;
    if dt > 0.0 {
        for w in &sys.walls {
            for &node in &w.nodes {
                let x = sys.xyz.get(node).copied().unwrap_or_default();
                let at = |k: usize| NODE_DOFS * node + k;
                let gap: f64 = (0..3).map(|k| (x[k] + u.get(at(k)).copied().unwrap_or(0.0) - w.point[k]) * w.normal[k]).sum();
                let normal: f64 = (0..3).map(|k| v.get(at(k)).copied().unwrap_or(0.0) * w.normal[k]).sum();
                if gap + dt * normal < 0.0 {
                    let change = -gap / dt - normal;
                    for k in 0..3 {
                        if let Some(c) = v.get_mut(at(k)) {
                            *c += change * w.normal[k];
                        }
                    }
                    if !w.friction.is_frictionless() {
                        dissipated += wall_friction(v, node, w, change, mass.get(at(0)).copied().unwrap_or(0.0));
                    }
                }
            }
        }
    }
    for (c, p) in v.iter_mut().zip(&sys.prescribed) {
        if let Some(p) = p {
            *c = *p;
        }
    }
    for row in &sys.planes {
        row.project(v);
    }
    dissipated
}

/// Coulomb friction of a node the wall `w` stopped with a normal velocity change `change`: its
/// tangential velocity changes by at most μ · change (stick when that stops it); returns the
/// kinetic energy taken, ½ m (|v_t|² − |v_t′|²).
fn wall_friction(v: &mut [f64], node: usize, w: &crate::system::Wall, change: f64, mass: f64) -> f64 {
    let at = |k: usize| NODE_DOFS * node + k;
    let vel = [0, 1, 2].map(|k| v.get(at(k)).copied().unwrap_or(0.0));
    let vn: f64 = (0..3).map(|k| vel[k] * w.normal[k]).sum();
    let tangential = [0, 1, 2].map(|k| vel[k] - vn * w.normal[k]);
    let slide = tangential.iter().map(|c| c * c).sum::<f64>().sqrt();
    if slide == 0.0 {
        return 0.0;
    }
    let cap = w.friction.coefficient(slide) * change;
    let keep = if slide <= cap { 0.0 } else { 1.0 - cap / slide };
    for k in 0..3 {
        if let Some(c) = v.get_mut(at(k)) {
            *c = vn * w.normal[k] + keep * tangential[k];
        }
    }
    0.5 * mass * slide * slide * (1.0 - keep * keep)
}

/// v_out = v_in + dt M⁻¹ (f_ext + f_contact − f_int).
fn accelerate(sys: &System, ws: &Workspace, mass: &[f64], v_in: &[f64], dt: f64, v_out: &mut [f64]) {
    let internal = ws.stress_force.iter().zip(&ws.hourglass_force).zip(&ws.viscous_force);
    for ((((out, v), m), e), ((s, h), q)) in v_out.iter_mut().zip(v_in).zip(mass).zip(&sys.external).zip(internal) {
        *out = v + dt * ((e - (s + h + q)) / m);
    }
    // Contact adds its forces (a separate pass, so a model without contact runs the same sums).
    if !ws.contact_force.is_empty() {
        for ((out, m), c) in v_out.iter_mut().zip(mass).zip(&ws.contact_force) {
            *out += dt * (c / m);
        }
    }
}

/// The contacts' forces at `u` (with the increment `du` over `dt_prev` that led there, and the
/// stable increment `dt`) into `ws.contact_force`, and what each found.
#[allow(clippy::too_many_arguments)]
fn contact_forces(
    contacts: &mut [Interface],
    sys: &System,
    u: &[f64],
    du: &[f64],
    dt_prev: f64,
    dt: f64,
    mass: &[f64],
    state: &State,
    ws: &mut Workspace,
) -> Vec<Evaluation> {
    ws.contact_force.fill(0.0);
    contacts.iter_mut().map(|c| c.forces(&sys.xyz, u, du, dt_prev, dt, mass, state, &mut ws.contact_force)).collect()
}

/// A contact's record over the run.
fn track(summary: &mut ContactSummary, e: &Evaluation, t: f64) {
    if e.active > 0 {
        summary.first_contact.get_or_insert(t);
        summary.last_contact = Some(t);
    }
    if e.normal_force > summary.peak_force {
        summary.peak_force = e.normal_force;
        summary.peak_time = t;
    }
    summary.max_penetration = summary.max_penetration.max(e.penetration);
    summary.friction += e.dissipation;
}

/// The constraint reactions M (v_constrained − v_free) / Δt, into `r`.
fn impulses(mass: &[f64], free: &[f64], constrained: &[f64], dt: f64, r: &mut [f64]) {
    for (r, (m, (f, c))) in r.iter_mut().zip(mass.iter().zip(free.iter().zip(constrained))) {
        *r = if c == f { 0.0 } else { m * (c - f) / dt };
    }
}

fn per_node(x: &[f64]) -> Vec<[f64; 3]> {
    x.as_chunks::<NODE_DOFS>().0.to_vec()
}

/// Solve an explicit dynamic step: the step named `step`, or the model's only one.
pub fn solve_explicit(model: &Model, step: Option<&str>, limits: &Limits) -> Result<ExplicitResult> {
    let resolved = model.resolve()?;
    let (step, settings) = select_step(model, step)?;
    let mut sys = System::new(model, &resolved, step, settings, limits)?;
    let end = settings.end_time;
    let scale = settings.dt_scale;
    let bulk = settings.bulk_viscosity;
    let max_increments = settings.max_increments.min(limits.max_increments);
    let mut warnings = sys.warnings.clone();

    // The controlling elements: the smallest critical increments at rest.
    let mut critical: Vec<ControllingElement> = sys
        .hexes
        .iter()
        .map(|h| (h.label, h.length, h.material, h.mass_factor))
        .chain(sys.tets.iter().map(|t| (t.label, t.length, t.material, t.mass_factor)))
        .chain(sys.rods.iter().map(|r| (r.label, r.length, r.material, 1.0)))
        .map(|(label, length, m, f)| ControllingElement {
            label,
            // Mass scaling slows the element's waves by √f.
            increment: length * f.sqrt() / sys.materials.get(m).map_or(f64::INFINITY, |m| m.wave_speed),
        })
        .collect();
    critical.sort_by(|a, b| a.increment.total_cmp(&b.increment).then(a.label.cmp(&b.label)));
    critical.truncate(CONTROLLING);
    let n = sys.dofs();
    let mass = sys.dof_mass();
    if let Some(d) = mass.iter().position(|m| !(*m > 0.0 && m.is_finite())) {
        let label = sys.labels.get(crate::discrete::node_direction(d).0).copied().unwrap_or(0);
        return Err(Error::Invalid(format!(
            "node {label} has no mass in direction {}; explicit dynamics needs a mass at every node (give the nodes only springs or dashpots use a point mass)",
            crate::discrete::node_direction(d).1 + 1
        )));
    }
    // Springs and dashpots have their own stable increment on the lumped mass.
    let elements_smallest = critical.first().map(|c| c.increment);
    if let Some((l, name)) = &sys.discrete_increment
        && elements_smallest.is_none_or(|e| *l < e)
    {
        warnings.push(format!("spring or dashpot `{name}` limits the time increment to {l:.3e} (before the scale factor)"));
    }
    let smallest = elements_smallest
        .into_iter()
        .chain(sys.discrete_increment.as_ref().map(|l| l.0))
        .reduce(f64::min)
        .filter(|d| *d > 0.0 && d.is_finite())
        .ok_or_else(|| Error::Invalid("no element has a finite, positive critical time increment (check the material constants)".into()))?;
    let estimate = end / (scale * smallest);
    if estimate > max_increments as f64 {
        let first = critical.first().map_or(0, |c| c.label);
        return Err(Error::Limit(format!(
            "the step needs about {estimate:.3e} increments of {:.3e} (element {first} controls it), more than the cap of {max_increments}",
            scale * smallest
        )));
    }

    // Output sizes.
    let output = Output {
        displacement: step.output.contains(&OutputVariable::Displacement),
        stress: step.output.contains(&OutputVariable::Stress),
        strain: step.output.contains(&OutputVariable::Strain),
    };
    let (nodes, elements) = (sys.labels.len(), sys.elements.len());
    let frames = (end / settings.output_interval()).ceil() + 2.0;
    let samples = (end / settings.history_interval()).ceil() + 2.0;
    let per_frame = (if output.displacement { 6 * nodes } else { 0 }) + 6 * elements * (usize::from(output.stress) + usize::from(output.strain));
    let per_sample = 16 + 6 * sys.history_nodes.len();
    let values = frames * per_frame as f64 + samples * per_sample as f64;
    if values > limits.max_output_values as f64 {
        return Err(Error::Limit(format!(
            "the output would hold {values:.3e} numbers; the limit is {} (lengthen the output or history interval, or request less output)",
            limits.max_output_values
        )));
    }

    // Point masses count with their mean over the three directions.
    let total_mass: f64 = sys.mass.iter().sum::<f64>() + sys.discrete.mass.iter().sum::<f64>() / 3.0;
    let mut ws = Workspace::new(&sys);
    let mut previous = Workspace::new(&sys);
    let mut contacts = std::mem::take(&mut sys.contacts);
    let mut summaries: Vec<ContactSummary> =
        contacts.iter().map(|c| ContactSummary { name: c.name.clone(), adjusted: c.adjusted.clone(), ..ContactSummary::default() }).collect();
    let mut state = State::new(&sys);
    let mut eroded: Vec<ErodedElement> = Vec::new();
    let mut u = vec![0.0; n];
    let mut du = vec![0.0; n];
    let v0 = sys.initial_velocity.clone();
    let initial_kinetic = 0.5 * v0.iter().zip(&mass).map(|(v, m)| m * v * v).sum::<f64>();
    let inelastic = sys.materials.iter().any(|m| m.j2.as_ref().is_some_and(|j| j.flow.is_some()));

    // Start: forces at rest with the initial velocities, then the half increment.
    let start = Motion { u: &u, du: &du, v: &v0, dt: 0.0, time: 0.0, bulk: &bulk };
    let crit = internal_forces(&sys, &start, &mut state, &mut ws)?;
    let mut dt_half = (scale * crit).min(end);
    let initial_increment = dt_half;
    let mut found = contact_forces(&mut contacts, &sys, &u, &du, 0.0, scale * crit, &mass, &state, &mut ws);
    for (summary, e) in summaries.iter_mut().zip(&found) {
        track(summary, e, 0.0);
    }
    let mut v_half = vec![0.0; n];
    accelerate(&sys, &ws, &mass, &v0, 0.5 * dt_half, &mut v_half);
    let mut free = v_half.clone();
    let wall_friction = constrain(&sys, &u, &mut v_half, dt_half, &mass);
    let mut reaction = vec![0.0; n];
    impulses(&mass, &free, &v_half, 0.5 * dt_half, &mut reaction);
    // The constraints' work is that of their impulses on the mean of the velocities before
    // and after the increment, Δt r · ½(v^{n−½} + v^{n+½}): exact for the kinetic energy a
    // kinematic constraint removes or adds (a node stopped by a wall loses ½ m v²).
    let impulse_work = |r: &[f64], before: &[f64], after: &[f64], dt: f64| -> f64 {
        dt * r.iter().zip(before.iter().zip(after)).map(|(r, (a, b))| r * 0.5 * (a + b)).sum::<f64>()
    };
    let mut energy = Energies {
        kinetic: initial_kinetic,
        constraint: impulse_work(&reaction, &v0, &v_half, 0.5 * dt_half),
        friction: wall_friction,
        ..Energies::default()
    };

    let mut history = History { nodes: sys.history_nodes.iter().filter_map(|i| sys.labels.get(*i).copied()).collect(), ..History::default() };
    let mut out_frames = Vec::new();
    let sample = |history: &mut History, t: f64, dt: f64, e: Energies, u: &[f64], v: &[f64], r: &[f64], found: &[Evaluation]| {
        history.time.push(t);
        history.contact_force.push(found.iter().map(|e| e.normal_force).collect());
        history.increment.push(dt);
        history.energies.push(e);
        let sum = |x: &[f64], weight: &dyn Fn(usize, usize) -> f64| {
            let mut s = [0.0; 3];
            for (node, c) in x.as_chunks::<NODE_DOFS>().0.iter().enumerate() {
                for k in 0..3 {
                    s[k] += weight(node, k) * c[k];
                }
            }
            s
        };
        history.momentum.push(sum(v, &|node, k| mass.get(crate::discrete::dof(node, k)).copied().unwrap_or(0.0)));
        history.constraint_force.push(sum(r, &|_, _| 1.0));
        let pick = |x: &[f64]| sys.history_nodes.iter().map(|i| [0, 1, 2].map(|k| x.get(NODE_DOFS * i + k).copied().unwrap_or(0.0))).collect();
        history.displacement.push(pick(u));
        history.velocity.push(pick(v));
    };
    let frame = |t: f64, u: &[f64], v: &[f64], state: &State| {
        let each = |on: bool, f: &dyn Fn(&crate::elements::ElementState) -> f64| -> Vec<f64> {
            if on { state.elements.iter().map(f).collect() } else { Vec::new() }
        };
        let inelastic = output.stress && inelastic;
        Frame {
            time: t,
            displacement: if output.displacement { per_node(u) } else { Vec::new() },
            velocity: if output.displacement { per_node(v) } else { Vec::new() },
            stress: if output.stress { state.elements.iter().map(|e| e.point.stress).collect() } else { Vec::new() },
            strain: if output.strain { state.elements.iter().map(|e| e.strain).collect() } else { Vec::new() },
            plastic_strain: each(inelastic, &|e| e.point.plastic_strain),
            damage: each(inelastic, &|e| e.point.damage),
            // An element whose material has no temperature reports none (not 0 K).
            temperature: each(inelastic, &|e| if e.has_temperature { e.point.temperature } else { f64::NAN }),
            failed: if output.stress { state.elements.iter().map(|e| e.failed.is_some()).collect() } else { Vec::new() },
        }
    };
    let plastic = |state: &State| state.elements.iter().map(|e| e.plastic_work).sum::<f64>();
    let mut histories = Schedule { interval: settings.history_interval(), next: settings.history_interval() };
    let mut frames_due = Schedule { interval: settings.output_interval(), next: settings.output_interval() };
    sample(&mut history, 0.0, dt_half, energy, &u, &v0, &reaction, &found);
    out_frames.push(frame(0.0, &u, &v0, &state));

    let mut t = 0.0;
    let mut increments = 0u64;
    let mut peak_internal: f64 = 0.0;
    let mut v_now = v0;
    let mut completed = true;
    let mut v_new = vec![0.0; n];
    fieldcraft_model::progress::stage("Time", Some(end), "");
    loop {
        if increments.is_multiple_of(PROGRESS_EVERY) {
            fieldcraft_model::progress::count(t);
            fieldcraft_model::progress::check()?;
        }
        if increments >= max_increments {
            completed = false;
            warnings.push(format!("stopped at time {t:.6e} after {increments} increments (the cap); the step ends at {end:.6e}"));
            break;
        }
        for ((x, dx), v) in u.iter_mut().zip(du.iter_mut()).zip(&v_half) {
            *dx = dt_half * v;
            *x += *dx;
        }
        t += dt_half;
        increments += 1;
        std::mem::swap(&mut ws, &mut previous);
        let motion = Motion { u: &u, du: &du, v: &v_half, dt: dt_half, time: t, bulk: &bulk };
        let crit = internal_forces(&sys, &motion, &mut state, &mut ws)?;
        for &(index, failure) in &state.eroded {
            let element = state.elements.get(index);
            let internal = element.map_or(0.0, |e| e.energy);
            energy.eroded += internal;
            eroded.push(ErodedElement { label: sys.elements.get(index).copied().unwrap_or(0), time: t, failure, internal_energy: internal });
        }
        let remaining = end - t;
        let last = remaining <= 1e-12 * end;
        let dt_next = if last { 0.0 } else { (scale * crit).min(remaining) };
        let dt_n = 0.5 * (dt_half + dt_next);
        if !contacts.is_empty() {
            found = contact_forces(&mut contacts, &sys, &u, &du, dt_half, scale * crit, &mass, &state, &mut ws);
            for (summary, e) in summaries.iter_mut().zip(&found) {
                track(summary, e, t);
                energy.friction += e.dissipation;
            }
        }
        accelerate(&sys, &ws, &mass, &v_half, dt_n, &mut v_new);
        free.clone_from(&v_new);
        energy.friction += constrain(&sys, &u, &mut v_new, dt_next, &mass);
        impulses(&mass, &free, &v_new, dt_n, &mut reaction);
        // One pass over the degrees of freedom: v^n = v^{n−½} + ½Δt^{n−½} a^n with the constrained
        // a^n = (v^{n+½} − v^{n−½}) / Δt^n, and the works (trapezoidal over the displacement for
        // the internal and external forces).
        let weight = 0.5 * dt_half / dt_n;
        let [mut internal, mut hourglass, mut viscous, mut external, mut constraint, mut kinetic] = [0.0; 6];
        let forces = ws.stress_force.iter().zip(&ws.hourglass_force).zip(&ws.viscous_force);
        let before = previous.stress_force.iter().zip(&previous.hourglass_force).zip(&previous.viscous_force);
        let per_dof = v_now.iter_mut().zip(&v_half).zip(&v_new).zip(&mass).zip(&du).zip(&sys.external).zip(&reaction);
        for ((((((((now, vh), vn), m), dx), e), r), ((s, h), q)), ((s0, h0), q0)) in per_dof.zip(forces).zip(before) {
            *now = vh + weight * (vn - vh);
            internal += dx * (s0 + s);
            hourglass += dx * (h0 + h);
            viscous += dx * (q0 + q);
            external += dx * e;
            constraint += r * 0.5 * (vh + vn);
            kinetic += m * *now * *now;
        }
        energy.internal += 0.5 * internal;
        energy.hourglass += 0.5 * hourglass;
        energy.viscous += 0.5 * viscous;
        energy.external += external;
        energy.constraint += dt_n * constraint;
        energy.kinetic = 0.5 * kinetic;
        if !sys.added_mass.is_empty() {
            let per_node = v_now.as_chunks::<NODE_DOFS>().0.iter().zip(&sys.added_mass);
            energy.added_kinetic = 0.5 * per_node.map(|(v, m)| m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])).sum::<f64>();
        }
        if !ws.contact_force.is_empty() {
            let work: f64 = du.iter().zip(ws.contact_force.iter().zip(&previous.contact_force)).map(|(d, (c, c0))| d * (c0 + c)).sum();
            energy.contact -= 0.5 * work;
        }
        peak_internal = peak_internal.max(energy.internal.abs());
        energy.balance = energy.kinetic + energy.internal + energy.hourglass + energy.viscous + energy.contact
            - energy.external
            - energy.constraint
            - initial_kinetic;
        // Every element eroded leaves an infinite critical increment: the nodes fly freely.
        if !(energy.balance.is_finite() && !crit.is_nan()) {
            return Err(Error::Unstable(format!("the solution is not finite at time {t:.6e} (increment {increments})")));
        }
        if histories.due(t) || last {
            energy.plastic = plastic(&state);
            sample(&mut history, t, dt_half, energy, &u, &v_now, &reaction, &found);
        }
        if frames_due.due(t) || last {
            out_frames.push(frame(t, &u, &v_now, &state));
        }
        std::mem::swap(&mut v_half, &mut v_new);
        dt_half = dt_next;
        if last {
            break;
        }
    }
    if !completed {
        if history.time.last() != Some(&t) {
            sample(&mut history, t, dt_half, energy, &u, &v_now, &reaction, &found);
        }
        if out_frames.last().map(|f| f.time) != Some(t) {
            out_frames.push(frame(t, &u, &v_now, &state));
        }
    }
    energy.plastic = plastic(&state);
    if !eroded.is_empty() {
        warnings.push(format!(
            "{} elements eroded (their material failed); their internal energy at erosion was {:.4e}",
            eroded.len(),
            energy.eroded
        ));
    }

    let error = energy.relative_error(initial_kinetic);
    if error > BALANCE_TOLERANCE {
        warnings.push(format!("the energy balance is off by {:.2}% of the largest energy (kinetic, internal or work)", 100.0 * error));
    }
    // Round-off hourglass energy in a model at rest or in rigid motion is not worth a warning.
    let floor =
        1e-9 * [initial_kinetic, energy.kinetic, peak_internal, energy.external.abs(), energy.constraint.abs()].into_iter().fold(0.0, f64::max);
    if energy.hourglass.abs() > HOURGLASS_TOLERANCE * peak_internal && energy.hourglass.abs() > floor {
        warnings.push(format!(
            "hourglass energy is {:.1}% of the largest internal energy; refine the mesh or change the hourglass control",
            100.0 * energy.hourglass.abs() / peak_internal.max(f64::MIN_POSITIVE)
        ));
    }
    Ok(ExplicitResult {
        step: step.name.clone(),
        nodes: sys.labels.clone(),
        elements: sys.elements.clone(),
        frames: out_frames,
        history,
        controlling: critical,
        initial_increment,
        increments,
        time: t,
        completed,
        total_mass,
        energy,
        initial_kinetic,
        eroded,
        contacts: summaries,
        mass_scaling: sys.mass_scaling.clone(),
        warnings,
    })
}
