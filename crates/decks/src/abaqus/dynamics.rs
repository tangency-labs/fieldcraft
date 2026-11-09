//! Native dynamics export, from public DYNAMIC, MODAL DYNAMIC, RESPONSE SPECTRUM and
//! AMPLITUDE keyword references. No commercial-generated decks or results are used.
//! The target spectrum interpolates log-log, unlike our linear tables: only a common,
//! constant ordinate is exported until exact per-mode ordinates can be represented.
use super::*;
use fieldcraft_model::amplitude::Amplitude;
use fieldcraft_model::dynamics::{Combination, DirectionCombination, Integration as TimeIntegration, Method, Transient};
use fieldcraft_model::harmonic::Damping;

const NODES: &str = "fieldcraft-dynamic-nodes";
const AMPLITUDE: &str = "fieldcraft-dynamic-load";
const SPECTRUM: &str = "fieldcraft-dynamic-spectrum";
const MODES: &str = "fieldcraft-dynamic-modes";

fn refuse(step: &Step, option: impl std::fmt::Display) -> Error {
    Error::Unsupported { format: FORMAT, what: format!("step `{}`: {option}", step.name) }
}
fn common_curve<'a>(step: &Step, t: &'a Transient) -> Result<Option<&'a Amplitude>> {
    let first = step.loads.first().and_then(|name| t.load_amplitudes.get(name));
    for name in &step.loads {
        if t.load_amplitudes.get(name) != first {
            return Err(refuse(step, format!("load `{name}` has a different amplitude; separate per-load amplitudes are not written yet")));
        }
    }
    Ok(first)
}

/// Only mappings whose mass operators are known to agree are written. In particular
/// T3D2 is lumped, not a representation of FieldCraft's default consistent rod.
/// Public Abaqus Element Library / Mass Formulation and rod dynamic-loading benchmark.
pub(super) fn check_mass(model: &Model, mass: fieldcraft_model::MassKind) -> Result<()> {
    for step in model.steps().iter().filter(|s| matches!(s.procedure, Procedure::Transient(_) | Procedure::ResponseSpectrum(_))) {
        for (label, element) in model.elements() {
            if element.kind.name() != "rod2" || mass != fieldcraft_model::MassKind::Lumped {
                return Err(refuse(
                    step,
                    format!(
                        "mass={mass:?}, element {label} `{}`: native dynamic mass mapping is not supported; rod2 requires lumped (T3D2)",
                        element.kind.name()
                    ),
                ));
            }
        }
    }
    Ok(())
}

pub(super) fn check(model: &Model) -> Result<()> {
    for step in model.steps().iter().filter(|s| matches!(s.procedure, Procedure::Transient(_) | Procedure::ResponseSpectrum(_))) {
        if model.steps().len() != 1 {
            return Err(refuse(step, "independent dynamics alongside other steps; select this step alone with deck.write steps"));
        }
        if model.sets().keys().any(|n| n.eq_ignore_ascii_case(NODES)) || step.name.eq_ignore_ascii_case(MODES) {
            return Err(refuse(step, "generated dynamics node-set or extraction-step name conflicts with the model"));
        }
        for name in &step.boundary_conditions {
            if model.boundary_condition(name).is_some_and(|bc| !holds_at_zero(bc)) {
                return Err(refuse(step, format!("nonzero boundary `{name}` requires prescribed-motion export")));
            }
        }
        if !model.reference_points().is_empty() || !model.contacts().is_empty() || !model.initial_fields().is_empty() {
            return Err(refuse(step, "reference points, contact or initial fields in a linear dynamics export"));
        }
        if matches!(&step.procedure, Procedure::Transient(t) if t.adaptive.is_some()) {
            return Err(refuse(step, "adaptive transient error estimator has no equivalent native export control"));
        }
        let projected = !matches!(&step.procedure, Procedure::Transient(t) if matches!(t.method, Method::Direct));
        for (name, m) in model.materials() {
            if m.plastic.is_some() || m.hyperelastic.is_some() {
                return Err(refuse(step, format!("nonlinear material `{name}` in a linear dynamics export")));
            }
            if m.damping.is_some_and(|d| d.loss_factor != 0.0 || (projected && d != Damping::default())) {
                return Err(refuse(step, format!("material `{name}` damping needs an explicit modal projection or loss-factor convention")));
            }
        }
        for (name, f) in model.foundations() {
            if f.loss_factor != 0.0 {
                return Err(refuse(step, format!("foundation `{name}` loss factor has no time-domain convention")));
            }
        }
        for output in &step.output {
            if matches!(output, OutputVariable::Temperature | OutputVariable::HeatFlux) {
                return Err(refuse(step, format!("thermal output {output:?} in structural dynamics")));
            }
        }
        if projected {
            for (name, s) in model.springs() {
                if matches!(s, fieldcraft_model::Spring::Dashpot { .. }) {
                    return Err(refuse(step, format!("dashpot `{name}` requires nondiagonal damping projection")));
                }
            }
            for (name, f) in model.foundations() {
                if f.normal_damping != 0.0 || f.tangential_damping != 0.0 || f.loss_factor != 0.0 {
                    return Err(refuse(step, format!("foundation `{name}` damping requires modal projection")));
                }
            }
        }
        match &step.procedure {
            Procedure::Transient(t) => {
                if t.damping != Damping::default() {
                    return Err(refuse(
                        step,
                        "step Rayleigh/loss damping is not written (material and point-mass contributions cannot be conflated)",
                    ));
                }
                if t.base_acceleration.is_some() || t.base_displacement.is_some() {
                    return Err(refuse(step, "base_acceleration/base_displacement export is not implemented"));
                }
                for (kind, values) in [("initial_displacement", &t.initial_displacement), ("initial_velocity", &t.initial_velocity)] {
                    if let Some((node, _)) = values.iter().find(|(_, v)| v.iter().any(|x| *x != 0.0)) {
                        return Err(refuse(step, format!("{kind} at node {node}: nonzero initial-state export is not implemented")));
                    }
                }
                if projected && !matches!(t.integration, TimeIntegration::AverageAcceleration) {
                    return Err(refuse(step, "modal HHT/linear-acceleration integration has no MODAL DYNAMIC counterpart"));
                }
                if let Some(Amplitude::Periodic { origin, .. }) = common_curve(step, t)?
                    && *origin != 0.0
                {
                    return Err(refuse(step, "periodic amplitude with nonzero origin: target pre-origin behavior differs"));
                }
            }
            Procedure::ResponseSpectrum(s) => {
                if !step.loads.is_empty() {
                    return Err(refuse(step, "response-spectrum spatial loads"));
                }
                if !s.tables.is_empty()
                    || !s.period_tables.is_empty()
                    || s.mode_damping(0) != s.damping_ratio
                    || s.modal_damping_ratios.iter().any(|z| *z != s.damping_ratio)
                {
                    return Err(refuse(step, "response-spectrum damping/period tables or unequal modal damping: interpolation conventions differ"));
                }
                let first = s.pseudo_acceleration.first().ok_or_else(|| refuse(step, "empty pseudo_acceleration table"))?;
                if first[1] <= 0.0 || s.pseudo_acceleration.iter().any(|p| p[1] != first[1]) {
                    return Err(refuse(
                        step,
                        "nonpositive/nonconstant response spectrum: Abaqus log-log interpolation differs from linear interpolation",
                    ));
                }
            }
            _ => {}
        }
    }
    Ok(())
}

impl Deck {
    pub(super) fn dynamics_data(&mut self, model: &Model) -> Result<()> {
        let Some(step) = model.steps().first().filter(|s| matches!(s.procedure, Procedure::Transient(_) | Procedure::ResponseSpectrum(_))) else {
            return Ok(());
        };
        self.line(&format!("*NSET, NSET={NODES}"));
        // One label per line also works for arbitrarily sparse model-wide labels.
        for node in model.nodes().keys() {
            self.line(&node.to_string());
        }
        match &step.procedure {
            Procedure::Transient(t) => {
                match common_curve(step, t)? {
                    Some(Amplitude::Tabular { points } | Amplitude::SmoothStep { points }) => {
                        let definition = if matches!(common_curve(step, t)?, Some(Amplitude::SmoothStep { .. })) { "SMOOTH STEP" } else { "TABULAR" };
                        self.line(&format!("*AMPLITUDE, NAME={AMPLITUDE}, DEFINITION={definition}, TIME=STEP TIME"));
                        for point in points {
                            self.numbers_line("", point)?;
                        }
                    }
                    Some(Amplitude::Periodic { angular_frequency, origin, offset, harmonics }) => {
                        self.line(&format!("*AMPLITUDE, NAME={AMPLITUDE}, DEFINITION=PERIODIC, TIME=STEP TIME"));
                        self.numbers_line(&harmonics.len().to_string(), &[*angular_frequency, *origin, *offset])?;
                        for chunk in harmonics.chunks(4) {
                            let values: Vec<_> = chunk.iter().flatten().copied().collect();
                            self.numbers_line("", &values)?;
                        }
                    }
                    None => {
                        // Explicit constant curve prevents a target default ramp changing f(0+).
                        self.line(&format!("*AMPLITUDE, NAME={AMPLITUDE}, DEFINITION=TABULAR, TIME=STEP TIME"));
                        self.numbers_line("", &[0.0, 1.0])?;
                        self.numbers_line("", &[t.time_step * f64::from(t.increments), 1.0])?;
                    }
                }
            }
            Procedure::ResponseSpectrum(s) => {
                self.comment("Spectrum ordinates are pseudo-acceleration; verify extracted modes lie within the declared frequency band.");
                self.line(&format!("*SPECTRUM, NAME={SPECTRUM}, TYPE=ACCELERATION"));
                for p in &s.pseudo_acceleration {
                    self.numbers_line("", &[p[1], p[0], s.damping_ratio])?;
                }
            }
            _ => {}
        }
        Ok(())
    }
    pub(super) fn dynamics_step(&mut self, model: &Model, r: &Resolved, step: &Step, everything: &(Vec<&str>, Vec<u32>)) -> Result<()> {
        let modes = match &step.procedure {
            Procedure::Transient(t) => match t.method {
                Method::Modal { modes, .. } => Some(modes),
                Method::Direct => None,
            },
            Procedure::ResponseSpectrum(s) => Some(s.modes),
            _ => return Err(refuse(step, "not a dynamics procedure")),
        };
        if let Some(n) = modes {
            self.line(&format!("*STEP, NAME={MODES}, PERTURBATION"));
            self.line("*FREQUENCY, EIGENSOLVER=LANCZOS, NORMALIZATION=MASS");
            self.line(&n.to_string());
            self.line("*END STEP");
        }
        match &step.procedure {
            Procedure::Transient(t) => {
                if let Method::Modal { modes, damping_ratios } = &t.method {
                    self.comment("Target modal integration is not Newmark; compare against the continuous response with time-step tolerance.");
                    self.line(&format!("*STEP, NAME={}, PERTURBATION", step.name));
                    self.line("*MODAL DYNAMIC, CONTINUE=NO");
                    self.numbers_line("", &[t.time_step, t.time_step * f64::from(t.increments)])?;
                    self.line("*MODAL DAMPING, DEFINITION=MODE NUMBERS");
                    for i in 0..*modes {
                        let z = damping_ratios.get(if damping_ratios.len() == 1 { 0 } else { i as usize }).copied().unwrap_or(0.0);
                        self.numbers_line(&format!("{}, {}", i + 1, i + 1), &[z])?;
                    }
                } else {
                    let (alpha, beta, gamma) = match t.integration {
                        TimeIntegration::AverageAcceleration => (0.0, 0.25, 0.5),
                        TimeIntegration::LinearAcceleration => (0.0, 1.0 / 6.0, 0.5),
                        TimeIntegration::Hht { alpha } => (alpha, (1.0 - alpha).powi(2) / 4.0, 0.5 - alpha),
                    };
                    self.line(&format!("*STEP, NAME={}, NLGEOM=NO, INC={}", step.name, t.increments));
                    self.line(&format!(
                        "*DYNAMIC, DIRECT, APPLICATION=TRANSIENT FIDELITY, ALPHA={}, BETA={}, GAMMA={}",
                        number(alpha)?,
                        number(beta)?,
                        number(gamma)?
                    ));
                    self.numbers_line("", &[t.time_step, t.time_step * f64::from(t.increments)])?;
                }
                let loads = StepLoads::new(model, r, step, everything);
                self.component_loads(&loads, Written::default(), &format!(", AMPLITUDE={AMPLITUDE}"))?;
            }
            Procedure::ResponseSpectrum(s) => {
                self.line(&format!("*STEP, NAME={}, PERTURBATION", step.name));
                let sum = match s.combination {
                    Combination::Srss => "SRSS",
                    Combination::Cqc => "CQC",
                };
                let comp = match s.direction_combination {
                    DirectionCombination::Srss => "SRSS",
                    DirectionCombination::Percent1004040 => "R40",
                };
                self.line(&format!("*RESPONSE SPECTRUM, SUM={sum}, COMP={comp}"));
                for e in s.excitations() {
                    self.numbers_line(SPECTRUM, &[e.direction[0], e.direction[1], e.direction[2], e.scale])?;
                }
                self.line("*MODAL DAMPING, DEFINITION=MODE NUMBERS");
                self.numbers_line(&format!("1, {}", s.modes), &[s.damping_ratio])?;
            }
            _ => {}
        }
        let (stride, vars) = match &step.procedure {
            Procedure::Transient(t) => (t.output_every, "U, V, A"),
            _ => (1, "U, A"),
        };
        self.line(&format!("*OUTPUT, FIELD, OP=NEW, FREQUENCY={stride}"));
        self.line(&format!("*NODE OUTPUT, NSET={NODES}"));
        self.line(&format!("{vars}{}", if step.output.contains(&OutputVariable::ReactionForce) { ", RF" } else { "" }));
        let element: Vec<_> = [(OutputVariable::Stress, "S"), (OutputVariable::Strain, "E")]
            .into_iter()
            .filter(|(v, _)| step.output.contains(v))
            .map(|(_, name)| name)
            .collect();
        if !element.is_empty() {
            self.line("*ELEMENT OUTPUT");
            self.line(&element.join(", "));
        }
        if matches!(step.procedure, Procedure::Transient(_)) {
            self.line(&format!("*OUTPUT, HISTORY, FREQUENCY={stride}"));
            self.line(&format!("*NODE OUTPUT, NSET={NODES}"));
            self.line("U, V, A");
        }
        self.line("*END STEP");
        Ok(())
    }
}
