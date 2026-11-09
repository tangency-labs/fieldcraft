//! Canonical inverse of the supported reader subset, not source replay.
//! Public keyword references: RANDOM RESPONSE, PSD-DEFINITION, CORRELATION,
//! MODAL DAMPING, MATRIX GENERATE and MATRIX OUTPUT:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-randomresponse.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-matrixgenerate.htm
use super::*;
use fieldcraft_model::{random::Excitation, steps::MatrixKind};

fn refuse(step: &Step, what: &str) -> Error {
    Error::Unsupported { format: FORMAT, what: format!("step `{}`: {what}", step.name) }
}

impl Deck {
    pub(super) fn random_data(&mut self, model: &Model) -> Result<()> {
        for (index, step) in model.steps().iter().enumerate() {
            let Procedure::RandomVibration(random) = &step.procedure else { continue };
            let modal = random.modal.as_ref().ok_or_else(|| refuse(step, "random explicit-frequency sampling is not written"))?;
            let previous = index.checked_sub(1).and_then(|i| model.steps().get(i));
            if !previous.is_some_and(|s| {
                s.name == modal.source_step
                    && matches!(s.procedure, Procedure::Frequency { modes, min_hz: None, max_hz: None } if modes == modal.modes)
            }) || modal.bands.len() != 1
                || !random.frequencies.is_empty()
                || random.damping != fieldcraft_model::harmonic::Damping::default()
                || step.base_state.is_some()
            {
                return Err(refuse(
                    step,
                    "random export requires its immediately preceding unwindowed frequency step, one band and only modal damping",
                ));
            }
            let [Excitation::CorrelatedLoads { spectrum }] = random.excitations.as_slice() else {
                return Err(refuse(step, "random export requires one fully correlated spatial load spectrum"));
            };
            if spectrum.windows(2).any(|w| w[0].psd != w[1].psd) {
                return Err(refuse(step, "nonconstant random PSD interpolation is not written"));
            }
            self.line(&format!("*PSD-DEFINITION, NAME=fieldcraft-random-{index}, TYPE=FORCE"));
            for point in spectrum {
                self.numbers_line("", &[point.psd, 0.0, point.frequency])?;
            }
        }
        Ok(())
    }

    pub(super) fn random_step(
        &mut self,
        model: &Model,
        r: &Resolved,
        step: &Step,
        index: usize,
        base: &BTreeSet<String>,
        everything: &(Vec<&str>, Vec<u32>),
    ) -> Result<()> {
        let Procedure::RandomVibration(random) = &step.procedure else { return Err(refuse(step, "expected random procedure")) };
        let modal = random.modal.as_ref().ok_or_else(|| refuse(step, "missing random modal definition"))?;
        let band = modal.bands.first().ok_or_else(|| refuse(step, "missing random band"))?;
        let loads = StepLoads::new(model, r, step, everything);
        let concentrated = !loads.concentrated.is_empty();
        let distributed = !loads.pressure.is_empty() || !loads.element.is_empty() || !loads.centrifugal.is_empty();
        if concentrated && distributed {
            return Err(refuse(step, "mixed concentrated/distributed random load cases are not written"));
        }
        self.line(&format!("*STEP, NAME={}, PERTURBATION", step.name));
        self.line("*RANDOM RESPONSE");
        self.numbers_line("", &[band.lower, band.upper, band.points as f64, band.bias, if band.logarithmic { 0.0 } else { 1.0 }])?;
        self.step_boundaries(model, r, base, &step.boundary_conditions)?;
        self.line("*MODAL DAMPING, DEFINITION=MODE NUMBERS");
        for mode in 0..modal.modes as usize {
            let ratio = if modal.damping_ratios.len() == 1 { modal.damping_ratios.first() } else { modal.damping_ratios.get(mode) };
            self.numbers_line(&format!("{}, {}", mode + 1, mode + 1), &[ratio.copied().unwrap_or(0.0)])?;
        }
        self.component_loads(&loads, Written::default(), if concentrated { ", LOAD CASE=2" } else { "" })?;
        self.line(&format!("*CORRELATION, TYPE=CORRELATED, PSD=fieldcraft-random-{index}"));
        self.line(if concentrated { "2, 1.0" } else { "1, 1.0" });
        Ok(())
    }

    pub(super) fn matrix_step(
        &mut self,
        model: &Model,
        r: &Resolved,
        step: &Step,
        base: &BTreeSet<String>,
        matrices: &BTreeSet<MatrixKind>,
    ) -> Result<()> {
        if !step.loads.is_empty()
            || step.base_state.is_some()
            || step.boundary_conditions.iter().any(|name| model.boundary_condition(name).is_some_and(|bc| !holds_at_zero(bc)))
        {
            return Err(refuse(step, "matrix export requires unloaded, homogeneous linear boundaries without a base state"));
        }
        let flags = matrices
            .iter()
            .map(|kind| match kind {
                MatrixKind::Stiffness => "STIFFNESS",
                MatrixKind::Mass => "MASS",
                MatrixKind::ViscousDamping => "VISCOUS DAMPING",
                MatrixKind::StructuralDamping => "STRUCTURAL DAMPING",
            })
            .collect::<Vec<_>>()
            .join(", ");
        self.line(&format!("*STEP, NAME={}", step.name));
        self.line(&format!("*MATRIX GENERATE, {flags}"));
        self.step_boundaries(model, r, base, &step.boundary_conditions)?;
        self.line(&format!("*MATRIX OUTPUT, {flags}, FORMAT=MATRIX INPUT"));
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const MODEL: &str = "*NODE\n1,0,0,0\n*ELEMENT,TYPE=MASS,ELSET=M\n1,1\n*MASS,ELSET=M\n2\n*ELEMENT,TYPE=SPRING1,ELSET=K\n2,1\n*SPRING,ELSET=K\n1\n200\n*BOUNDARY\n1,2,3\n";
    #[test]
    fn random_and_matrix_are_canonical_steps_not_source_replay() {
        let input = format!(
            "{MODEL}*PSD-DEFINITION,NAME=P,TYPE=FORCE\n3,0,0\n3,0,10\n*STEP,NAME=Modes\n*FREQUENCY\n1\n*END STEP\n*STEP,NAME=Random,PERTURBATION\n*RANDOM RESPONSE\n0,10,20,3,1\n*MODAL DAMPING,DEFINITION=MODE NUMBERS\n1,1,.02\n*CLOAD,LOAD CASE=2\n1,1,5\n*CORRELATION,TYPE=CORRELATED,PSD=P\n2,1\n*END STEP\n*STEP,NAME=Matrices\n*MATRIX GENERATE,STIFFNESS,MASS\n*END STEP\n"
        );
        let first = read::read_analysis_limited(&input, Preset::MKgS, Limits::default()).unwrap();
        assert!(first.analysis_complete, "{:?}", first.warnings);
        let deck = write(&first.model, Preset::MKgS).unwrap();
        assert!(deck.contains("fieldcraft-random-") && deck.contains("*MATRIX OUTPUT"));
        let second = read::read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert!(second.analysis_complete, "{:?}", second.warnings);
        assert_eq!(first.model.steps().len(), second.model.steps().len());
        for (a, b) in first.model.steps().iter().zip(second.model.steps()) {
            assert_eq!(a.procedure, b.procedure);
        }
        let mut unsupported = first.model.clone();
        let mut step = unsupported.step("Random").unwrap().clone();
        if let Procedure::RandomVibration(random) = &mut step.procedure {
            random.damping.alpha = 1.0;
        }
        unsupported.replace_step(step).unwrap();
        assert!(write(&unsupported, Preset::MKgS).unwrap_err().to_string().contains("only modal damping"));
    }
}
