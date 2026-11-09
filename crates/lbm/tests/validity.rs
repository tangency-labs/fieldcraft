mod common;

use fieldcraft_lbm::{Backend, Fields, ForceHistory, Limits, validity::*};
use fieldcraft_model::{BoxFace, FaceCondition};

#[test]
fn injected_fields_keep_window_extrema_and_spatial_reference() {
    let mut fields = Fields {
        density: vec![0.0, 0.98, 1.02, 1.2],
        velocity: vec![[f64::NAN; 3], [0.03, 0.04, 0.0], [0.0, 0.0, 0.1], [0.0; 3]],
        temperature: vec![],
        fill: vec![],
    };
    let sample = measure(1, &fields, |i| i != 0, |i| if i == 3 { 1.2 } else { 1.0 }).unwrap();
    assert!((sample.density_min - 0.98).abs() < 1e-12);
    assert!((sample.density_max - 1.02).abs() < 1e-12);
    assert!((sample.density_departure - 0.02).abs() < 1e-12);
    assert_eq!(sample.mach_min, 0.0);
    assert!((sample.mach_max - 0.1 * 3.0_f64.sqrt()).abs() < 1e-12);
    fields.density[1] = 0.72;
    fields.velocity[2] = [0.2, 0.0, 0.0];
    let excursion = measure(2, &fields, |i| i != 0, |_| 1.0).unwrap();
    let final_sample = Sample { step: 3, ..sample.clone() };
    let report = Report::new(vec![sample, excursion, final_sample], 1, 1, &[]).unwrap();
    let extremes = report.extrema.unwrap();
    assert!((extremes.density_min - 0.72).abs() < 1e-12);
    assert!((extremes.mach_max - 0.2 * 3.0_f64.sqrt()).abs() < 1e-12);
    assert_eq!(report.within_study_envelope, Some(false));
}

#[test]
fn hostile_fields_are_refused_without_panics() {
    let mut fields = Fields { density: vec![1.0], velocity: vec![[0.0; 3]], temperature: vec![], fill: vec![] };
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(measure(0, &fields, |_| true, |_| bad).is_err());
        fields.density[0] = bad;
        assert!(measure(0, &fields, |_| true, |_| 1.0).is_err());
        fields.density[0] = 1.0;
    }
    fields.velocity[0][0] = f64::INFINITY;
    assert!(measure(0, &fields, |_| true, |_| 1.0).is_err());
    fields.velocity.clear();
    assert!(measure(0, &fields, |_| true, |_| 1.0).is_err());
}

#[test]
fn force_spikes_use_signed_mean_and_have_no_zero_mean_ratio() {
    let mut h = ForceHistory { name: "body".into(), body: None, force: vec![[1.0, 0.0, 0.0], [9.0, 0.0, 0.0], [-4.0, 0.0, 0.0]] };
    let r = ForceIndicator::from_history(&h);
    assert_eq!(r.peak_over_mean, Some(4.5));
    assert_eq!(r.largest_jump_over_mean, Some(6.5));
    h.force = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    assert_eq!(ForceIndicator::from_history(&h).peak_over_mean, None);
    h.force = vec![[f64::MAX, 0.0, 0.0], [-0.9 * f64::MAX, 0.0, 0.0]];
    assert_eq!(ForceIndicator::from_history(&h).largest_jump_over_mean, None);
}

#[test]
fn audit_endpoint_summaries_are_outside_the_study_envelope() {
    // Own FieldCraft audit on cb976117, published in #793/#797. These are endpoints,
    // not fabricated time histories; the true window extrema can only be more extreme.
    for (lo, hi, mach) in [(0.721, 1.374, 0.348), (0.750, 1.327, 0.163), (0.882, 1.143, 0.157), (0.920, 1.106, 0.261)] {
        let fields = Fields { density: vec![lo, hi], velocity: vec![[mach / 3.0_f64.sqrt(), 0.0, 0.0]; 2], temperature: vec![], fill: vec![] };
        let report = Report::new(vec![measure(0, &fields, |_| true, |_| 1.0).unwrap()], 1, 1, &[]).unwrap();
        assert_eq!(report.within_study_envelope, Some(false));
    }
    // #791's negative drag is not detectable from density/Mach alone. Its signed force
    // history must remain available; field-envelope acceptance is never validation.
    let h = ForceHistory { name: "body".into(), body: None, force: vec![[-3.84, 0.0, 0.0]] };
    assert_eq!(ForceIndicator::from_history(&h).peak_over_mean, Some(1.0));
    assert_eq!(h.force[0][0], -3.84);
}

#[test]
fn runtime_report_is_bounded_and_does_not_change_fields_or_forces() {
    let mut d = common::domain([0.0; 3], [4.0, 2.0, 2.0], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    let m = common::model(1.0, 0.1, d, common::step(600.0, 1.0, None));
    let (setup, mut monitored) = common::prepare(&m, &[]);
    monitored.set_equilibrium(|i| (1.0 + (i as f64 * 0.1).sin() * 0.001, [0.02, 0.0, 0.0]));
    let (_, mut control) = common::prepare(&m, &[]);
    control.set_equilibrium(|i| (1.0 + (i as f64 * 0.1).sin() * 0.001, [0.02, 0.0, 0.0]));
    let result = fieldcraft_lbm::run(&setup, &mut monitored, &Limits::default()).unwrap();
    control.advance(setup.steps).unwrap();
    assert_eq!(monitored.fields().unwrap(), control.fields().unwrap());
    assert_eq!(monitored.forces(), control.forces());
    let r = result.stats.validity.unwrap();
    assert_eq!(r.field_every_steps, 3);
    // Measured only on reads the run makes anyway (#830 review): here the three frames.
    let steps: Vec<u64> = r.samples.iter().map(|s| s.step).collect();
    assert_eq!(steps, result.frames.iter().map(|f| f.step).collect::<Vec<_>>());
    assert_eq!(steps, [0, 300, 600]);
    assert_eq!(r.within_study_envelope, Some(true));
}

#[test]
fn frequent_reads_are_thinned_to_the_sample_cap() {
    let mut d = common::domain([0.0; 3], [4.0, 2.0, 2.0], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    let mut p = common::step(600.0, 1.0, None);
    if let fieldcraft_model::Procedure::Fluid { frames, .. } = &mut p {
        *frames = 600;
    }
    let (setup, mut flow) = common::prepare(&common::model(1.0, 0.1, d, p), &[]);
    flow.set_equilibrium(|i| (1.0 + (i as f64 * 0.1).sin() * 0.001, [0.02, 0.0, 0.0]));
    let result = fieldcraft_lbm::run(&setup, &mut flow, &Limits::default()).unwrap();
    assert_eq!(result.frames.len(), 601);
    let r = result.stats.validity.unwrap();
    // A frame every step, one sample per 3-step interval: 0, 3, …, 600.
    assert_eq!(r.samples.len(), 201);
    assert!(r.samples.len() <= r.max_field_samples);
    assert!(r.samples.iter().all(|s| s.step % 3 == 0));
    assert_eq!(r.samples.last().unwrap().step, 600);
}

#[test]
fn reports_refuse_forged_samples_and_bad_cadence() {
    assert!(Report::new(vec![], 0, 1, &[]).is_err());
    assert!(Report::new(vec![Sample::default()], 1, 1, &[]).is_err());
    let valid = Sample { density_min: 1.0, density_max: 1.0, ..Sample::default() };
    assert!(Report::new(vec![valid.clone(); MAX_SAMPLES + 1], 1, 1, &[]).is_err());
    assert!(Report::new(vec![valid.clone(), valid.clone()], 1, 1, &[]).is_err());
    assert!(Report::new(vec![Sample { mach_max: f64::NAN, ..valid }], 1, 1, &[]).is_err());
}

#[test]
fn study_envelope_warns_but_existing_blowup_guard_still_errors() {
    let mut d = common::domain([0.0; 3], [2.0; 3], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    let model = common::model(1.0, 0.1, d, common::step(2.0, 1.0, None));
    let (setup, mut flow) = common::prepare(&model, &[]);
    flow.set_equilibrium(|_| (1.05, [0.02, 0.0, 0.0]));
    let report = fieldcraft_lbm::run(&setup, &mut flow, &Limits::default()).unwrap();
    assert_eq!(report.stats.steps, 2);
    assert_eq!(report.stats.validity.unwrap().within_study_envelope, Some(false));
    assert!(report.warnings.iter().any(|w| w.contains("study envelope")));
    let (_, mut unstable) = common::prepare(&model, &[]);
    unstable.set_equilibrium(|_| (1.0, [0.6, 0.0, 0.0]));
    assert!(matches!(fieldcraft_lbm::run(&setup, &mut unstable, &Limits::default()), Err(fieldcraft_lbm::Error::Unstable(_))));
}

#[test]
fn unusable_samples_cannot_become_a_passing_envelope() {
    let report = Report::new(vec![Sample { step: 3, invalid: true, ..Sample::default() }], 1, 1, &[]).unwrap();
    assert_eq!(report.invalid_field_samples, 1);
    assert_eq!(report.extrema, None);
    assert_eq!(report.within_study_envelope, None);
}

#[test]
fn original_audit_histories_reproduce_force_spikes() {
    for (text, count, peak, jump) in [
        (include_str!("data/validity/cylinder-late.csv"), 2001, 14.767437269143057, 18.153969107215243),
        (include_str!("data/validity/cylinder-ma100.csv"), 1001, 4.342400090183929, 3.4231096380122272),
    ] {
        let force = text
            .lines()
            .skip(1)
            .map(|line| {
                let mut values = line.split(',').skip(1).map(|v| v.parse::<f64>().unwrap());
                [values.next().unwrap(), values.next().unwrap(), values.next().unwrap()]
            })
            .collect();
        let r = ForceIndicator::from_history(&ForceHistory { name: "cylinder".into(), body: None, force });
        assert_eq!(r.samples, count);
        assert!((r.peak_over_mean.unwrap() - peak).abs() < 1e-12);
        assert!((r.largest_jump_over_mean.unwrap() - jump).abs() < 1e-12);
    }
}
