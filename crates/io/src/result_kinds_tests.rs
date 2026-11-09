//! Original synthetic result records, not solver/commercial fixtures.
use super::*;
use fieldcraft_solver::{HarmonicResult, MassKind, StepResult};
use serde_json::json;

fn fatigue() -> fieldcraft_solver::fatigue::field::FatigueField {
    fieldcraft_solver::fatigue::field::FatigueField {
        damage: vec![Some(0.25)],
        life_seconds: vec![Some(40.)],
        inputs_json: "{}".into(),
        dependencies: vec![("H".into(), [3; 32])],
    }
}
fn harmonic() -> HarmonicResult {
    HarmonicResult {
        mass: MassKind::Consistent,
        operator_diagonal_bounds: None,
        fatigue: None,
        step: "H".into(),
        nodes: vec![1],
        reference_points: vec![],
        natural_frequencies: vec![],
        residual_vectors: 0,
        stats: Default::default(),
        warnings: vec![],
        frames: vec![fieldcraft_solver::harmonic::Frame {
            frequency: 1.,
            real: vec![[1., 0., 0.]],
            imaginary: vec![[0.; 3]],
            rotation_real: vec![],
            rotation_imaginary: vec![],
            points_real: vec![],
            points_imaginary: vec![],
            stress_real: Some(vec![[2., 0., 0., 0., 0., 0.]]),
            stress_imaginary: Some(vec![[0.; 6]]),
            relative_residual: 0.,
            residual_norm_si: 0.,
            projected_residual: None,
            input_power: 0.,
            dissipated_power: 0.,
        }],
    }
}
fn static_result() -> fieldcraft_solver::StaticResult {
    serde_json::from_value(json!({"step":"H","nodes":[1],"displacement":[[1,0,0]],"rotation":[],"reaction":[[0,0,0]],"reaction_moment":[],
        "stress":[[2,0,0,0,0,0]],"von_mises":[2],"elements":[],"shells":[],"beams":[],"reference_points":[],
        "applied_total":[0,0,0],"reaction_total":[0,0,0],"stats":fieldcraft_solver::Stats::default(),"warnings":[],"temperature":null,"nonlinear":null,"contacts":[]})).unwrap()
}
fn original_results() -> Vec<Solved> {
    vec![
        Solved::Implicit(StepResult::Static(static_result())),
        Solved::Implicit(StepResult::Frequency(fieldcraft_solver::FrequencyResult {
            step:"H".into(),nodes:vec![1],reference_points:vec![],modes:vec![],total_mass:1.,stats:Default::default(),warnings:vec![],buckling:None,
        })),
        Solved::Implicit(StepResult::Heat(serde_json::from_value(json!({"step":"H","nodes":[1],"times":[0],"frames":[[300]],"flux":[],"elements":[],"reaction":[],"loads":[],"applied_total":0,"reaction_total":0,"stats":fieldcraft_solver::Stats::default(),"warnings":[]})).unwrap())),
        Solved::Explicit(ExplicitResult {
            step:"H".into(), nodes:vec![1],elements:vec![],frames:vec![],history:Default::default(),controlling:vec![],initial_increment:0.01,increments:0,
            time:0.,completed:true,total_mass:1.,energy:Default::default(),initial_kinetic:0.,eroded:vec![],contacts:vec![],mass_scaling:None,warnings:vec![],
        }),
        Solved::Fluid(serde_json::from_value(json!({"port_properties":null,"step":"H","grid":{"origin":[0,0,0],"spacing":1,"cells":[1,1,1]},"solid":[false],"open_faces":[false,false,false,false,false,false],
            "frames":[],"history_times":[],"forces":[],"heat":[],"zones":[],"warnings":[],"steady":null,"temperature_reference":0,"solids":[],"surface":null,
            "stats":{"cells":1,"fluid_cells":1,"steps":0,"time_step":1,"end_time":0,"tau_plus":1,"tau_minus":1,"magic":0.25,"collision":fieldcraft_model::Collision::default(),
            "velocity":0,"lattice_velocity":0,"mach":0,"cell_updates":0,"backend":"cpu","tau_thermal":null,"eddy_viscosity_ratio":null}})).unwrap()),
        Solved::Implicit(StepResult::Harmonic(harmonic())),
    ]
}
fn record(magic: &[u8; 8], kind: u8, payload: &[u8]) -> Vec<u8> {
    let mut b = magic.to_vec();
    if magic == MAGIC {
        b.extend_from_slice(&LAYOUT_SENTINEL.to_le_bytes());
        b.extend_from_slice(LAYOUT);
    }
    b.extend_from_slice(&1u32.to_le_bytes());
    b.extend_from_slice(&1u64.to_le_bytes());
    b.push(b'H');
    b.extend_from_slice(&[3; 32]);
    b.push(kind);
    b.extend_from_slice(&(payload.len() as u64).to_le_bytes());
    b.extend_from_slice(payload);
    b
}
fn decoded(b: &[u8]) -> Solved {
    let Decoded::Results(mut r) = decode(b).unwrap() else { panic!("wrong version") };
    assert_eq!(r.len(), 1);
    assert_eq!(r[0].definition, [3; 32]);
    r.remove(0).result
}
fn round_trip(r: &Solved, kind: u8) {
    let b = encode(&[("H", StoredRef { definition: [3; 32], result: r.as_ref() })]).unwrap();
    let file = crate::write_with(
        &Model::default(),
        &Geometry::default(),
        &[],
        &crate::JournalHeader::default(),
        &crate::View::default(),
        &crate::Extras { results: vec![("H", StoredRef { definition: [3; 32], result: r.as_ref() })], ..Default::default() },
    )
    .unwrap();
    if let Some(directory) = std::env::var_os("FIELDCRAFT_RESULT_FIXTURES") {
        let directory = std::path::PathBuf::from(directory);
        std::fs::create_dir_all(&directory).unwrap();
        let hash = crate::sha256_hex(&file);
        std::fs::write(directory.join(format!("kind-{kind}-{hash}.fieldcraft")), &file).unwrap();
    }
    let opened = crate::read(&file).unwrap();
    assert_eq!(opened.results.len(), 1);
    assert_eq!(opened.results[0].result, *r);
    assert_eq!(&b[..8], MAGIC);
    assert_eq!(b[73], kind);
    assert_eq!(decoded(&b), *r);
    let back = decoded(&b);
    assert_eq!(encode(&[("H", StoredRef { definition: [3; 32], result: back.as_ref() })]).unwrap(), b);
}
#[test]
fn original_zero_through_five_read_version_05() {
    for (kind, r) in original_results().into_iter().enumerate() {
        round_trip(&r, kind as u8);
        let mut payload = vec![];
        // Explicit original payload schema, independent of the writer's kind dispatch.
        match &r {
            Solved::Implicit(StepResult::Static(x)) => {
                let old: LegacyStatic = serde_json::from_value(serde_json::to_value(x).unwrap()).unwrap();
                codec::encode_into(&old, &mut payload).unwrap();
            }
            Solved::Implicit(StepResult::Frequency(x)) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Implicit(StepResult::Heat(x)) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Implicit(StepResult::Harmonic(x)) => codec::encode_into(&LegacyHarmonicRef::from(x), &mut payload).unwrap(),
            Solved::Implicit(StepResult::Dynamics(_)) => panic!("not an original kind"),
            Solved::Explicit(x) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Fluid(x) => codec::encode_into(x, &mut payload).unwrap(),
        }
        let old = record(LEGACY_MAGIC, kind as u8, &payload);
        assert_eq!(decoded(&old), r);
        let new = encode(&[("H", StoredRef { definition: [3; 32], result: r.as_ref() })]).unwrap();
        if kind != 0 {
            assert_eq!(&new[28..], &old[8..], "original kind {kind} payload changed");
        }
    }
}
#[test]
fn allocated_seven_through_ten_round_trip_without_identity_loss() {
    let mut s = static_result();
    s.pretension =
        vec![fieldcraft_solver::pretension::SectionResult { name: "Bolt".into(), adjustment: 0.01, section_force: 10., applied_force: 10. }];
    s.fatigue = Some(fatigue());
    round_trip(&Solved::Implicit(StepResult::Static(s)), 7);
    let mut h = harmonic();
    h.fatigue = Some(fatigue());
    round_trip(&Solved::Implicit(StepResult::Harmonic(h)), 8);
    for mass in [MassKind::Consistent, MassKind::Lumped] {
        let mut h = harmonic();
        h.mass = mass;
        h.operator_diagonal_bounds = Some([[1000., 1000.], [0.5, 0.5]]);
        round_trip(&Solved::Implicit(StepResult::Harmonic(h.clone())), 9);
        h.fatigue = Some(fatigue());
        round_trip(&Solved::Implicit(StepResult::Harmonic(h)), 10);
    }
    let mut h = harmonic();
    h.mass = MassKind::Lumped;
    round_trip(&Solved::Implicit(StepResult::Harmonic(h)), 9);
}
#[test]
fn old_overlapping_kinds_are_never_guessed() {
    let mut old_static = vec![];
    codec::encode_into(&(static_result(), fatigue()), &mut old_static).unwrap();
    let error = decode(&record(MAGIC, 6, &old_static)).unwrap_err().to_string();
    assert!(error.contains("damaged"), "{error}");
    let mut old_harmonic = vec![];
    codec::encode_into(&(harmonic(), fatigue()), &mut old_harmonic).unwrap();
    for kind in 6..=10 {
        // Even decodable fatigue bytes may not be guessed from an old overlapping number.
        for payload in [&old_static, &old_harmonic] {
            let error = decode(&record(LEGACY_MAGIC, kind, payload)).unwrap_err().to_string();
            assert!(error.contains("conflicting producer layouts"), "{kind}: {error}");
        }
    }
}

#[test]
fn unmarked_six_from_both_producers_requires_verified_recovery() {
    for kind in 0..=10 {
        let current = record(MAGIC, kind, &[]);
        let mut old = MAGIC.to_vec();
        old.extend_from_slice(&current[28..]);
        assert!(decode(&old).unwrap_err().to_string().contains("unmarked version 06"));
    }
    let mut unknown = record(MAGIC, 6, &[]);
    unknown[12] ^= 1;
    assert!(decode(&unknown).unwrap_err().to_string().contains("unknown version 06 layout"));
}

#[test]
fn six_preserves_complete_dynamics_and_energy() {
    use fieldcraft_solver::transient::{Energy, Frame, Response};
    let r = Response {
        step: "H".into(),
        nodes: vec![1],
        spectrum: false,
        natural_frequencies: vec![2.],
        warnings: vec!["sample".into()],
        frames: vec![Frame {
            time: Some(0.0),
            displacement: vec![[1., 2., 3.]],
            rotation: vec![[4., 5., 6.]],
            velocity: vec![[7., 8., 9.]],
            acceleration: vec![[10., 11., 12.]],
            stress: Some(vec![[1., 2., 3., 4., 5., 6.]]),
            kinetic_energy: Some(2.),
            strain_energy: Some(3.),
        }],
        energy: Some(vec![Energy { damping_dissipation: 4., external_work: 5., energy_balance: 6. }]),
    };
    round_trip(&Solved::Implicit(StepResult::Dynamics(r.clone())), 6);
    let mut spectrum = r;
    spectrum.frames[0].velocity.clear();
    spectrum.frames[0].kinetic_energy = None;
    spectrum.frames[0].strain_energy = None;
    spectrum.spectrum = true;
    spectrum.frames[0].time = None;
    spectrum.energy = None;
    round_trip(&Solved::Implicit(StepResult::Dynamics(spectrum)), 6);
}
