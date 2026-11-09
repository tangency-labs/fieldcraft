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
    assert_eq!(&b[..8], MAGIC);
    assert_eq!(b[53], kind);
    assert_eq!(decoded(&b), *r);
    let back = decoded(&b);
    assert_eq!(encode(&[("H", StoredRef { definition: [3; 32], result: back.as_ref() })]).unwrap(), b);
}
#[test]
fn original_zero_through_five_keep_payloads_and_read_version_05() {
    for (kind, r) in original_results().into_iter().enumerate() {
        round_trip(&r, kind as u8);
        let mut payload = vec![];
        // Explicit original payload schema, independent of the writer's kind dispatch.
        match &r {
            Solved::Implicit(StepResult::Static(x)) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Implicit(StepResult::Frequency(x)) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Implicit(StepResult::Heat(x)) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Implicit(StepResult::Harmonic(x)) => codec::encode_into(&LegacyHarmonicRef::from(x), &mut payload).unwrap(),
            Solved::Explicit(x) => codec::encode_into(x, &mut payload).unwrap(),
            Solved::Fluid(x) => codec::encode_into(x, &mut payload).unwrap(),
        }
        let old = record(LEGACY_MAGIC, kind as u8, &payload);
        assert_eq!(decoded(&old), r);
        let new = encode(&[("H", StoredRef { definition: [3; 32], result: r.as_ref() })]).unwrap();
        assert_eq!(&new[8..], &old[8..], "original kind {kind} payload changed");
    }
}
#[test]
fn allocated_seven_through_ten_round_trip_without_identity_loss() {
    let mut s = static_result();
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
fn six_is_reserved_and_old_overlapping_kinds_are_never_guessed() {
    let mut old_static = vec![];
    codec::encode_into(&(static_result(), fatigue()), &mut old_static).unwrap();
    let error = decode(&record(MAGIC, 6, &old_static)).unwrap_err().to_string();
    assert!(error.contains("reserved for linear dynamics"), "{error}");
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
