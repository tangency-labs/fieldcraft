//! Version and kind, not payload shape, determine the result contract.
use fieldcraft_io::{
    codec,
    results::{self, Decoded, Solved, StoredRef},
};
use fieldcraft_solver::{
    MassKind, StepResult,
    transient::{Energy, Frame, Response},
};

fn dynamics(mass: MassKind, spectrum: bool) -> Response {
    Response {
        mass,
        step: "Step".into(),
        nodes: vec![1],
        frames: vec![Frame {
            time: if spectrum { None } else { Some(0.) },
            displacement: vec![[0.; 3]],
            rotation: vec![[0.; 3]],
            velocity: if spectrum { vec![] } else { vec![[0.; 3]] },
            acceleration: vec![[0.; 3]],
            stress: None,
            kinetic_energy: if spectrum { None } else { Some(0.) },
            strain_energy: if spectrum { None } else { Some(0.) },
        }],
        spectrum,
        natural_frequencies: vec![1.],
        warnings: vec![],
        energy: if spectrum { None } else { Some(vec![Energy { damping_dissipation: 0., external_work: 0., energy_balance: 0. }]) },
    }
}
fn envelope(magic: &[u8; 8], kind: u8, payload: &[u8]) -> Vec<u8> {
    let mut bytes = magic.to_vec();
    bytes.extend(1u32.to_le_bytes());
    bytes.extend(4u64.to_le_bytes());
    bytes.extend(b"Step");
    bytes.extend([0; 32]);
    bytes.push(kind);
    bytes.extend((payload.len() as u64).to_le_bytes());
    bytes.extend(payload);
    bytes
}
#[test]
fn current_dynamic_kind6_keeps_mass_and_optional_energy() {
    for mass in [MassKind::Consistent, MassKind::Lumped] {
        for spectrum in [false, true] {
            let r = dynamics(mass, spectrum);
            let result = Solved::Implicit(StepResult::Dynamics(r.clone()));
            let bytes = results::encode(&[("Step", StoredRef { definition: [0; 32], result: result.as_ref() })]).unwrap();
            assert_eq!(&bytes[..8], b"MCRSLT06");
            assert_eq!(bytes[56], 6);
            // Payload begins after the u64 byte count: field presence, then enum's LE u32 index.
            assert_eq!(&bytes[65..70], if mass == MassKind::Consistent { &[1, 0, 0, 0, 0] } else { &[1, 1, 0, 0, 0] });
            let Decoded::Results(decoded) = results::decode(&bytes).unwrap() else { panic!() };
            assert_eq!(decoded[0].result, result);
        }
    }
}
#[test]
fn overlapping_legacy_kinds_are_refused_even_for_decodable_payloads() {
    let payload = codec::encode(&dynamics(MassKind::Lumped, false)).unwrap();
    for kind in 6..=10 {
        for content in [payload.as_slice(), &[255][..]] {
            let error = results::decode(&envelope(b"MCRSLT05", kind, content)).unwrap_err().to_string();
            assert!(error.contains("conflicting producer layouts"), "{error}");
        }
    }
}
#[test]
fn foreign_allocations_are_reserved_without_trial_decoding() {
    let payload = codec::encode(&dynamics(MassKind::Consistent, true)).unwrap();
    for kind in 7..=10 {
        let error = results::decode(&envelope(b"MCRSLT06", kind, &payload)).unwrap_err().to_string();
        assert!(error.contains("reserved for fatigue/harmonic"), "{error}");
    }
}
#[test]
fn version06_does_not_accept_the_original_type6_field_order() {
    #[derive(serde::Serialize)]
    struct Original {
        step: String,
        nodes: Vec<u32>,
        frames: Vec<Frame>,
        spectrum: bool,
        natural_frequencies: Vec<f64>,
        warnings: Vec<String>,
    }
    let r = dynamics(MassKind::Consistent, false);
    let old = Original {
        step: r.step,
        nodes: r.nodes,
        frames: r.frames,
        spectrum: r.spectrum,
        natural_frequencies: r.natural_frequencies,
        warnings: r.warnings,
    };
    let error = results::decode(&envelope(b"MCRSLT06", 6, &codec::encode(&old).unwrap())).unwrap_err().to_string();
    assert!(error.contains("damaged"), "{error}");
}
