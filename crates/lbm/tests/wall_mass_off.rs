//! The walls' mass correction is off unless the flow's API switches it on, and off it changes
//! nothing: the populations of a short open-channel run (the geometry of `interpolated_open.rs`,
//! 300 steps, CPU f64) hash to what the code gave before the correction existed (FNV-1a over the
//! bits of every population of every cell). Recorded on the commit before the correction; a
//! change that moves it without meaning to fails here.

mod common;

use common::*;
use fieldcraft_model::WallTreatment;

fn hash_after(walls: WallTreatment) -> u64 {
    let steps = 300u64;
    let (m, obstacles) = open_channel_model(32, steps, walls);
    let (_, mut flow) = prepare(&m, &obstacles);
    for _ in 0..steps {
        flow.step().unwrap();
    }
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for c in 0..flow.grid().len() {
        for q in 0..19 {
            for b in flow.population(q, c).to_bits().to_le_bytes() {
                h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    h
}

#[test]
fn without_the_correction_the_populations_are_what_they_were() {
    let (lin, quad) = (hash_after(WallTreatment::Interpolated), hash_after(WallTreatment::InterpolatedQuadratic));
    println!("linear {lin:#018x} quadratic {quad:#018x}");
    assert_eq!(lin, LINEAR, "linear {lin:#018x}");
    assert_eq!(quad, QUADRATIC, "quadratic {quad:#018x}");
}

const LINEAR: u64 = 0xc3326ee557b0db00;
const QUADRATIC: u64 = 0xe0f21da528905b5b;
