#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use super::check;
use crate::quality::{GateOptions, validity_gate_elements};
use crate::tet::{TaggedFace, boundary_volume, detached_tets, face_toward, face_use, linear_gate, sorted3};
use fieldcraft_model::ElementKind;
use std::collections::BTreeSet;

fn boundary(tets: &[[u32; 4]]) -> Vec<TaggedFace> {
    let count = face_use(tets);
    tets.iter()
        .enumerate()
        .flat_map(|(i, t)| {
            (0..4).filter_map({
                let count = &count;
                move |k| (count[&sorted3(face_toward(*t, k))] == 1).then_some(TaggedFace { tet: i as u32, opposite: k as u8, face: 0 })
            })
        })
        .collect()
}

/// Assert the actual old checks independently of the new certificate: positive/admissible
/// tets, exact face multiplicities/listed boundary, shared faces and the winding-weighted volume.
fn old_checks(points: &[[f64; 3]], tets: &[[u32; 4]], boundary: &[TaggedFace]) -> f64 {
    let count = face_use(tets);
    assert!(count.values().all(|n| *n <= 2));
    assert!(detached_tets(tets, &count).is_empty());
    let listed: BTreeSet<_> = boundary.iter().map(|b| sorted3(face_toward(tets[b.tet as usize], b.opposite as usize))).collect();
    assert_eq!(listed.len(), boundary.len());
    assert_eq!(listed, count.iter().filter(|(_, n)| **n == 1).map(|(f, _)| *f).collect());
    let volume = boundary_volume(points, tets, boundary).unwrap();
    let report = validity_gate_elements(
        tets.iter().enumerate().map(|(i, t)| (i as u32 + 1, ElementKind::Tet4, t.map(|v| points[v as usize]))),
        &GateOptions { expected_volume: Some(volume), ..GateOptions::default() },
    )
    .unwrap();
    assert!(report.passed(), "{report:?}");
    volume
}

fn ring(star: bool) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let p = vec![[0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0], [-2.0, 1.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -2.0, 0.0]];
    let order = if star { [2, 4, 6, 3, 5] } else { [2, 3, 4, 5, 6] };
    let t = (0..5).map(|i| [0, 1, order[i], order[(i + 1) % 5]]).collect();
    (p, t)
}

#[test]
fn a_two_turn_star_passes_every_old_check_but_not_the_edge_certificate() {
    for star in [false, true] {
        let (p, t) = ring(star);
        let b = boundary(&t);
        let volume = old_checks(&p, &t, &b);
        let result = linear_gate(&p, &t, &b, volume);
        if star {
            let error = result.unwrap_err().to_string();
            assert!(error.contains("interior edge (0, 1)") && error.contains("12.566370614359") && error.contains("(0, 0, -1)"), "{error}");
        } else {
            assert!(result.unwrap().passed());
        }
    }
}

#[test]
fn a_boundary_gap_must_be_distinguishable_from_a_full_turn() {
    // Four sectors spanning nearly one turn; the first and last ring vertices have different
    // labels. A zero-width seam is still a listed boundary and passes every old gate check.
    // Four incidences have a 4e-8 radian budget. Exercise actual measured geometry on both
    // sides of that cutoff, rather than only manipulating the accumulation's fields.
    for gap in [1e-4, 8e-8, 1e-8, 0.0] {
        let mut p = vec![[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]];
        p.extend([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [1.0, -gap, 0.0]]);
        let t: Vec<_> = (2..6).map(|i| [0, 1, i, i + 1]).collect();
        let b = boundary(&t);
        let volume = old_checks(&p, &t, &b);
        let result = linear_gate(&p, &t, &b, volume);
        if gap >= 8e-8 {
            assert!(result.unwrap().passed());
        } else {
            let error = result.unwrap_err().to_string();
            assert!(error.contains("boundary edge (0, 1)") && error.contains("budget 4.000e-8 radians"), "gap {gap}: {error}");
        }
    }
}

#[test]
fn canonical_reduction_and_unit_aware_location_survive_order_and_scale() {
    struct Halves;
    impl crate::message::Units for Halves {
        fn value(&self, _: fieldcraft_model::units::Quantity, x: f64) -> f64 {
            x / 2.0
        }
        fn symbol(&self, _: fieldcraft_model::units::Quantity) -> String {
            "u".into()
        }
    }
    for star in [false, true] {
        let (original, original_tets) = ring(star);
        for (scale, offset) in [(1.0, [0.0; 3]), (1e-3, [2.0, -3.0, 4.0]), (1e3, [-2e6, 3e6, -4e6])] {
            let p: Vec<_> = original.iter().map(|x| std::array::from_fn(|k| x[k] * scale + offset[k])).collect();
            let mut t = original_tets.clone();
            let b = boundary(&t);
            let volume = old_checks(&p, &t, &b);
            let first = linear_gate(&p, &t, &b, volume);
            t.reverse();
            for tet in &mut t {
                tet.rotate_left(2);
            } // an even permutation keeps positive orientation
            let b = boundary(&t);
            let permuted = linear_gate(&p, &t, &b, volume);
            if star {
                assert_eq!(first.as_ref().unwrap_err(), permuted.as_ref().unwrap_err());
                let error = first.unwrap_err().shown(&Halves).to_string();
                assert!(error.contains("edge (0, 1)") && error.contains(" u"), "{error}");
            } else {
                assert!(first.unwrap().passed());
                assert!(permuted.unwrap().passed());
            }
        }
    }
}

#[test]
fn invalid_or_unmeasurable_input_is_never_certified_or_panics() {
    let (p, t) = ring(false);
    let b = boundary(&t);
    let volume = old_checks(&p, &t, &b);
    for bad in [f64::NAN, f64::INFINITY, fieldcraft_model::MAX_COORD * 2.0] {
        let mut points = p.clone();
        points[0][0] = bad;
        assert!(linear_gate(&points, &t, &b, volume).is_err());
    }
    let mut missing = t.clone();
    missing[0][0] = u32::MAX;
    assert!(linear_gate(&p, &missing, &boundary(&missing), volume).is_err());
    let mut flat = p.clone();
    flat[1] = flat[0];
    assert!(!linear_gate(&flat, &t, &b, volume).unwrap().passed());
    // The private certificate also refuses an unmeasurable incidence if its precondition is
    // violated; zero fallback values from the measure helpers must not make a gap look valid.
    assert!(check(&flat, &t, &[]).is_err());
}
