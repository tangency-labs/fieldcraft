//! Stress linearisation along a line through a wall: the membrane, bending and peak parts of a
//! through-thickness stress distribution, as pressure-equipment design checks classify them.
//!
//! The method is the stress-integration procedure of ASME Boiler and Pressure Vessel Code,
//! Section VIII, Division 2, Part 5, Annex 5-A "Linearization of stress results for stress
//! classification" (2023 edition), paragraph 5-A.4 (stress integration method), for continuum
//! elements (5-A.4.1). It is described here in our own words.
//!
//! A stress classification line runs from a point A on one surface of the wall to a point B on
//! the opposite surface. Its length is the thickness t, and x measures the distance from A. The
//! stresses are expressed in a local frame whose first axis is the line itself (the
//! through-thickness direction); the second and third axes lie in the section. For each
//! component σᵢⱼ of that local tensor:
//!
//! - the **membrane** stress is the mean through the wall, σᵐᵢⱼ = (1/t) ∫₀ᵗ σᵢⱼ dx (5-A.4.1.2);
//! - the **bending** stress is the linear part with zero mean, scaled to its value at the
//!   surfaces, σᵇᵢⱼ = (6/t²) ∫₀ᵗ σᵢⱼ (t/2 − x) dx (5-A.4.1.3): +σᵇ at A and −σᵇ at B. The
//!   annex computes a bending part only for the two normal stresses in the section and the
//!   shear between them; the normal stress along the line and the two shears on its planes are
//!   carried as membrane and peak only, with no bending part;
//! - the **peak** stress at each surface is what is left of the total there once the membrane
//!   and the bending stress at that surface are taken away (5-A.4.1.4).
//!
//! The integrals are evaluated with the composite Simpson rule on an odd number of equally
//! spaced samples, which is exact for polynomial distributions up to the third degree, so a
//! constant, a linear and a quadratic distribution are linearised exactly (the tests). The
//! equivalent (von Mises) stresses of the membrane tensor, of membrane plus bending and of the
//! total at both surfaces are reported as 5-A.4.1.5 defines them; whether they meet an allowable
//! stress is left to the engineer, since allowables depend on the code edition and the material.
//!
//! Background on choosing lines and integrating through the wall: J. L. Hechmer and
//! G. L. Hollinger, *3D Stress Criteria Guidelines for Application*, WRC Bulletin 429, Welding
//! Research Council, 1998.

use fieldcraft_fem::frame::{cross, dot, sub, unit};
use fieldcraft_fem::rotation::{Mat3, rotate};

use crate::{EngineError, Result};

/// Fewest samples a linearisation integrates (Simpson needs an odd number of at least three; five
/// resolve a bending profile beyond its end values).
pub const MIN_POINTS: usize = 5;
/// Most samples a linearisation integrates.
pub const MAX_POINTS: usize = 1001;

/// The components of a symmetric tensor in the order it is stored: xx, yy, zz, xy, xz, yz.
pub const COMPONENTS: [&str; 6] = ["xx", "yy", "zz", "xy", "xz", "yz"];

/// Whether a local component has a bending part: the two normal stresses in the section (yy,
/// zz) and the shear between them (yz) do; the normal stress along the line (xx) and the shears
/// on its planes (xy, xz) do not.
const BENDING: [bool; 6] = [false, true, true, false, false, true];

/// The linearised stresses of one line, in its local frame, all tensors as (xx, yy, zz, xy, xz, yz).
#[derive(Debug, Clone, PartialEq)]
pub struct Linearised {
    /// The through-wall mean of each component.
    pub membrane: [f64; 6],
    /// The bending stress at the start of the line; the end's is its negative. Zero for the
    /// components the annex carries without a bending part.
    pub bending: [f64; 6],
    /// The sampled total stress at the start and at the end.
    pub total: [[f64; 6]; 2],
    /// The peak stress at the start and at the end: total less membrane less the bending there.
    pub peak: [[f64; 6]; 2],
    /// Von Mises of the membrane tensor.
    pub membrane_equivalent: f64,
    /// Von Mises of membrane plus bending at the start and at the end.
    pub membrane_plus_bending_equivalent: [f64; 2],
    /// Von Mises of the total at the start and at the end.
    pub total_equivalent: [f64; 2],
}

fn error(why: impl Into<String>) -> EngineError {
    EngineError::invalid("points", why)
}

/// The local frame of a line from `from` to `to`, as the rows of a rotation matrix, and the
/// line's length: x along the line; z the global axis least aligned with it (ties broken in
/// favour of global Z, then Y), made orthogonal to the line; y completing a right-handed frame
/// (z × x). A radial line in a horizontal plane therefore keeps the hoop direction as its y and
/// the vessel axis Z as its z. `None` when the points coincide or are not finite.
pub fn frame(from: [f64; 3], to: [f64; 3]) -> Option<(Mat3, f64)> {
    let d = sub(to, from);
    let x = unit(d)?;
    // Strictly smaller wins, so the first candidate (Z) keeps a tie.
    let mut axis = 2;
    for k in [1, 0] {
        if x[k].abs() < x[axis].abs() {
            axis = k;
        }
    }
    let mut global = [0.0; 3];
    global[axis] = 1.0;
    let along = dot(global, x);
    let z = unit(sub(global, x.map(|v| along * v)))?;
    Some(([x, cross(z, x), z], fieldcraft_fem::frame::norm(d)))
}

/// Linearise `samples`, the global stress tensors at equally spaced points along the line
/// from its start to its end (length `t`), an odd count from [`MIN_POINTS`] to [`MAX_POINTS`].
/// `frame` holds the local axes as rows ([`frame`]). The stresses are rotated into the frame
/// and integrated with the composite Simpson rule.
pub fn linearise(samples: &[[f64; 6]], t: f64, frame: Mat3) -> Result<Linearised> {
    let n = samples.len();
    if !(MIN_POINTS..=MAX_POINTS).contains(&n) || n.is_multiple_of(2) {
        return Err(error(format!("use an odd number of points from {MIN_POINTS} to {MAX_POINTS}")));
    }
    if !(t > 0.0 && t.is_finite()) {
        return Err(error("the line must have a finite, positive length"));
    }
    let h = t / (n - 1) as f64;
    let mut local = Vec::with_capacity(n);
    for (i, stress) in samples.iter().enumerate() {
        if stress.iter().any(|v| !v.is_finite()) {
            return Err(error(format!("the stress at point {i} is not finite")));
        }
        local.push(rotate(&frame, stress));
    }
    // Composite Simpson weights 1, 4, 2, 4, …, 4, 1 times h/3.
    let weight = |i: usize| {
        if i == 0 || i + 1 == n {
            1.0
        } else if i.is_multiple_of(2) {
            2.0
        } else {
            4.0
        }
    };
    let mut membrane = [0.0; 6];
    let mut bending = [0.0; 6];
    for (i, stress) in local.iter().enumerate() {
        let w = weight(i) * h / 3.0;
        let lever = t / 2.0 - i as f64 * h;
        for c in 0..6 {
            membrane[c] += w * stress[c];
            if BENDING[c] {
                bending[c] += w * lever * stress[c];
            }
        }
    }
    for c in 0..6 {
        membrane[c] /= t;
        bending[c] *= 6.0 / (t * t);
    }
    let [first, .., last] = local.as_slice() else { return Err(error("no samples")) };
    let total = [*first, *last];
    let sign = [1.0, -1.0];
    let mut peak = [[0.0; 6]; 2];
    let mut membrane_plus_bending = [[0.0; 6]; 2];
    for end in 0..2 {
        for c in 0..6 {
            membrane_plus_bending[end][c] = membrane[c] + sign[end] * bending[c];
            peak[end][c] = total[end][c] - membrane_plus_bending[end][c];
        }
    }
    let out = Linearised {
        membrane,
        bending,
        total,
        peak,
        membrane_equivalent: fieldcraft_fem::von_mises(&membrane),
        membrane_plus_bending_equivalent: membrane_plus_bending.map(|s| fieldcraft_fem::von_mises(&s)),
        total_equivalent: total.map(|s| fieldcraft_fem::von_mises(&s)),
    };
    let finite = out.membrane.iter().chain(&out.bending).chain(out.peak.iter().flatten()).all(|v| v.is_finite())
        && out.membrane_equivalent.is_finite()
        && out.membrane_plus_bending_equivalent.iter().chain(&out.total_equivalent).all(|v| v.is_finite());
    if !finite {
        return Err(error("the linearised stresses overflow"));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_fem::rotation::IDENTITY;

    fn samples(n: usize, t: f64, f: impl Fn(f64) -> [f64; 6]) -> Vec<[f64; 6]> {
        (0..n).map(|i| f(i as f64 * t / (n - 1) as f64)).collect()
    }

    fn about_z(degrees: f64) -> Mat3 {
        let a = degrees.to_radians();
        [[a.cos(), -a.sin(), 0.0], [a.sin(), a.cos(), 0.0], [0.0, 0.0, 1.0]]
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * (1.0 + a.abs().max(b.abs()))
    }

    #[test]
    fn constant_distribution_is_all_membrane() {
        let c = [3.0, -7.0, 11.0, 0.5, -0.25, 2.0];
        let out = linearise(&samples(5, 2.0, |_| c), 2.0, IDENTITY).unwrap();
        assert!((0..6).all(|k| close(out.membrane[k], c[k])), "{out:?}");
        assert!(out.bending.iter().all(|v| close(*v, 0.0)));
        assert!(out.peak.iter().flatten().all(|v| close(*v, 0.0)));
        assert!(close(out.membrane_equivalent, fieldcraft_fem::von_mises(&c)));
        assert!(close(out.membrane_plus_bending_equivalent[0], out.membrane_equivalent));
        assert!(close(out.total_equivalent[1], out.membrane_equivalent));
    }

    #[test]
    fn linear_distribution_is_pure_bending_in_the_section_components_only() {
        let t = 4.0;
        // σ = c (t/2 − x) in every component: membrane 0, bending the surface value c t/2 at the
        // start (its negative at the end), peak 0 where the component has a bending part and the
        // whole total where it has not.
        let c = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let out = linearise(&samples(41, t, |x| c.map(|c| c * (t / 2.0 - x))), t, IDENTITY).unwrap();
        for k in 0..6 {
            assert!(close(out.membrane[k], 0.0), "{out:?}");
            if BENDING[k] {
                assert!(close(out.bending[k], c[k] * t / 2.0), "{k}: {out:?}");
                assert!(close(out.peak[0][k], 0.0) && close(out.peak[1][k], 0.0), "{k}: {out:?}");
            } else {
                assert!(close(out.bending[k], 0.0), "{k}: {out:?}");
                assert!(close(out.peak[0][k], c[k] * t / 2.0) && close(out.peak[1][k], -c[k] * t / 2.0), "{k}: {out:?}");
            }
            assert!(close(out.total[0][k], c[k] * t / 2.0) && close(out.total[1][k], -c[k] * t / 2.0));
        }
    }

    #[test]
    fn quadratic_distribution_has_the_closed_form_peak() {
        let t = 3.0;
        // σ_yy = x²: membrane t²/3, bending (6/t²) ∫ x² (t/2 − x) dx = (6/t²)(t⁴/6 − t⁴/4) = −t²/2,
        // peak at A = 0 − t²/3 + t²/2 = t²/6, peak at B = t² − t²/3 − t²/2 = t²/6.
        let out = linearise(&samples(5, t, |x| [0.0, x * x, 0.0, 0.0, 0.0, 0.0]), t, IDENTITY).unwrap();
        assert!(close(out.membrane[1], t * t / 3.0), "{out:?}");
        assert!(close(out.bending[1], -t * t / 2.0), "{out:?}");
        assert!(close(out.peak[0][1], t * t / 6.0) && close(out.peak[1][1], t * t / 6.0), "{out:?}");
        assert!(close(out.total_equivalent[1], t * t));
        // Simpson is exact for a cubic, so more points change nothing.
        let fine = linearise(&samples(1001, t, |x| [0.0, x * x, 0.0, 0.0, 0.0, 0.0]), t, IDENTITY).unwrap();
        assert!(close(fine.bending[1], out.bending[1]) && close(fine.membrane[1], out.membrane[1]));
    }

    #[test]
    fn frame_follows_the_line_and_keeps_z_up_when_it_can() {
        let (f, t) = frame([10.0, 0.0, 0.5], [20.0, 0.0, 0.5]).unwrap();
        assert!(close(t, 10.0));
        assert_eq!(f, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        // A radial line at 30° in the horizontal plane: y is the hoop direction, z stays up.
        let (c, s) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
        let (f, _) = frame([10.0 * c, 10.0 * s, 0.0], [20.0 * c, 20.0 * s, 0.0]).unwrap();
        assert!((0..3).all(|k| close(f[0][k], [c, s, 0.0][k]) && close(f[1][k], [-s, c, 0.0][k]) && close(f[2][k], [0.0, 0.0, 1.0][k])), "{f:?}");
        // A vertical line: the least aligned axes tie, Y is chosen as z, and y = z × x = X.
        let (f, _) = frame([0.0; 3], [0.0, 0.0, 2.0]).unwrap();
        assert_eq!(f, [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        // Rows are orthonormal for a skew line.
        let (f, _) = frame([1.0, 2.0, 3.0], [-2.0, 5.0, 0.5]).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| f[i][k] * f[j][k]).sum();
                assert!(close(dot, if i == j { 1.0 } else { 0.0 }), "{f:?}");
            }
        }
        assert!(frame([1.0; 3], [1.0; 3]).is_none());
        assert!(frame([0.0; 3], [f64::NAN, 0.0, 0.0]).is_none());
        assert!(frame([-f64::MAX; 3], [f64::MAX; 3]).is_none());
    }

    #[test]
    fn rotating_the_line_and_the_stresses_together_leaves_the_local_result_unchanged() {
        let t = 10.0;
        let hoop = |x: f64| 33.0 + 13_333.0 / ((10.0 + x) * (10.0 + x));
        let global = |x: f64| [-hoop(x) + 66.0, hoop(x), 20.0, 0.0, 0.0, 0.0];
        let (f0, _) = frame([10.0, 0.0, 0.5], [20.0, 0.0, 0.5]).unwrap();
        let reference = linearise(&samples(41, t, global), t, f0).unwrap();
        let q = about_z(30.0);
        let (f1, _) = frame(fieldcraft_fem::rotation::apply(&q, &[10.0, 0.0, 0.5]), fieldcraft_fem::rotation::apply(&q, &[20.0, 0.0, 0.5])).unwrap();
        let rotated = linearise(&samples(41, t, |x| rotate(&q, &global(x))), t, f1).unwrap();
        for k in 0..6 {
            assert!(close(rotated.membrane[k], reference.membrane[k]), "{rotated:?} vs {reference:?}");
            assert!(close(rotated.bending[k], reference.bending[k]));
            assert!(close(rotated.peak[0][k], reference.peak[0][k]) && close(rotated.peak[1][k], reference.peak[1][k]));
        }
    }

    #[test]
    fn hostile_samples_are_refused() {
        let good = samples(5, 1.0, |_| [1.0; 6]);
        assert!(linearise(&good, 1.0, IDENTITY).is_ok());
        assert!(linearise(&good[..4], 1.0, IDENTITY).is_err(), "even count");
        assert!(linearise(&good[..3], 1.0, IDENTITY).is_err(), "too few");
        assert!(linearise(&samples(1003, 1.0, |_| [1.0; 6]), 1.0, IDENTITY).is_err(), "too many");
        assert!(linearise(&good, 0.0, IDENTITY).is_err());
        assert!(linearise(&good, f64::INFINITY, IDENTITY).is_err());
        let mut nan = good.clone();
        nan[1][3] = f64::NAN;
        assert!(linearise(&nan, 1.0, IDENTITY).is_err());
        assert!(linearise(&samples(5, 1.0, |_| [f64::MAX; 6]), 1.0, IDENTITY).is_err(), "overflow");
    }
}
