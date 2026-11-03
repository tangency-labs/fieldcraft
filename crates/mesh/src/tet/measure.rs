//! Tetrahedron measures for the quality stage, on floating-point coordinates.
//!
//! Decisions use only `+ − × ÷ √` (correctly rounded in IEEE 754), so they are the same on every
//! platform: the minimum sine of the dihedral angles (the improvement objective, Klingner and
//! Shewchuk 2007, §2), the circumcentre and the radius–edge ratio (the refinement criterion,
//! Shewchuk 1998). Angles in degrees are for reports only.
//!
//! Sources: J. R. Shewchuk, *What is a good linear finite element? Interpolation, conditioning,
//! anisotropy, and quality measures*, 2002 (the dihedral-angle sine formula
//! `sin θ_ij = 3 V ℓ_ij / (2 A_k A_l)` and the volume–length and radius-ratio measures);
//! B. M. Klingner, J. R. Shewchuk, *Aggressive Tetrahedral Mesh Improvement*, IMR 2007, §2.

pub(crate) use crate::vec3::{P, add, cross, dist, dot, lerp, mid, norm, scale, sub};

/// Six times the signed volume (positive in tet4 order).
pub fn volume6(x: &[P; 4]) -> f64 {
    let [a, b, c, d] = *x;
    dot(sub(b, a), cross(sub(c, a), sub(d, a)))
}

/// The cross product of two edges of the face opposite vertex `k` (its length is twice the
/// face's area).
fn face_cross(x: &[P; 4], k: usize) -> P {
    let [a, b, c] = match k {
        0 => [x[1], x[2], x[3]],
        1 => [x[0], x[2], x[3]],
        2 => [x[0], x[1], x[3]],
        _ => [x[0], x[1], x[2]],
    };
    cross(sub(b, a), sub(c, a))
}

/// Twice the area of the face opposite vertex `k`.
fn face_area2(x: &[P; 4], k: usize) -> f64 {
    norm(face_cross(x, k))
}

/// The edges as `(i, j, k, l)`: edge `ij`, and the two other vertices (whose opposite faces meet
/// at the edge).
const EDGE_FACES: [(usize, usize, usize, usize); 6] = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2), (1, 2, 0, 3), (1, 3, 0, 2), (2, 3, 0, 1)];

/// The sines of the six dihedral angles, signed by orientation (all negative for an inverted
/// tetrahedron), in [`EDGE_FACES`] order; zeros for a degenerate one.
pub fn dihedral_sines(x: &[P; 4]) -> [f64; 6] {
    let v6 = volume6(x);
    let areas = [0, 1, 2, 3].map(|k| face_area2(x, k));
    EDGE_FACES.map(|(i, j, k, l)| {
        let len = dist(x[i], x[j]);
        // sin θ = 3 V ℓ / (2 A_k A_l) with V = v6 / 6 and A = area2 / 2.
        let den = areas[k] * areas[l];
        let s = v6 * len / den;
        if s.is_finite() { s } else { 0.0 }
    })
}

/// The cosines of the six dihedral angles (from the outward face normals).
pub fn dihedral_cosines(x: &[P; 4]) -> [f64; 6] {
    let normal = |k: usize| -> P {
        let f: Vec<P> = (0..4).filter(|&i| i != k).filter_map(|i| x.get(i).copied()).collect();
        match f.as_slice() {
            [a, b, c] => {
                let n = cross(sub(*b, *a), sub(*c, *a));
                // Point it away from vertex k.
                let o = x.get(k).copied().unwrap_or(*a);
                if dot(n, sub(o, *a)) > 0.0 { scale(n, -1.0) } else { n }
            }
            _ => [0.0; 3],
        }
    };
    let n = [0, 1, 2, 3].map(normal);
    EDGE_FACES.map(|(_, _, k, l)| {
        let c = -dot(n[k], n[l]) / (norm(n[k]) * norm(n[l]));
        if c.is_finite() { c.clamp(-1.0, 1.0) } else { 1.0 }
    })
}

/// The quality objective of a valid tetrahedron: the smaller of its smallest dihedral sine and
/// its weighted scaled Jacobian ([`JACOBIAN_WEIGHT`]; the regular tetrahedron scores
/// √(8/9) ≈ 0.943). The dihedral sine catches slivers and caps; it does not see needles, whose
/// dihedral angles can all be fine while a corner is nearly flat (Shewchuk 2002), and the scaled
/// Jacobian `√2 · 6V / (ℓ₁ℓ₂ℓ₃)` at the worst corner, the measure analysts check (the skill's
/// limit, [`crate::quality`]), does. 0 for a flat tetrahedron; for an inverted one the negative
/// sine of largest magnitude, so worse inversions score lower. Extremes are found on squared
/// lengths and areas, with one square root each.
pub fn quality(x: &[P; 4]) -> f64 {
    let v6 = volume6(x);
    if v6 == 0.0 || !v6.is_finite() {
        return 0.0;
    }
    // Loops rather than array maps, which were not inlined in the meshing loops; the squared
    // edge lengths serve the dihedral sines and the scaled Jacobian both.
    let mut area2sq = [0.0; 4];
    for (k, s) in area2sq.iter_mut().enumerate() {
        let a = face_cross(x, k);
        *s = dot(a, a);
    }
    let mut len2 = [0.0; 6];
    let mut pick = if v6 > 0.0 { f64::INFINITY } else { 0.0 };
    for (l2, (i, j, k, l)) in len2.iter_mut().zip(EDGE_FACES) {
        let e = sub(x[i], x[j]);
        *l2 = dot(e, e);
        let r = *l2 / (area2sq[k] * area2sq[l]);
        pick = if v6 > 0.0 { pick.min(r) } else { pick.max(r) };
    }
    let q = v6 * pick.sqrt();
    let q = if v6 > 0.0 { q.min(JACOBIAN_WEIGHT * scaled_jacobian(&len2, v6)) } else { q };
    if q.is_finite() { q } else { 0.0 }
}

/// The weight of the scaled Jacobian in [`quality`]: the meshing skill's limit, a scaled
/// Jacobian of 0.3, counts as much as a 30° dihedral angle (sin 30°, the improvement's target,
/// [`super::improve::ImproveOptions::target`]), so the Jacobian decides only for tetrahedra near
/// that limit (needles) and the dihedral sine everywhere else.
pub const JACOBIAN_WEIGHT: f64 = 0.5 / 0.3;

/// The edges at each corner, as indices in [`EDGE_FACES`]: corner c's edges to the other three
/// corners in increasing order.
const CORNER_EDGES: [[usize; 3]; 4] = [[0, 1, 2], [0, 3, 4], [1, 3, 5], [2, 4, 5]];

/// The scaled Jacobian at the worst corner of a tetrahedron with six times the volume `v6` and
/// squared edge lengths `len2` (in [`EDGE_FACES`] order), as [`crate::quality`] reports it: the
/// largest product of a corner's three squared edge lengths. A length squared does not depend on
/// the edge's direction (x_j − x_c is −(x_c − x_j), and the squares lose the sign), so each is the
/// one taken from the corner.
fn scaled_jacobian(len2: &[f64; 6], v6: f64) -> f64 {
    let mut worst = 0.0_f64;
    for edges in CORNER_EDGES {
        let mut product = 1.0;
        for e in edges {
            product *= len2[e];
        }
        worst = worst.max(product);
    }
    std::f64::consts::SQRT_2 * v6 / worst.sqrt()
}

/// The smallest and largest dihedral angle in degrees (reports only).
pub fn dihedral_range(x: &[P; 4]) -> (f64, f64) {
    let (s, c) = (dihedral_sines(x), dihedral_cosines(x));
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for (sv, cv) in s.iter().zip(c) {
        let a = crate::vec3::atan2(sv.abs(), cv).to_degrees();
        lo = lo.min(a);
        hi = hi.max(a);
    }
    (lo, hi)
}

/// The circumcentre, or `None` for a flat tetrahedron.
pub fn circumcenter(x: &[P; 4]) -> Option<P> {
    let [a, b, c, d] = *x;
    let (b, c, d) = (sub(b, a), sub(c, a), sub(d, a));
    let den = 2.0 * dot(b, cross(c, d));
    if den == 0.0 || !den.is_finite() {
        return None;
    }
    let num = add(add(scale(cross(c, d), dot(b, b)), scale(cross(d, b), dot(c, c))), scale(cross(b, c), dot(d, d)));
    let o = add(a, scale(num, 1.0 / den));
    o.iter().all(|v| v.is_finite()).then_some(o)
}

/// The circumradius (infinite for a flat tetrahedron).
pub fn circumradius(x: &[P; 4]) -> f64 {
    circumcenter(x).map_or(f64::INFINITY, |o| dist(o, x[0]))
}

/// The shortest and longest edge.
pub fn edge_range(x: &[P; 4]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = 0.0_f64;
    for (i, j, _, _) in EDGE_FACES {
        let l = x.get(i).zip(x.get(j)).map_or(0.0, |(a, b)| dist(*a, *b));
        lo = lo.min(l);
        hi = hi.max(l);
    }
    (lo, hi)
}

/// Circumradius over shortest edge (√6/4 ≈ 0.612 for the regular tetrahedron).
pub fn radius_edge(x: &[P; 4]) -> f64 {
    circumradius(x) / edge_range(x).0
}

/// Three times the inradius over the circumradius, signed (1 for the regular tetrahedron).
pub fn radius_ratio(x: &[P; 4]) -> f64 {
    let v6 = volume6(x);
    let area2: f64 = (0..4).map(|k| face_area2(x, k)).sum();
    // r = 3V / A_total = (v6 / 2) / (area2 / 2) = v6 / area2.
    let r = v6 / area2;
    let q = 3.0 * r / circumradius(x);
    if q.is_finite() { q } else { 0.0 }
}

/// The centroid.
pub fn centroid(x: &[P; 4]) -> P {
    let s = x.iter().fold([0.0; 3], |acc, p| add(acc, *p));
    scale(s, 0.25)
}

/// The circumcentre of a triangle (in its plane), or `None` when degenerate.
pub fn triangle_circumcenter(a: P, b: P, c: P) -> Option<P> {
    let (u, v) = (sub(b, a), sub(c, a));
    let w = cross(u, v);
    let den = 2.0 * dot(w, w);
    if den == 0.0 || !den.is_finite() {
        return None;
    }
    let num = add(scale(cross(w, u), dot(v, v)), scale(cross(v, w), dot(u, u)));
    let o = add(a, scale(num, 1.0 / den));
    o.iter().all(|x| x.is_finite()).then_some(o)
}

/// The unit normal of a triangle (zero when degenerate).
pub fn triangle_normal(a: P, b: P, c: P) -> P {
    let n = cross(sub(b, a), sub(c, a));
    let l = norm(n);
    if l > 0.0 && l.is_finite() { scale(n, 1.0 / l) } else { [0.0; 3] }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// [`quality`] as written before it shared the squared edge lengths: the reference its
    /// loops reorganise.
    fn quality_by_maps(x: &[P; 4]) -> f64 {
        let v6 = volume6(x);
        if v6 == 0.0 || !v6.is_finite() {
            return 0.0;
        }
        let area2sq = [0, 1, 2, 3].map(|k| {
            let a = face_cross(x, k);
            dot(a, a)
        });
        let mut pick = if v6 > 0.0 { f64::INFINITY } else { 0.0 };
        for (i, j, k, l) in EDGE_FACES {
            let e = sub(x[i], x[j]);
            let r = dot(e, e) / (area2sq[k] * area2sq[l]);
            pick = if v6 > 0.0 { pick.min(r) } else { pick.max(r) };
        }
        let corners = [(0, [1, 2, 3]), (1, [0, 2, 3]), (2, [0, 1, 3]), (3, [0, 1, 2])].map(|(c, o): (usize, [usize; 3])| {
            o.iter().fold(1.0, |m, &j| {
                let e = sub(x[j], x[c]);
                m * dot(e, e)
            })
        });
        let jacobian = std::f64::consts::SQRT_2 * v6 / corners.iter().fold(0.0_f64, |m, &p| m.max(p)).sqrt();
        let q = v6 * pick.sqrt();
        let q = if v6 > 0.0 { q.min(JACOBIAN_WEIGHT * jacobian) } else { q };
        if q.is_finite() { q } else { 0.0 }
    }

    /// The quality is the same bit for bit as before: random, inverted, flat, needle-like and
    /// axis-aligned tetrahedra (whose equal coordinates make zero differences).
    #[test]
    fn quality_is_unchanged_bit_for_bit() {
        let mut state = 3u64;
        let mut random = move || {
            state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let mut tets: Vec<[P; 4]> = vec![
            [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0; 3], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0; 3], [1e-9, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ];
        for _ in 0..500 {
            tets.push([0, 1, 2, 3].map(|_| [random(), random(), random()]));
        }
        for x in tets {
            assert_eq!(quality(&x).to_bits(), quality_by_maps(&x).to_bits(), "{x:?}");
        }
    }

    const REGULAR: [P; 4] = [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];

    fn positive(mut x: [P; 4]) -> [P; 4] {
        if volume6(&x) < 0.0 {
            x.swap(0, 1);
        }
        x
    }

    #[test]
    fn regular_tetrahedron() {
        let x = positive(REGULAR);
        assert!((quality(&x) - (8.0f64 / 9.0).sqrt()).abs() < 1e-12, "{}", quality(&x));
        assert!((radius_ratio(&x) - 1.0).abs() < 1e-12);
        assert!((radius_edge(&x) - 6f64.sqrt() / 4.0).abs() < 1e-12);
        let (lo, hi) = dihedral_range(&x);
        assert!((lo - 70.528_779).abs() < 1e-5 && (hi - lo).abs() < 1e-9, "{lo} {hi}");
        let o = circumcenter(&x).unwrap();
        assert!(norm(o) < 1e-12);
    }

    #[test]
    fn right_corner_and_inverted() {
        let x = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(volume6(&x) > 0.0);
        let (lo, hi) = dihedral_range(&x);
        assert!((hi - 90.0).abs() < 1e-9 && (lo - 54.735_610).abs() < 1e-5, "{lo} {hi}");
        let mut y = x;
        y.swap(1, 2);
        assert!(quality(&y) < 0.0);
        let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        assert_eq!(quality(&flat), 0.0);
        assert!(circumcenter(&flat).is_none());
    }

    /// A needle: every dihedral angle above 20°, yet one corner nearly flat; the objective sees
    /// it through the scaled Jacobian.
    #[test]
    fn a_needle_scores_low() {
        let x = positive([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.866, 0.0], [0.5, 0.29, 12.0]]);
        let (lo, hi) = dihedral_range(&x);
        assert!(lo > 20.0 && hi < 120.0, "{lo} {hi}");
        assert!(quality(&x) < 0.15 * JACOBIAN_WEIGHT, "{}", quality(&x));
    }

    #[test]
    fn sliver_scores_low() {
        let x = positive([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.01]]);
        assert!(quality(&x) < 0.05 && quality(&x) > 0.0, "{}", quality(&x));
        assert!(radius_ratio(&x) < 0.1);
    }
}
