//! Differential tests against independent oracles: exact rationals (`num-rational`, Leibniz
//! determinants written here from the definitions), `i128` on integer grids, and georust's
//! `robust` (a port of Shewchuk's public-domain `predicates.c`).

use std::cmp::Ordering;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{FromPrimitive, Signed as _, ToPrimitive, Zero};

use super::adaptive::exact;
use super::*;

// ---------------------------------------------------------------------------------------------
// Oracles and generators
// ---------------------------------------------------------------------------------------------

fn q(x: f64) -> BigRational {
    BigRational::from_f64(x).unwrap()
}

fn rsign(x: &BigRational) -> Sign {
    if x.is_zero() {
        Sign::Zero
    } else if x.is_positive() {
        Sign::Positive
    } else {
        Sign::Negative
    }
}

/// Laplace expansion along the first row.
fn rdet(m: &[Vec<BigRational>]) -> BigRational {
    if m.len() == 1 {
        return m[0][0].clone();
    }
    let mut total = BigRational::zero();
    for j in 0..m.len() {
        let minor: Vec<Vec<BigRational>> =
            m[1..].iter().map(|row| row.iter().enumerate().filter(|(k, _)| *k != j).map(|(_, v)| v.clone()).collect()).collect();
        let term = m[0][j].clone() * rdet(&minor);
        if j % 2 == 0 { total += term } else { total -= term }
    }
    total
}

fn rdiff(p: &[f64], o: &[f64]) -> Vec<BigRational> {
    p.iter().zip(o).map(|(a, b)| q(*a) - q(*b)).collect()
}

fn lifted(mut row: Vec<BigRational>) -> Vec<BigRational> {
    let lift = row.iter().fold(BigRational::zero(), |acc, v| acc + v.clone() * v.clone());
    row.push(lift);
    row
}

fn o2(a: &Point2, b: &Point2, c: &Point2) -> Sign {
    rsign(&rdet(&[rdiff(&a.0, &c.0), rdiff(&b.0, &c.0)]))
}
fn o3(a: &Point3, b: &Point3, c: &Point3, d: &Point3) -> Sign {
    rsign(&rdet(&[rdiff(&a.0, &d.0), rdiff(&b.0, &d.0), rdiff(&c.0, &d.0)]))
}
fn ic(a: &Point2, b: &Point2, c: &Point2, d: &Point2) -> Sign {
    rsign(&rdet(&[lifted(rdiff(&a.0, &d.0)), lifted(rdiff(&b.0, &d.0)), lifted(rdiff(&c.0, &d.0))]))
}
fn is(a: &Point3, b: &Point3, c: &Point3, d: &Point3, e: &Point3) -> Sign {
    rsign(&rdet(&[lifted(rdiff(&a.0, &e.0)), lifted(rdiff(&b.0, &e.0)), lifted(rdiff(&c.0, &e.0)), lifted(rdiff(&d.0, &e.0))]))
}

/// SplitMix64: a deterministic generator for reproducible tests.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [-1, 1).
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }
    fn int(&mut self, n: i64) -> i64 {
        (self.next_u64() % (2 * n as u64 + 1)) as i64 - n
    }
    /// `x` moved by up to `k` ulps either way.
    fn nudge(&mut self, mut x: f64, k: i64) -> f64 {
        let steps = self.int(k);
        for _ in 0..steps.abs() {
            x = if steps > 0 { x.next_up() } else { x.next_down() };
        }
        x
    }
}

fn p2(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}

fn rc2(p: &Point2) -> robust::Coord<f64> {
    robust::Coord { x: p.0[0], y: p.0[1] }
}
fn rc3(p: &Point3) -> robust::Coord3D<f64> {
    robust::Coord3D { x: p.0[0], y: p.0[1], z: p.0[2] }
}

#[test]
fn finite_volume_orientation_handles_subnormals_and_zero() {
    let origin = [0.; 3];
    let x = [1., 0., 0.];
    let y = [0., 1., 0.];
    let z = [0., 0., f64::from_bits(1)];
    assert_eq!(tet_volume_sign_finite(origin, x, y, z).unwrap(), Sign::Positive);
    assert_eq!(tet_volume_sign_finite(origin, y, x, z).unwrap(), Sign::Negative);
    assert_eq!(tet_volume_sign_finite(origin, x, y, y).unwrap(), Sign::Zero);
    assert!(tet_volume_sign_finite(origin, x, y, [0., 0., f64::INFINITY]).is_err());
}

#[test]
fn finite_planar_orientation_preserves_subnormal_coordinates() {
    let a = [0., 0.];
    let b = [1., 0.];
    let c = [0., f64::from_bits(1)];
    assert_eq!(orient2d_finite(a, b, c).unwrap(), Sign::Positive);
    assert_eq!(orient2d_finite(a, c, b).unwrap(), Sign::Negative);
    assert_eq!(orient2d_finite(a, c, c).unwrap(), Sign::Zero);
    assert!(orient2d_finite([f64::INFINITY, 0.], b, c).is_err());
}

/// Every stage agrees with every oracle, and the filter, when it answers, agrees too.
fn check2(a: &Point2, b: &Point2, c: &Point2, d: &Point2) {
    let want = o2(a, b, c);
    assert_eq!(orient2d(a, b, c), want, "orient2d {a:?} {b:?} {c:?}");
    assert_eq!(exact::orient2d(a, b, c), want);
    assert_eq!(Sign::of(robust::orient2d(rc2(a), rc2(b), rc2(c))), want);
    assert!(orient2d_filter(a, b, c).is_none_or(|s| s == want));
    let want = ic(a, b, c, d);
    assert_eq!(incircle(a, b, c, d), want, "incircle {a:?} {b:?} {c:?} {d:?}");
    assert_eq!(exact::incircle(a, b, c, d), want);
    assert_eq!(Sign::of(robust::incircle(rc2(a), rc2(b), rc2(c), rc2(d))), want);
    assert!(incircle_filter(a, b, c, d).is_none_or(|s| s == want));
}

fn check3(a: &Point3, b: &Point3, c: &Point3, d: &Point3, e: &Point3) {
    let want = o3(a, b, c, d);
    assert_eq!(orient3d(a, b, c, d), want, "orient3d {a:?} {b:?} {c:?} {d:?}");
    assert_eq!(exact::orient3d(a, b, c, d), want);
    assert_eq!(Sign::of(robust::orient3d(rc3(a), rc3(b), rc3(c), rc3(d))), want);
    assert!(orient3d_filter(a, b, c, d).is_none_or(|s| s == want));
    let want = is(a, b, c, d, e);
    assert_eq!(insphere(a, b, c, d, e), want, "insphere {a:?} {b:?} {c:?} {d:?} {e:?}");
    assert_eq!(exact::insphere(a, b, c, d, e), want);
    assert_eq!(Sign::of(robust::insphere(rc3(a), rc3(b), rc3(c), rc3(d), rc3(e))), want);
    assert!(insphere_filter(a, b, c, d, e).is_none_or(|s| s == want));
}

// ---------------------------------------------------------------------------------------------
// Explicit predicates
// ---------------------------------------------------------------------------------------------

#[test]
fn random_points_match_the_oracles() {
    let mut r = Rng(1);
    for _ in 0..300 {
        let pt2 = |r: &mut Rng| p2(r.unit(), r.unit());
        let (a, b, c, d) = (pt2(&mut r), pt2(&mut r), pt2(&mut r), pt2(&mut r));
        check2(&a, &b, &c, &d);
        let pt3 = |r: &mut Rng| p3(r.unit(), r.unit(), r.unit());
        let (a, b, c, d, e) = (pt3(&mut r), pt3(&mut r), pt3(&mut r), pt3(&mut r), pt3(&mut r));
        check3(&a, &b, &c, &d, &e);
    }
}

#[test]
fn near_degenerate_points_match_the_oracles() {
    // Points on a line, plane, circle or sphere, rounded and then moved by a few ulps: the
    // filter fails and every adaptive stage is exercised.
    let mut r = Rng(2);
    for i in 0..400 {
        let k = if i % 4 == 0 { 0 } else { 2 };
        let scale = [1.0, 1e-3, 1e5, 3.0e-20][i % 4];
        let (ax, ay, bx, by) = (r.unit() * scale, r.unit() * scale, r.unit() * scale, r.unit() * scale);
        let t = r.unit();
        let c = p2(r.nudge(ax + t * (bx - ax), k), r.nudge(ay + t * (by - ay), k));
        let (a, b) = (p2(ax, ay), p2(bx, by));
        let on_circle = |r: &mut Rng| {
            let th = r.unit() * 3.2;
            p2(r.nudge(scale * th.cos() + scale * 0.25, k), r.nudge(scale * th.sin() - scale * 0.5, k))
        };
        let (ca, cb, cc, cd) = (on_circle(&mut r), on_circle(&mut r), on_circle(&mut r), on_circle(&mut r));
        check2(&a, &b, &c, &c);
        check2(&ca, &cb, &cc, &cd);

        let pt = |r: &mut Rng| [r.unit() * scale, r.unit() * scale, r.unit() * scale];
        let (a, b, c) = (pt(&mut r), pt(&mut r), pt(&mut r));
        let (s, t) = (r.unit(), r.unit());
        let d: Vec<f64> = (0..3).map(|j| r.nudge(a[j] + s * (b[j] - a[j]) + t * (c[j] - a[j]), k)).collect();
        let on_sphere = |r: &mut Rng| {
            let (th, ph) = (r.unit() * 3.2, r.unit() * 1.6);
            p3(r.nudge(scale * ph.cos() * th.cos() + scale, k), r.nudge(scale * ph.cos() * th.sin(), k), r.nudge(scale * ph.sin() - scale * 0.5, k))
        };
        let (a, b, c, d) = (p3(a[0], a[1], a[2]), p3(b[0], b[1], b[2]), p3(c[0], c[1], c[2]), p3(d[0], d[1], d[2]));
        let s: Vec<Point3> = (0..5).map(|_| on_sphere(&mut r)).collect();
        check3(&a, &b, &c, &d, &s[0]);
        check3(&s[0], &s[1], &s[2], &s[3], &s[4]);
    }
}

#[test]
fn ulp_grid_near_a_line_matches_the_oracle() {
    // The classic failure of naive orient2d (Kettner et al., "Classroom examples of robustness
    // problems"): p on a 64×64 grid of ulps around (0.5, 0.5), q = (12, 12), r = (24, 24).
    // Variants with points off the diagonal make the coordinate differences inexact.
    let ulp = 2f64.powi(-53);
    for (q2, r2) in [(p2(12.0, 12.0), p2(24.0, 24.0)), (p2(0.1, 0.3), p2(17.3, 0.7)), (p2(1e-9, 3.0), p2(0.9, 1e-7))] {
        for i in 0..64 {
            for j in 0..64 {
                let p = p2(0.5 + i as f64 * ulp, 0.5 + j as f64 * ulp);
                let want = o2(&p, &q2, &r2);
                assert_eq!(orient2d(&p, &q2, &r2), want);
                assert_eq!(orient2d(&q2, &r2, &p), want);
                assert_eq!(orient2d(&r2, &p, &q2), want);
            }
        }
    }
}

fn i128_det(m: &[Vec<i128>]) -> i128 {
    if m.len() == 1 {
        return m[0][0];
    }
    (0..m.len())
        .map(|j| {
            let minor: Vec<Vec<i128>> = m[1..].iter().map(|row| row.iter().enumerate().filter(|(k, _)| *k != j).map(|(_, v)| *v).collect()).collect();
            let t = m[0][j] * i128_det(&minor);
            if j % 2 == 0 { t } else { -t }
        })
        .sum()
}

#[test]
fn integer_grids_match_i128() {
    // Small grids are full of exact degeneracies (collinear, coplanar, cocircular,
    // cospherical); large coordinates (2^30 for orient, 2^14 for in-sphere) stress the filter.
    let mut r = Rng(3);
    for &(n, big) in &[(2i64, 1i64), (3, 1 << 30), (3, 1 << 14)] {
        for _ in 0..400 {
            let v: Vec<[i64; 3]> = (0..5).map(|_| [r.int(n) * big, r.int(n) * big, r.int(n) * big]).collect();
            let pt = |p: &[i64; 3]| p3(p[0] as f64, p[1] as f64, p[2] as f64);
            let d = |p: &[i64; 3], o: &[i64; 3]| (0..3).map(|k| (p[k] - o[k]) as i128).collect::<Vec<_>>();
            let o3 = i128_det(&[d(&v[0], &v[3]), d(&v[1], &v[3]), d(&v[2], &v[3])]);
            assert_eq!(orient3d(&pt(&v[0]), &pt(&v[1]), &pt(&v[2]), &pt(&v[3])), Sign::of(o3 as f64));
            let o2v = |p: &[i64; 3]| p2(p[0] as f64, p[1] as f64);
            let d2 = |p: &[i64; 3], o: &[i64; 3]| d(p, o)[..2].to_vec();
            let want = i128_det(&[d2(&v[0], &v[2]), d2(&v[1], &v[2])]);
            assert_eq!(orient2d(&o2v(&v[0]), &o2v(&v[1]), &o2v(&v[2])), Sign::of(want as f64));
            if big <= 1 << 14 || n == 2 {
                let lift = |mut row: Vec<i128>| {
                    row.push(row.iter().map(|x| x * x).sum());
                    row
                };
                let s = i128_det(&[lift(d(&v[0], &v[4])), lift(d(&v[1], &v[4])), lift(d(&v[2], &v[4])), lift(d(&v[3], &v[4]))]);
                assert_eq!(insphere(&pt(&v[0]), &pt(&v[1]), &pt(&v[2]), &pt(&v[3]), &pt(&v[4])), Sign::of(s as f64));
                let c = i128_det(&[lift(d2(&v[0], &v[3])), lift(d2(&v[1], &v[3])), lift(d2(&v[2], &v[3]))]);
                assert_eq!(incircle(&o2v(&v[0]), &o2v(&v[1]), &o2v(&v[2]), &o2v(&v[3])), Sign::of(c as f64));
            }
        }
    }
}

#[test]
fn extreme_magnitudes_inside_the_domain_stay_exact() {
    // Coordinates at both ends of the domain, mixed, near-degenerate: no overflow or underflow
    // may corrupt a sign.
    let mut r = Rng(4);
    let tiny = MIN_MAGNITUDE;
    let huge = MAX_MAGNITUDE;
    for i in 0..150 {
        let s = [tiny, huge, 1.0, tiny * 1e10][i % 4];
        let pick = |r: &mut Rng| {
            let m = r.unit();
            if m.abs() < 0.2 { 0.0 } else { (m * s).clamp(-huge, huge) }
        };
        let fix = |x: f64| if is_admissible(x) { x } else { 0.0 };
        let v: Vec<Point3> = (0..5)
            .map(|_| {
                let (x, y, z) = (pick(&mut r), pick(&mut r), pick(&mut r));
                p3(fix(r.nudge(x, 1)), fix(y), fix(z))
            })
            .collect();
        check3(&v[0], &v[1], &v[2], &v[3], &v[4]);
        // Collinear and coplanar with huge and tiny offsets mixed.
        let base = p3(huge / 2.0, tiny, 0.0);
        let shifted = p3(huge / 2.0, tiny * 2.0, tiny);
        check3(&base, &shifted, &v[2], &v[3], &v[4]);
    }
}

#[test]
fn permutations_follow_the_determinant() {
    let mut r = Rng(5);
    for _ in 0..200 {
        // Grid points so that zero results also occur.
        let pt = |r: &mut Rng| p3(r.int(2) as f64, r.int(2) as f64, r.int(2) as f64 * 0.5);
        let (a, b, c, d, e) = (pt(&mut r), pt(&mut r), pt(&mut r), pt(&mut r), pt(&mut r));
        let s = orient3d(&a, &b, &c, &d);
        assert_eq!(orient3d(&b, &a, &c, &d), -s);
        assert_eq!(orient3d(&a, &b, &d, &c), -s);
        assert_eq!(orient3d(&b, &c, &a, &d), s);
        assert_eq!(orient3d(&b, &c, &d, &a), -s); // a 4-cycle is odd
        let s = insphere(&a, &b, &c, &d, &e);
        assert_eq!(insphere(&b, &a, &c, &d, &e), -s);
        assert_eq!(insphere(&a, &b, &c, &e, &d), -s);
        assert_eq!(insphere(&c, &a, &b, &d, &e), s);
        let (a2, b2, c2, d2) = (p2(a.0[0], a.0[1]), p2(b.0[0], b.0[1]), p2(c.0[0], c.0[1]), p2(d.0[0], d.0[2]));
        let s = orient2d(&a2, &b2, &c2);
        assert_eq!(orient2d(&b2, &a2, &c2), -s);
        assert_eq!(orient2d(&b2, &c2, &a2), s);
        let s = incircle(&a2, &b2, &c2, &d2);
        assert_eq!(incircle(&b2, &a2, &c2, &d2), -s);
        assert_eq!(incircle(&a2, &b2, &d2, &c2), -s);
    }
}

#[test]
fn conventions_match_the_documentation() {
    let (o, x, y, z) = (p3(0.0, 0.0, 0.0), p3(1.0, 0.0, 0.0), p3(0.0, 1.0, 0.0), p3(0.0, 0.0, 1.0));
    // o, x, y counterclockwise seen from +z; a point below gives Positive.
    assert_eq!(orient3d(&o, &x, &y, &p3(0.0, 0.0, -1.0)), Sign::Positive);
    assert_eq!(orient3d(&o, &x, &y, &z), Sign::Negative);
    let tet = [o, x, y, p3(0.0, 0.0, -1.0)];
    assert_eq!(insphere(&tet[0], &tet[1], &tet[2], &tet[3], &p3(0.1, 0.1, -0.1)), Sign::Positive);
    assert_eq!(insphere(&tet[0], &tet[1], &tet[2], &tet[3], &p3(5.0, 5.0, 5.0)), Sign::Negative);
    assert_eq!(insphere(&tet[0], &tet[1], &tet[2], &tet[3], &p3(1.0, 1.0, 0.0)), Sign::Zero);
    let (a, b, c) = (p2(0.0, 0.0), p2(1.0, 0.0), p2(0.0, 1.0));
    assert_eq!(orient2d(&a, &b, &c), Sign::Positive);
    assert_eq!(incircle(&a, &b, &c, &p2(0.5, 0.5)), Sign::Positive);
    assert_eq!(incircle(&a, &b, &c, &p2(1.0, 1.0)), Sign::Zero);
}

#[test]
fn the_domain_is_enforced_at_construction() {
    for bad in [f64::NAN, f64::INFINITY, -f64::INFINITY, 1e-40, -1e-300, 5e-324, 1e49, f64::MAX] {
        assert!(Point3::new(0.0, bad, 1.0).is_err(), "{bad}");
        assert!(Point2::try_from([bad, 0.0]).is_err(), "{bad}");
    }
    for good in [0.0, -0.0, MIN_MAGNITUDE, -MAX_MAGNITUDE, 1e-36, 1e48] {
        assert!(Point3::new(good, good, good).is_ok(), "{good}");
    }
}

#[test]
fn hostile_bit_patterns_never_panic() {
    let mut r = Rng(6);
    let mut accepted = Vec::new();
    for _ in 0..20_000 {
        let c = [f64::from_bits(r.next_u64()), f64::from_bits(r.next_u64()), f64::from_bits(r.next_u64())];
        if let Ok(p) = Point3::try_from(c) {
            accepted.push(p);
        }
    }
    accepted.push(p3(0.0, 0.0, 0.0));
    for w in accepted.windows(5).take(200) {
        let s = orient3d(&w[0], &w[1], &w[2], &w[3]);
        assert_eq!(s, o3(&w[0], &w[1], &w[2], &w[3]));
        let s = insphere(&w[0], &w[1], &w[2], &w[3], &w[4]);
        assert_eq!(s, is(&w[0], &w[1], &w[2], &w[3], &w[4]));
        let _ = insphere_sos(&w[0], &w[1], &w[2], &w[3], &w[4]);
    }
}

// ---------------------------------------------------------------------------------------------
// Symbolic perturbation
// ---------------------------------------------------------------------------------------------

/// The perturbed determinant evaluated with an actual tiny ε, lifts `|p|² + ε^rank` (rank 1 for
/// the lexicographically largest point): an independent check of the cofactor signs.
fn perturbed_sign(points: &[Vec<f64>]) -> Sign {
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by(|&i, &j| points[j].partial_cmp(&points[i]).unwrap_or(Ordering::Equal));
    let eps = BigRational::new(BigInt::from(1), BigInt::from(1u64 << 40));
    let rows: Vec<Vec<BigRational>> = points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let rank = order.iter().position(|&k| k == i).unwrap() + 1;
            let mut row = lifted(p.iter().map(|x| q(*x)).collect());
            let mut e = BigRational::from_integer(BigInt::from(1));
            for _ in 0..rank {
                e *= eps.clone();
            }
            if let Some(lift) = row.last_mut() {
                *lift += e;
            }
            row.push(BigRational::from_integer(BigInt::from(1)));
            row
        })
        .collect();
    rsign(&rdet(&rows))
}

#[test]
fn perturbed_insphere_matches_an_infinitesimal_lift() {
    let mut r = Rng(7);
    let mut zeros_resolved = 0;
    for _ in 0..600 {
        let v: Vec<Vec<f64>> = (0..5).map(|_| vec![r.int(1) as f64, r.int(1) as f64, r.int(1) as f64]).collect();
        if (0..5).any(|i| (0..i).any(|j| v[i] == v[j])) {
            continue;
        }
        let p: Vec<Point3> = v.iter().map(|c| p3(c[0], c[1], c[2])).collect();
        if orient3d(&p[0], &p[1], &p[2], &p[3]).is_zero() {
            continue;
        }
        let got = insphere_sos(&p[0], &p[1], &p[2], &p[3], &p[4]);
        assert_ne!(got, Sign::Zero);
        assert_eq!(got, perturbed_sign(&v), "{v:?}");
        zeros_resolved += insphere(&p[0], &p[1], &p[2], &p[3], &p[4]).is_zero() as usize;
    }
    assert!(zeros_resolved > 20, "the grid must produce cospherical cases ({zeros_resolved})");
}

#[test]
fn perturbed_incircle_matches_an_infinitesimal_lift() {
    let mut r = Rng(8);
    let mut zeros_resolved = 0;
    for _ in 0..600 {
        let v: Vec<Vec<f64>> = (0..4).map(|_| vec![r.int(1) as f64, r.int(1) as f64]).collect();
        if (0..4).any(|i| (0..i).any(|j| v[i] == v[j])) {
            continue;
        }
        let p: Vec<Point2> = v.iter().map(|c| p2(c[0], c[1])).collect();
        if orient2d(&p[0], &p[1], &p[2]).is_zero() {
            continue;
        }
        let got = incircle_sos(&p[0], &p[1], &p[2], &p[3]);
        assert_ne!(got, Sign::Zero);
        assert_eq!(got, perturbed_sign(&v), "{v:?}");
        zeros_resolved += incircle(&p[0], &p[1], &p[2], &p[3]).is_zero() as usize;
    }
    assert!(zeros_resolved > 20, "{zeros_resolved}");
}

// ---------------------------------------------------------------------------------------------
// Indirect predicates
// ---------------------------------------------------------------------------------------------

type RPoint = [BigRational; 3];

fn rp(p: &Point3) -> RPoint {
    p.0.map(q)
}
fn rsub(a: &RPoint, b: &RPoint) -> RPoint {
    std::array::from_fn(|k| a[k].clone() - b[k].clone())
}
fn rdot(a: &RPoint, b: &RPoint) -> BigRational {
    (0..3).fold(BigRational::zero(), |acc, k| acc + a[k].clone() * b[k].clone())
}
fn rcross(a: &RPoint, b: &RPoint) -> RPoint {
    [
        a[1].clone() * b[2].clone() - a[2].clone() * b[1].clone(),
        a[2].clone() * b[0].clone() - a[0].clone() * b[2].clone(),
        a[0].clone() * b[1].clone() - a[1].clone() * b[0].clone(),
    ]
}

/// Line–plane intersection by the parametric form (independent of the determinant formula).
fn r_lpi(q1: &Point3, q2: &Point3, r: &Point3, s: &Point3, t: &Point3) -> RPoint {
    let (q1, q2, r, s, t) = (rp(q1), rp(q2), rp(r), rp(s), rp(t));
    let n = rcross(&rsub(&s, &r), &rsub(&t, &r));
    let dir = rsub(&q2, &q1);
    let u = rdot(&n, &rsub(&r, &q1)) / rdot(&n, &dir);
    std::array::from_fn(|k| q1[k].clone() + u.clone() * dir[k].clone())
}

/// Three-plane intersection by Gaussian elimination.
fn r_tpi(planes: &[[Point3; 3]; 3]) -> RPoint {
    let mut m: Vec<Vec<BigRational>> = planes
        .iter()
        .map(|[r, s, t]| {
            let (r, s, t) = (rp(r), rp(s), rp(t));
            let n = rcross(&rsub(&s, &r), &rsub(&t, &r));
            let c = rdot(&n, &r);
            vec![n[0].clone(), n[1].clone(), n[2].clone(), c]
        })
        .collect();
    for col in 0..3 {
        let piv = (col..3).find(|&i| !m[i][col].is_zero()).unwrap();
        m.swap(col, piv);
        for i in 0..3 {
            if i != col {
                let f = m[i][col].clone() / m[col][col].clone();
                let pivot_row = m[col].clone();
                for (cell, p) in m[i].iter_mut().zip(&pivot_row) {
                    *cell -= p.clone() * f.clone();
                }
            }
        }
    }
    std::array::from_fn(|k| m[k][3].clone() / m[k][k].clone())
}

fn r_orient3d(p: &[RPoint; 4]) -> Sign {
    let rows: Vec<Vec<BigRational>> = (0..3).map(|i| rsub(&p[i], &p[3]).to_vec()).collect();
    rsign(&rdet(&rows))
}

fn rand_pt(r: &mut Rng, scale: f64) -> Point3 {
    p3(r.unit() * scale, r.unit() * scale, r.unit() * scale)
}

#[test]
fn indirect_orient3d_matches_rational_constructions() {
    let mut r = Rng(9);
    let mut checked = 0;
    for i in 0..150 {
        let pts: Vec<Point3> = (0..14).map(|_| rand_pt(&mut r, 1.0)).collect();
        let Ok(lpi) = Lpi::new(&pts[0], &pts[1], &pts[2], &pts[3], &pts[4]) else { continue };
        let Ok(tpi) = Tpi::new([[&pts[5], &pts[6], &pts[7]], [&pts[8], &pts[9], &pts[10]], [&pts[11], &pts[12], &pts[13]]]) else {
            continue;
        };
        let rl = r_lpi(&pts[0], &pts[1], &pts[2], &pts[3], &pts[4]);
        let rt = r_tpi(&[[pts[5], pts[6], pts[7]], [pts[8], pts[9], pts[10]], [pts[11], pts[12], pts[13]]]);
        let g = [GenericPoint::Lpi(lpi.clone()), GenericPoint::Tpi(tpi.clone())];
        let rg = [rl.clone(), rt.clone()];
        let ex: Vec<Point3> = (0..4).map(|_| rand_pt(&mut r, 1.0)).collect();
        // Every mix of implicit and explicit points in every position.
        for mask in 0u32..81 {
            let mut gp = Vec::new();
            let mut rpv = Vec::new();
            let mut m = mask;
            for e in &ex {
                match m % 3 {
                    0 => {
                        gp.push(GenericPoint::Explicit(*e));
                        rpv.push(rp(e));
                    }
                    k => {
                        gp.push(g[k as usize - 1].clone());
                        rpv.push(rg[k as usize - 1].clone());
                    }
                }
                m /= 3;
            }
            if !(i + mask as usize).is_multiple_of(7) {
                continue; // a deterministic sample of the masks keeps the test fast
            }
            let want = r_orient3d(&[rpv[0].clone(), rpv[1].clone(), rpv[2].clone(), rpv[3].clone()]);
            assert_eq!(generic_orient3d(&gp[0], &gp[1], &gp[2], &gp[3]), want, "mask {mask}");
            checked += 1;
        }
        // The constructed points lie exactly on their defining planes.
        let lp = GenericPoint::Lpi(lpi);
        let tp = GenericPoint::Tpi(tpi);
        let e = |p: &Point3| GenericPoint::Explicit(*p);
        assert_eq!(generic_orient3d(&e(&pts[2]), &e(&pts[3]), &e(&pts[4]), &lp), Sign::Zero);
        assert_eq!(generic_orient3d(&lp, &e(&pts[2]), &e(&pts[3]), &e(&pts[4])), Sign::Zero);
        for k in 0..3 {
            let pl = &pts[5 + 3 * k..8 + 3 * k];
            assert_eq!(generic_orient3d(&e(&pl[0]), &e(&pl[1]), &tp, &e(&pl[2])), Sign::Zero);
        }
        // ... and on their line.
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            assert_eq!(generic_orient2d(&e(&pts[0]), &e(&pts[1]), &lp, axis), Sign::Zero);
        }
        // Two different constructions of the same point compare equal on every axis.
        let again = GenericPoint::Lpi(Lpi::new(&pts[1], &pts[0], &pts[3], &pts[4], &pts[2]).unwrap());
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            assert_eq!(generic_compare(&lp, &again, axis), Ordering::Equal);
            let want = rl[axis.index()].cmp(&rt[axis.index()]);
            assert_eq!(generic_compare(&lp, &tp, axis), want);
            let want = rl[axis.index()].cmp(&q(ex[0].0[axis.index()]));
            assert_eq!(generic_compare(&lp, &e(&ex[0]), axis), want);
        }
        let approx = lp.approx();
        for k in 0..3 {
            assert!((approx[k] - rl[k].to_f64().unwrap()).abs() <= 1e-9 * (1.0 + approx[k].abs()));
        }
    }
    assert!(checked > 500, "{checked}");
}

#[test]
fn indirect_orient2d_matches_rational_constructions() {
    let mut r = Rng(10);
    for _ in 0..200 {
        let pts: Vec<Point3> = (0..8).map(|_| rand_pt(&mut r, 10.0)).collect();
        let Ok(lpi) = Lpi::new(&pts[0], &pts[1], &pts[2], &pts[3], &pts[4]) else { continue };
        let rl = r_lpi(&pts[0], &pts[1], &pts[2], &pts[3], &pts[4]);
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let (u, v) = axis.projected();
            for slot in 0..3 {
                let mut g: Vec<GenericPoint> = pts[5..8].iter().map(|p| GenericPoint::Explicit(*p)).collect();
                let mut rr: Vec<RPoint> = pts[5..8].iter().map(rp).collect();
                g[slot] = GenericPoint::Lpi(lpi.clone());
                rr[slot] = rl.clone();
                let d = |p: &RPoint, o: &RPoint| vec![p[u].clone() - o[u].clone(), p[v].clone() - o[v].clone()];
                let want = rsign(&rdet(&[d(&rr[0], &rr[2]), d(&rr[1], &rr[2])]));
                assert_eq!(generic_orient2d(&g[0], &g[1], &g[2], axis), want);
            }
        }
    }
}

#[test]
fn undefined_constructions_are_rejected() {
    let o = p3(0.0, 0.0, 0.0);
    let (x, y, z) = (p3(1.0, 0.0, 0.0), p3(0.0, 1.0, 0.0), p3(0.0, 0.0, 1.0));
    // Line along x, plane z = 0 containing it: parallel.
    assert!(Lpi::new(&o, &x, &o, &x, &y).is_err());
    // Degenerate line.
    assert!(Lpi::new(&z, &z, &o, &x, &y).is_err());
    // Collinear plane points.
    assert!(Lpi::new(&o, &z, &o, &x, &p3(2.0, 0.0, 0.0)).is_err());
    // A nearly parallel line is fine: z rises by one ulp over the segment.
    assert!(Lpi::new(&p3(0.0, 0.0, -1.0), &p3(1.0, 0.0, -1.0 + f64::EPSILON), &o, &x, &y).is_ok());
    // Two parallel planes.
    let lifted = |p: &Point3| p3(p.0[0], p.0[1], 1.0);
    assert!(Tpi::new([[&o, &x, &y], [&lifted(&o), &lifted(&x), &lifted(&y)], [&o, &x, &z]]).is_err());
    let t = Tpi::new([[&o, &x, &y], [&o, &y, &z], [&o, &z, &x]]).unwrap();
    let t = GenericPoint::Tpi(t);
    for axis in [Axis::X, Axis::Y, Axis::Z] {
        assert_eq!(generic_compare(&t, &GenericPoint::Explicit(o), axis), Ordering::Equal);
    }
}

// ---------------------------------------------------------------------------------------------
// LNC points and the generic in-sphere (Diazzi et al. 2023, §4.3)
// ---------------------------------------------------------------------------------------------

fn r_lnc(a: &Point3, b: &Point3, t: f64) -> RPoint {
    let (a, b) = (rp(a), rp(b));
    std::array::from_fn(|k| a[k].clone() + q(t) * (b[k].clone() - a[k].clone()))
}

fn r_insphere(p: &[RPoint; 5]) -> Sign {
    let rows: Vec<Vec<BigRational>> = (0..4).map(|i| lifted(rsub(&p[i], &p[4]).to_vec())).collect();
    rsign(&rdet(&rows))
}

#[test]
fn lnc_points_lie_on_their_segment_in_parameter_order() {
    let mut r = Rng(11);
    for _ in 0..200 {
        let (a, b) = (rand_pt(&mut r, 10.0), rand_pt(&mut r, 10.0));
        let (t0, t1) = (r.unit().clamp(1e-9, 0.5), 0.5 + r.unit().clamp(1e-9, 0.49));
        let (p0, p1) = (GenericPoint::Lnc(Lnc::new(&a, &b, t0).unwrap()), GenericPoint::Lnc(Lnc::new(&a, &b, t1).unwrap()));
        let (ga, gb) = (GenericPoint::Explicit(a), GenericPoint::Explicit(b));
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            assert_eq!(generic_orient2d(&ga, &gb, &p0, axis), Sign::Zero);
            assert_eq!(generic_orient2d(&p0, &p1, &gb, axis), Sign::Zero);
            let k = axis.index();
            let want = r_lnc(&a, &b, t0)[k].cmp(&r_lnc(&a, &b, t1)[k]);
            assert_eq!(generic_compare(&p0, &p1, axis), want);
        }
        let exact = r_lnc(&a, &b, t0);
        let rounded = p0.approx();
        let (lo, hi) = p0.bounds();
        for k in 0..3 {
            let x = exact[k].to_f64().unwrap();
            assert!((rounded[k] - x).abs() <= 2.0 * f64::EPSILON * x.abs().max(1e-300), "{} vs {x}", rounded[k]);
            assert!(q(lo[k]) <= exact[k] && exact[k] <= q(hi[k]));
        }
    }
    let (a, b) = (p3(0.0, 0.0, 0.0), p3(1.0, 1.0, 1.0));
    for bad in [0.0, 1.0, -0.5, 1.5, f64::NAN, f64::INFINITY, 1e-300] {
        assert!(Lnc::new(&a, &b, bad).is_err(), "{bad}");
    }
    assert!(Lnc::new(&a, &a, 0.5).is_err());
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn generic_orient3d_and_insphere_match_rationals_on_lnc_points() {
    let mut r = Rng(12);
    let mut checked = 0;
    for i in 0..120 {
        let pts: Vec<Point3> = (0..9).map(|_| rand_pt(&mut r, 1.0)).collect();
        let t = r.unit().clamp(1e-6, 1.0 - 1e-6);
        let lnc = Lnc::new(&pts[0], &pts[1], t).unwrap();
        let Ok(lpi) = Lpi::new(&pts[2], &pts[3], &pts[4], &pts[5], &pts[6]) else { continue };
        let g = [GenericPoint::Lnc(lnc), GenericPoint::Lpi(lpi)];
        let rg = [r_lnc(&pts[0], &pts[1], t), r_lpi(&pts[2], &pts[3], &pts[4], &pts[5], &pts[6])];
        let ex: Vec<Point3> = (0..5).map(|_| rand_pt(&mut r, 1.0)).collect();
        for mask in 0u32..243 {
            if !(i + mask as usize).is_multiple_of(5) {
                continue;
            }
            let (mut gp, mut rpv) = (Vec::new(), Vec::new());
            let mut m = mask;
            for e in &ex {
                match m % 3 {
                    0 => {
                        gp.push(GenericPoint::Explicit(*e));
                        rpv.push(rp(e));
                    }
                    k => {
                        gp.push(g[k as usize - 1].clone());
                        rpv.push(rg[k as usize - 1].clone());
                    }
                }
                m /= 3;
            }
            let want = r_orient3d(&[rpv[0].clone(), rpv[1].clone(), rpv[2].clone(), rpv[3].clone()]);
            assert_eq!(generic_orient3d(&gp[0], &gp[1], &gp[2], &gp[3]), want, "orient mask {mask}");
            let want = r_insphere(&[rpv[0].clone(), rpv[1].clone(), rpv[2].clone(), rpv[3].clone(), rpv[4].clone()]);
            assert_eq!(generic_insphere(&gp[0], &gp[1], &gp[2], &gp[3], &gp[4]), want, "insphere mask {mask}");
            checked += 1;
        }
    }
    assert!(checked > 3000, "{checked}");
    // The explicit fast path and the formula agree on the sign convention.
    let (a, b, c, d) = (p3(0.0, 0.0, 0.0), p3(1.0, 0.0, 0.0), p3(0.0, 1.0, 0.0), p3(0.0, 0.0, 1.0));
    let inside = GenericPoint::Lnc(Lnc::new(&a, &p3(1.0, 1.0, 1.0), 0.25).unwrap());
    let e = |p: &Point3| GenericPoint::Explicit(*p);
    let shewchuk_positive = orient3d(&b, &a, &c, &d) == Sign::Positive;
    let (x, y) = if shewchuk_positive { (b, a) } else { (a, b) };
    assert_eq!(insphere(&x, &y, &c, &d, &p3(0.25, 0.25, 0.25)), Sign::Positive);
    assert_eq!(generic_insphere(&e(&x), &e(&y), &e(&c), &e(&d), &inside), Sign::Positive);
}

/// On a grid, a segment's midpoint (an LNC point) is a grid point: the perturbed in-sphere of
/// the implicit form must give the explicit form's answer, including every cospherical tie.
#[test]
fn perturbed_generic_insphere_agrees_with_explicit_points() {
    let mut r = Rng(13);
    let mut ties = 0;
    for _ in 0..800 {
        let v: Vec<[f64; 3]> = (0..5).map(|_| [r.int(1) as f64, r.int(1) as f64, r.int(1) as f64]).collect();
        if (0..5).any(|i| (0..i).any(|j| v[i] == v[j])) {
            continue;
        }
        let p: Vec<Point3> = v.iter().map(|c| p3(c[0], c[1], c[2])).collect();
        if orient3d(&p[0], &p[1], &p[2], &p[3]).is_zero() {
            continue;
        }
        // Each point as the midpoint of a segment through it (direction from the generator).
        let g: Vec<GenericPoint> = p
            .iter()
            .map(|c| {
                if r.int(1) == 0 {
                    return GenericPoint::Explicit(*c);
                }
                let d = [r.int(2) as f64, r.int(2) as f64, (r.int(2) as f64) + 7.0];
                let x = c.coords();
                let a = p3(x[0] - d[0], x[1] - d[1], x[2] - d[2]);
                let b = p3(x[0] + d[0], x[1] + d[1], x[2] + d[2]);
                GenericPoint::Lnc(Lnc::new(&a, &b, 0.5).unwrap())
            })
            .collect();
        let want = insphere_sos(&p[0], &p[1], &p[2], &p[3], &p[4]);
        assert_eq!(generic_insphere_sos(&g[0], &g[1], &g[2], &g[3], &g[4]), want, "{v:?}");
        ties += insphere(&p[0], &p[1], &p[2], &p[3], &p[4]).is_zero() as usize;
        for i in 0..5 {
            for j in 0..5 {
                assert_eq!(generic_lex(&g[i], &g[j]), v[i].partial_cmp(&v[j]).unwrap());
            }
        }
    }
    assert!(ties > 20, "{ties}");
}
