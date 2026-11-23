//! Exact geometric predicates: the sign of `orient2d`, `orient3d`, `incircle` and `insphere`
//! for every admissible input, symbolic perturbation for the Delaunay tie-breaks, and indirect
//! predicates on points constructed as line–plane and three-plane intersections.
//!
//! # Sources
//! - J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
//!   Predicates*, DCG 18:305–363, 1997: expansions (§2), adaptive evaluation (§3), the four
//!   predicates and their error bounds (§4). [`orient2d`], [`orient3d`], [`incircle`],
//!   [`insphere`] and their `*_filter` stage-A variants.
//! - O. Devillers, M. Teillaud, *Perturbations for Delaunay and weighted Delaunay 3D
//!   triangulations*, CGTA 44(3):160–168, 2011; H. Edelsbrunner, E. P. Mücke, *Simulation of
//!   Simplicity*, ACM TOG 9(1), 1990: [`insphere_sos`], [`incircle_sos`].
//! - M. Attene, *Indirect Predicates for Geometric Constructions*, CAD 126:102856, 2020 (paper
//!   only; his code is LGPL and is not read): [`Lpi`], [`Tpi`], [`GenericPoint`] and the
//!   `generic_*` predicates.
//! - L. Diazzi, D. Panozzo, A. Vaxman, M. Attene, *Constrained Delaunay Tetrahedrization: A
//!   Robust and Practical Approach*, ACM TOG 42(6):181, 2023, §4.3 (paper only; the authors'
//!   code is GPL/LGPL and is not read): the [`Lnc`] point on a segment, [`generic_insphere`],
//!   and its perturbed form [`generic_insphere_sos`] with the same lexicographic lift as
//!   [`insphere_sos`].
//!
//! # Arithmetic assumptions
//! IEEE 754 binary64, round-to-nearest-even, no extended intermediate precision, no FMA
//! contraction. Rust's `f64` gives exactly that on x86_64 (SSE2), aarch64 and
//! wasm32-unknown-unknown: rustc never fuses `a * b + c` and never uses x87 on these targets, and
//! the code uses no platform intrinsics, so every predicate returns the same sign on every
//! platform. The results are deterministic functions of the inputs.
//!
//! # Domain
//! Shewchuk's error bounds and the exactness of expansion arithmetic need no overflow and no
//! underflow (§5). Every coordinate of a [`Point2`] / [`Point3`] is therefore either zero or has
//! magnitude in [[`MIN_MAGNITUDE`], [`MAX_MAGNITUDE`]] = [2⁻¹²⁰, 2¹⁶⁰]: then every intermediate
//! value of the degree-5 `insphere` is either zero or above 2⁻⁸⁶⁰ and below 2⁸¹², all normal
//! doubles. The constructors reject anything else (NaN, infinities, out-of-range magnitudes) with
//! [`crate::Error::Invalid`]; the predicates themselves are total and never fail. Callers snap
//! tessellation noise (e.g. `1e-40`) to zero before building points.

mod adaptive;
mod dyadic;
mod eft;
mod expansion;
mod indirect;
mod interval;
mod ring;
mod sos;

pub use adaptive::{incircle, incircle_filter, insphere, insphere_filter, orient2d, orient2d_filter, orient3d, orient3d_filter};
pub use indirect::{Axis, GenericPoint, Lnc, Lpi, Tpi, generic_compare, generic_insphere, generic_orient2d, generic_orient3d, orient2d_projected};
pub use sos::{generic_insphere_sos, generic_lex, incircle_sos, insphere_sos};

use crate::{Error, Result};
use dyadic::Dyadic;
use ring::{Ring, Signed};

/// The orientation of three finite points, including coordinates outside the expansion domain.
/// Use the adaptive predicate when its domain accepts the points and exact dyadics otherwise.
pub fn orient2d_finite(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Result<Sign> {
    if let (Ok(a), Ok(b), Ok(c)) = (Point2::try_from(a), Point2::try_from(b), Point2::try_from(c)) {
        return Ok(orient2d(&a, &b, &c));
    }
    if !a.into_iter().chain(b).chain(c).all(f64::is_finite) {
        return Err(Error::Invalid("orientation requires finite point coordinates".into()));
    }
    let [ax, ay] = a.map(Dyadic::from_f64);
    let [bx, by] = b.map(Dyadic::from_f64);
    let [cx, cy] = c.map(Dyadic::from_f64);
    let det = (bx - ax.clone()) * (cy - ay.clone()) - (by - ay) * (cx - ax);
    Ok(det.sign().unwrap_or(Sign::Zero))
}

/// The exact sign of the signed volume in FieldCraft's tetrahedron node order.
pub fn tet_volume_sign_finite(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Result<Sign> {
    let points = [a, b, c, d];
    if let [Ok(a), Ok(b), Ok(c), Ok(d)] = points.map(Point3::try_from) {
        return Ok(tet_volume_sign(&a, &b, &c, &d));
    }
    if !points.as_flattened().iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("orientation requires finite point coordinates".into()));
    }
    let origin = a.map(Dyadic::from_f64);
    let rows = [b, c, d].map(|p| std::array::from_fn(|i| Dyadic::from_f64(p[i]) - origin[i].clone()));
    Ok(ring::det3(&rows).sign().unwrap_or(Sign::Zero))
}

/// The smallest nonzero coordinate magnitude a point may have: 2⁻¹²⁰ ≈ 7.5e-37.
pub const MIN_MAGNITUDE: f64 = f64::from_bits((1023 - 120) << 52);
/// The largest coordinate magnitude a point may have: 2¹⁶⁰ ≈ 1.5e48.
pub const MAX_MAGNITUDE: f64 = f64::from_bits((1023 + 160) << 52);

/// The sign of a predicate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Sign {
    Negative = -1,
    Zero = 0,
    Positive = 1,
}

impl Sign {
    /// The sign of a double (`Zero` for both zeros; NaN, which the domain excludes, maps to
    /// `Zero`).
    #[inline(always)]
    pub fn of(x: f64) -> Sign {
        if x > 0.0 {
            Sign::Positive
        } else if x < 0.0 {
            Sign::Negative
        } else {
            Sign::Zero
        }
    }

    #[inline(always)]
    pub fn is_zero(self) -> bool {
        self == Sign::Zero
    }
}

impl std::ops::Neg for Sign {
    type Output = Sign;

    #[inline(always)]
    fn neg(self) -> Sign {
        match self {
            Sign::Negative => Sign::Positive,
            Sign::Zero => Sign::Zero,
            Sign::Positive => Sign::Negative,
        }
    }
}

impl std::ops::Mul for Sign {
    type Output = Sign;

    #[inline(always)]
    fn mul(self, rhs: Sign) -> Sign {
        match (self, rhs) {
            (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
            (a, b) if a == b => Sign::Positive,
            _ => Sign::Negative,
        }
    }
}

/// The sign of the volume `det[b − a, c − a, d − a]`: positive when `d` lies on the side of
/// face `a, b, c`'s right-hand normal (the tet4 node order). It is [`orient3d`] negated.
#[inline]
pub fn tet_volume_sign(a: &Point3, b: &Point3, c: &Point3, d: &Point3) -> Sign {
    -orient3d(a, b, c, d)
}

/// Whether `x` may be a point coordinate (see the module docs).
#[inline]
pub fn is_admissible(x: f64) -> bool {
    x == 0.0 || (x.abs() >= MIN_MAGNITUDE && x.abs() <= MAX_MAGNITUDE)
}

/// `x` with a magnitude below [`MIN_MAGNITUDE`] rounded to zero (tessellation and construction
/// noise); other values unchanged. The result is admissible when `|x| <= MAX_MAGNITUDE`.
#[inline]
pub fn snap_to_domain(x: f64) -> f64 {
    if x.abs() < MIN_MAGNITUDE { 0.0 } else { x }
}

fn check<const N: usize>(c: [f64; N]) -> Result<[f64; N]> {
    match c.iter().find(|x| !is_admissible(**x)) {
        None => Ok(c),
        Some(x) => {
            Err(Error::Invalid(format!("point coordinate {x:e} is outside the exact-predicate domain (zero, or magnitude in [2^-120, 2^160])")))
        }
    }
}

/// A 2D point whose coordinates are admissible for the exact predicates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2([f64; 2]);

impl Point2 {
    pub fn new(x: f64, y: f64) -> Result<Self> {
        check([x, y]).map(Self)
    }

    #[inline(always)]
    pub fn coords(&self) -> [f64; 2] {
        self.0
    }
}

impl TryFrom<[f64; 2]> for Point2 {
    type Error = Error;

    fn try_from(c: [f64; 2]) -> Result<Self> {
        check(c).map(Self)
    }
}

/// A 3D point whose coordinates are admissible for the exact predicates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3([f64; 3]);

impl Point3 {
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self> {
        check([x, y, z]).map(Self)
    }

    #[inline(always)]
    pub fn coords(&self) -> [f64; 3] {
        self.0
    }
}

impl TryFrom<[f64; 3]> for Point3 {
    type Error = Error;

    fn try_from(c: [f64; 3]) -> Result<Self> {
        check(c).map(Self)
    }
}

#[cfg(test)]
mod tests;
