//! The scalar of a sparse matrix: `f64` for the real factorizations and `Complex64` for the
//! complex-symmetric ones of a harmonic response (K − ω²M + iωC, with K, M and C real).
//!
//! The matrices are **complex symmetric**, A = Aᵀ, never Hermitian: every product, solve,
//! factorization and residual transposes and never conjugates. The bilinear form xᵀ y is the
//! one the factorization preserves; a conjugate would change the answer for any complex
//! damping. The trait promises only what `f64` and `Complex64` both have: arithmetic, zero and
//! one, a modulus for the pivot tests, and finiteness. It has no order and no sign, so the real
//! inertia (Sylvester's law, [`crate::ldlt::Inertia`]) is defined for `f64` alone.

use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_complex::Complex64;

/// An element of a sparse symmetric matrix and of the vectors it acts on.
pub trait Scalar:
    Copy
    + Send
    + Sync
    + std::fmt::Debug
    + Default
    + Sum
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    const ZERO: Self;
    const ONE: Self;

    /// |x|: the pivot tests compare moduli. A NaN stays NaN.
    fn modulus(self) -> f64;

    /// Both parts finite (every part of a complex value).
    fn is_finite(self) -> bool;
}

impl Scalar for f64 {
    const ZERO: f64 = 0.0;
    const ONE: f64 = 1.0;

    #[inline]
    fn modulus(self) -> f64 {
        self.abs()
    }

    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
}

impl Scalar for Complex64 {
    const ZERO: Complex64 = Complex64::new(0.0, 0.0);
    const ONE: Complex64 = Complex64::new(1.0, 0.0);

    #[inline]
    fn modulus(self) -> f64 {
        self.norm()
    }

    #[inline]
    fn is_finite(self) -> bool {
        Complex64::is_finite(self)
    }
}
