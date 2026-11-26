//! The scalar of a sparse matrix: `f64` for the real factorizations and `Complex64` for the
//! complex-symmetric ones of a harmonic response (K − ω²M + iωC, with K, M and C real).
//!
//! The matrices are **complex symmetric**, A = Aᵀ, never Hermitian: every product, solve,
//! factorization and residual transposes and never conjugates. The bilinear form xᵀ y is the
//! one the factorization preserves; a conjugate would change the answer for any complex
//! damping. The trait promises only what `f64` and `Complex64` both have: arithmetic, zero and
//! one, a modulus for the pivot tests, and finiteness. It has no order and no sign, so the real
//! inertia (Sylvester's law, [`crate::ldlt::Inertia`]) is defined for `f64` alone.

use std::cell::RefCell;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_complex::Complex64;

/// Rows packed together for the micro-tiles of the dense kernels ([`Scalar::with_packed`]).
pub const PACKED_ROWS: usize = 8;

/// Packed rows a thread keeps from one call of [`Scalar::with_packed`] to the next: 2¹⁴ of them,
/// which serve a tile over 2048 columns and are 1 MB of `f64` (2 MB of `Complex64`), the order
/// of a core's second-level cache. A larger buffer is given back when the call ends: it would
/// not stay in the cache, and taking it costs 1/wu of the work of a tile over wu columns. So
/// the memory the threads keep is at most this each.
pub const PACKED_KEPT: usize = 1 << 14;

/// The buffer of packed rows a thread keeps for one scalar type.
type Packed<T> = RefCell<Vec<[T; PACKED_ROWS]>>;

/// Why [`Scalar::with_packed`] had no buffer to give.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackedError {
    /// The memory for the rows could not be had.
    OutOfMemory,
    /// The thread's buffer is in use: a call from within a call, which no kernel makes.
    InUse,
}

/// [`Scalar::with_packed`] on a thread's buffer `kept`.
fn packed<T: Copy, R>(
    kept: &'static std::thread::LocalKey<Packed<T>>,
    rows: usize,
    zero: T,
    f: impl FnOnce(&mut [[T; PACKED_ROWS]]) -> R,
) -> Result<R, PackedError> {
    kept.with(|kept| {
        let mut buffer = kept.try_borrow_mut().map_err(|_| PackedError::InUse)?;
        let more = rows.saturating_sub(buffer.len());
        buffer.try_reserve_exact(more).map_err(|_| PackedError::OutOfMemory)?;
        buffer.resize(rows, [zero; PACKED_ROWS]);
        let result = f(&mut buffer);
        if buffer.capacity() > PACKED_KEPT {
            *buffer = Vec::new();
        }
        Ok(result)
    })
}

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

    /// Run `f` on `rows` packed rows in a buffer that the calling thread keeps from one call to
    /// the next (up to [`PACKED_KEPT`] rows), so a dense kernel's tiles take no memory of their
    /// own for them. What is in the rows when `f` gets them is of no meaning. An error, never an
    /// abort, when the memory cannot be had.
    fn with_packed<R>(rows: usize, f: impl FnOnce(&mut [[Self; PACKED_ROWS]]) -> R) -> Result<R, PackedError>;
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

    fn with_packed<R>(rows: usize, f: impl FnOnce(&mut [[f64; PACKED_ROWS]]) -> R) -> Result<R, PackedError> {
        thread_local! {
            static PACKED: Packed<f64> = const { RefCell::new(Vec::new()) };
        }
        packed(&PACKED, rows, Self::ZERO, f)
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

    fn with_packed<R>(rows: usize, f: impl FnOnce(&mut [[Complex64; PACKED_ROWS]]) -> R) -> Result<R, PackedError> {
        thread_local! {
            static PACKED: Packed<Complex64> = const { RefCell::new(Vec::new()) };
        }
        packed(&PACKED, rows, Self::ZERO, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A thread's packed rows are as many as asked for, are kept up to [`PACKED_KEPT`] and given
    /// back beyond; a call from within a call and a lack of memory are errors.
    #[test]
    fn a_thread_keeps_its_packed_rows_up_to_a_bound() {
        thread_local! {
            static KEPT: Packed<f64> = const { RefCell::new(Vec::new()) };
        }
        let capacity = || KEPT.with(|kept| kept.borrow().capacity());
        let filled = packed(&KEPT, 5, 0.0, |rows| {
            rows.fill([1.0; PACKED_ROWS]);
            rows.len()
        });
        assert_eq!(filled, Ok(5));
        assert!(capacity() >= 5);
        assert_eq!(packed(&KEPT, 3, 0.0, |rows| rows.len()), Ok(3));
        assert_eq!(packed(&KEPT, PACKED_KEPT, 0.0, |rows| rows.len()), Ok(PACKED_KEPT));
        assert_eq!(capacity(), PACKED_KEPT);
        assert_eq!(packed(&KEPT, PACKED_KEPT + 1, 0.0, |rows| rows.len()), Ok(PACKED_KEPT + 1));
        assert_eq!(capacity(), 0);
        assert_eq!(packed(&KEPT, 2, 0.0, |_| packed(&KEPT, 2, 0.0, |_| ())), Ok(Err(PackedError::InUse)));
        assert_eq!(packed(&KEPT, usize::MAX, 0.0, |_| ()), Err(PackedError::OutOfMemory));
        assert_eq!(f64::with_packed(2, |_| f64::with_packed(2, |_| ())), Ok(Err(PackedError::InUse)));
        assert_eq!(Complex64::with_packed(4, |rows| rows.len()), Ok(4));
    }
}
