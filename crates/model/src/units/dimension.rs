//! Physical dimensions as integer exponents of the SI base dimensions.

use std::fmt;

/// A physical dimension: the exponents of length L, mass M, time T, temperature Θ and electric
/// current I. Dimensionless quantities (Poisson's ratio, strain, angles) have every exponent 0.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct Dimension {
    pub length: i8,
    pub mass: i8,
    pub time: i8,
    pub temperature: i8,
    pub current: i8,
}

impl Dimension {
    pub const NONE: Dimension = Dimension::new(0, 0, 0, 0);
    pub const LENGTH: Dimension = Dimension::new(1, 0, 0, 0);
    pub const MASS: Dimension = Dimension::new(0, 1, 0, 0);
    pub const TIME: Dimension = Dimension::new(0, 0, 1, 0);
    pub const TEMPERATURE: Dimension = Dimension::new(0, 0, 0, 1);
    pub const CURRENT: Dimension = Dimension { current: 1, ..Dimension::NONE };

    /// L^length M^mass T^time Θ^temperature (no current: FieldCraft has no electrical quantities yet).
    pub const fn new(length: i8, mass: i8, time: i8, temperature: i8) -> Dimension {
        Dimension { length, mass, time, temperature, current: 0 }
    }

    fn exponents(self) -> [i8; 5] {
        [self.length, self.mass, self.time, self.temperature, self.current]
    }

    fn from_exponents([length, mass, time, temperature, current]: [i8; 5]) -> Dimension {
        Dimension { length, mass, time, temperature, current }
    }

    fn zip(self, other: Dimension, f: impl Fn(i8, i8) -> Option<i8>) -> Option<Dimension> {
        let (a, b) = (self.exponents(), other.exponents());
        let mut out = [0i8; 5];
        for ((o, x), y) in out.iter_mut().zip(a).zip(b) {
            *o = f(x, y)?;
        }
        Some(Dimension::from_exponents(out))
    }

    /// The dimension of a product, or `None` if an exponent overflows.
    pub fn checked_mul(self, other: Dimension) -> Option<Dimension> {
        self.zip(other, i8::checked_add)
    }

    /// The dimension of a quotient, or `None` if an exponent overflows.
    pub fn checked_div(self, other: Dimension) -> Option<Dimension> {
        self.zip(other, i8::checked_sub)
    }

    /// The dimension raised to an integer power, or `None` if an exponent overflows.
    pub fn checked_pow(self, power: i8) -> Option<Dimension> {
        self.zip(Dimension::NONE, |x, _| x.checked_mul(power))
    }
}

impl fmt::Display for Dimension {
    /// `L M T^-2` style; `1` for dimensionless.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for (symbol, e) in ["L", "M", "T", "Θ", "I"].into_iter().zip(self.exponents()) {
            if e == 0 {
                continue;
            }
            if !first {
                f.write_str(" ")?;
            }
            first = false;
            if e == 1 { write!(f, "{symbol}")? } else { write!(f, "{symbol}^{e}")? }
        }
        if first { f.write_str("1") } else { Ok(()) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_and_display() {
        let force = Dimension::new(1, 1, -2, 0);
        let area = Dimension::LENGTH.checked_pow(2).unwrap();
        let pressure = force.checked_div(area).unwrap();
        assert_eq!(pressure, Dimension::new(-1, 1, -2, 0));
        assert_eq!(pressure.to_string(), "L^-1 M T^-2");
        assert_eq!(Dimension::NONE.to_string(), "1");
        assert_eq!(pressure.checked_mul(area), Some(force));
        assert_eq!(Dimension::new(100, 0, 0, 0).checked_pow(2), None);
        assert_eq!(Dimension::new(127, 0, 0, 0).checked_mul(Dimension::LENGTH), None);
    }
}
