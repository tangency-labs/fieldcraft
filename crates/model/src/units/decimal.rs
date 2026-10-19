//! Scaling by powers of ten, and adding offsets, done on the decimal a value was written as.
//!
//! A number typed as `3.97` (mm) is the decimal 3.97, not the binary double nearest to it. Scaled
//! to metres in binary (`3.97 / 1000`) and back (`× 1000`) it comes out as `3.9699999999999998`
//! about once in a hundred values, and a stored `0.00397` m would differ from the double nearest
//! to the decimal 0.00397. Here the scaling shifts the exponent of the value's shortest decimal
//! form instead (Rust prints the shortest decimal that reads back as the same double), and the
//! result is the double nearest to the shifted decimal. So a value written with up to
//! [`EXACT_DIGITS`] significant digits converts to the double nearest the exact result, and back
//! to exactly what was written.
//!
//! A value with more digits (a computed coordinate, a result) is scaled in binary with a single
//! correctly rounded multiplication or division by an exact power of ten: within half an ulp of
//! the exact result, and the round trip within one ulp.

/// Significant decimal digits up to which every decimal survives the trip through a double: two
/// distinct decimals of this many digits are always two distinct doubles.
pub const EXACT_DIGITS: usize = 15;

/// The powers of ten a double holds exactly.
const POW10: [f64; 23] =
    [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22];

/// The shortest decimal form of `v` as an integer mantissa and a power of ten (`3.97` is
/// `(397, -2)`), when it has at most [`EXACT_DIGITS`] significant digits; `None` for longer ones,
/// zero and non-finite values.
pub(super) fn short_decimal(v: f64) -> Option<(i128, i32)> {
    if v == 0.0 || !v.is_finite() {
        return None;
    }
    // `{:e}` is the shortest round-trip form: `-3.97e0`, `1e-5`, `2.1e11`.
    let text = format!("{v:e}");
    let (mantissa, exponent) = text.split_once('e')?;
    let exponent: i32 = exponent.parse().ok()?;
    let (negative, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => (true, m),
        None => (false, mantissa),
    };
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = whole.len() + fraction.len();
    if digits > EXACT_DIGITS {
        return None;
    }
    let mut m: i128 = 0;
    for c in whole.chars().chain(fraction.chars()) {
        m = m.checked_mul(10)?.checked_add(i128::from(c.to_digit(10)?))?;
    }
    let fraction_len = i32::try_from(fraction.len()).ok()?;
    Some((if negative { -m } else { m }, exponent.checked_sub(fraction_len)?))
}

/// The double nearest to `mantissa × 10^exponent`.
pub(super) fn from_decimal(mantissa: i128, exponent: i32) -> f64 {
    // Rust's float parsing rounds correctly.
    format!("{mantissa}e{exponent}").parse().unwrap_or(f64::NAN)
}

/// `v × 10^k`.
pub(super) fn shift(v: f64, k: i32) -> f64 {
    if k == 0 || v == 0.0 || !v.is_finite() {
        return v;
    }
    if let Some((m, e)) = short_decimal(v)
        && let Some(e) = e.checked_add(k)
    {
        return from_decimal(m, e);
    }
    binary_shift(v, k)
}

/// `v × 10^k` in binary: one correctly rounded operation by an exact power of ten while `|k|`
/// allows, then steps of 10^22.
fn binary_shift(mut v: f64, mut k: i32) -> f64 {
    let step = POW10.len() as i32 - 1;
    while k != 0 && v.is_finite() && v != 0.0 {
        let n = k.clamp(-step, step);
        let p = POW10.get(n.unsigned_abs() as usize).copied().unwrap_or(1.0);
        v = if n > 0 { v * p } else { v / p };
        k -= n;
    }
    v
}

/// Largest integer a double holds exactly.
const EXACT_INTEGER: i128 = 1 << 53;

/// An absolute temperature `f` in °F, in kelvin: (F + 459.67) × 5/9 worked out on the decimal F
/// was written as (one correctly rounded division of two exact integers), so 68 °F is 293.15 K.
pub(super) fn fahrenheit_to_kelvin(f: f64) -> f64 {
    if let Some((m, e)) = short_decimal(f)
        && let Some(p) = u32::try_from(-e).ok().map(|p| p.max(2)).filter(|p| *p <= 6)
    {
        let scale = 10i128.pow(p);
        let n = m.checked_mul(10i128.pow(p - u32::try_from(-e).unwrap_or(0))).and_then(|m| m.checked_add(45967 * scale / 100));
        if let Some(n) = n.and_then(|n| n.checked_mul(5)).filter(|n| n.abs() < EXACT_INTEGER) {
            return n as f64 / (9 * scale) as f64;
        }
    } else if f.fract() == 0.0 && f.abs() < 1e12 {
        // A whole number written with no decimals (short_decimal's exponent is ≥ 0).
        return fahrenheit_to_kelvin_whole(f);
    }
    (f + 459.67) * 5.0 / 9.0
}

fn fahrenheit_to_kelvin_whole(f: f64) -> f64 {
    let n = (f as i128) * 100 + 45967;
    (n * 5) as f64 / 900.0
}

/// The inverse of [`fahrenheit_to_kelvin`]: a kelvin value that one gave from a decimal of up to
/// six places comes back as that decimal; any other in binary.
pub(super) fn kelvin_to_fahrenheit(k: f64) -> f64 {
    for p in 2..=6u32 {
        let scale = 10f64.powi(p as i32);
        let q = k * 9.0 * scale / 5.0;
        let n = q.round();
        if (q - n).abs() < 1e-6 && n.abs() < EXACT_INTEGER as f64 {
            let digits = n as i128 - 45967 * 10i128.pow(p) / 100;
            let f = from_decimal(digits, -(p as i32));
            if fahrenheit_to_kelvin(f) == k {
                return f;
            }
        }
    }
    k * 9.0 / 5.0 - 459.67
}

/// `a + b`, exact in decimal when both are short decimals (an offset such as 273.15 added to a
/// typed temperature), else in binary.
pub(super) fn add(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        return a + b;
    }
    let (Some((ma, ea)), Some((mb, eb))) = (short_decimal(a), short_decimal(b)) else { return a + b };
    let e = ea.min(eb);
    let scale = |m: i128, from: i32| -> Option<i128> {
        let n = u32::try_from(from - e).ok()?;
        m.checked_mul(10i128.checked_pow(n)?)
    };
    match (scale(ma, ea), scale(mb, eb)) {
        (Some(x), Some(y)) => match x.checked_add(y) {
            Some(0) => 0.0,
            Some(sum) => from_decimal(sum, e),
            None => a + b,
        },
        _ => a + b,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_forms() {
        assert_eq!(short_decimal(3.97), Some((397, -2)));
        assert_eq!(short_decimal(-2.1e11), Some((-21, 10)));
        assert_eq!(short_decimal(7.85e-9), Some((785, -11)));
        assert_eq!(short_decimal(0.1 + 0.2), None);
        assert_eq!(short_decimal(0.0), None);
        assert_eq!(short_decimal(f64::NAN), None);
    }

    #[test]
    fn typed_values_round_trip_exactly() {
        // Every value of the form i/100 and i/7 rounded to 12 digits, through ×10^-3 and back.
        for i in 1..50_000 {
            for v in [f64::from(i) / 100.0, format!("{:.12e}", f64::from(i) / 7.0).parse::<f64>().unwrap()] {
                for k in [-3, -6, 6, 9, 12, -12] {
                    let there = shift(v, k);
                    assert_eq!(shift(there, -k), v, "{v} by 1e{k}");
                }
            }
        }
        assert_eq!(shift(3.97, -3), 0.00397);
        assert_eq!(shift(7.85e-9, 12), 7850.0);
        assert_eq!(shift(210000.0, 6), 2.1e11);
    }

    #[test]
    fn long_values_scale_in_binary() {
        let v = 0.1 + 0.2;
        assert_eq!(shift(v, -3), v / 1e3);
        assert_eq!(shift(v, 3), v * 1e3);
        assert_eq!(shift(v, 30), v * 1e22 * 1e8);
        assert_eq!(shift(-0.0, 3).to_bits(), (-0.0f64).to_bits());
        assert!(shift(1e300, 30).is_infinite());
    }

    #[test]
    fn fahrenheit_converts_on_its_decimal() {
        assert_eq!(fahrenheit_to_kelvin(68.0), 293.15);
        assert_eq!(kelvin_to_fahrenheit(353.15), 176.0);
        assert_eq!(kelvin_to_fahrenheit(fahrenheit_to_kelvin(-196.0)), -196.0);
        assert_eq!(kelvin_to_fahrenheit(fahrenheit_to_kelvin(98.6)), 98.6);
    }

    #[test]
    fn offsets_add_in_decimal() {
        assert_eq!(add(20.0, 273.15), 293.15);
        assert_eq!(add(-200.0, 273.15), 73.15);
        assert_eq!(add(293.15, -273.15), 20.0);
        assert_eq!(add(73.15, -273.15), -200.0);
        assert_eq!(add(273.15, -273.15), 0.0);
        for i in -27315..100_000 {
            let c = f64::from(i) / 100.0;
            assert_eq!(add(add(c, 273.15), -273.15), c, "{c}");
        }
    }
}
