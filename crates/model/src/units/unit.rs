//! Named units and unit expressions (`MPa`, `t/mm³`, `W/(m·K)`), each reduced to a dimension and
//! a factor to SI.
//!
//! This table is the only place in FieldCraft where a conversion factor is written down.

use super::{Dimension, UnitError};

/// Exact definitions of the non-metric units (NIST SP 811).
const INCH: f64 = 0.0254;
const FOOT: f64 = 0.3048;
const POUND_MASS: f64 = 0.453_592_37;
const STANDARD_GRAVITY: f64 = 9.806_65;
const POUND_FORCE: f64 = POUND_MASS * STANDARD_GRAVITY;
/// Celsius zero in kelvin.
/// 0 °C in kelvin.
pub const CELSIUS_ZERO: f64 = 273.15;
/// One degree Fahrenheit (or Rankine) in kelvin.
const FAHRENHEIT_STEP: f64 = 5.0 / 9.0;
/// Fahrenheit zero in kelvin: 459.67 °R.
const FAHRENHEIT_ZERO: f64 = 459.67 * FAHRENHEIT_STEP;

const L: Dimension = Dimension::LENGTH;
const M: Dimension = Dimension::MASS;
const T: Dimension = Dimension::TIME;
const TH: Dimension = Dimension::TEMPERATURE;
const ONE: Dimension = Dimension::NONE;
const FORCE: Dimension = Dimension::new(1, 1, -2, 0);
const PRESSURE: Dimension = Dimension::new(-1, 1, -2, 0);
const ENERGY: Dimension = Dimension::new(2, 1, -2, 0);
const POWER: Dimension = Dimension::new(2, 1, -3, 0);
const FREQUENCY: Dimension = Dimension::new(0, 0, -1, 0);

/// A named unit: its spellings (the first is the canonical symbol), dimension, the size of one
/// unit in SI as `mantissa × 10^exp10` (prefixes are exact powers of ten, so `t/mm³` is exactly
/// 10^12 and its conversions are exact in decimal, see [`super::decimal`]), the SI value of its
/// zero (non-zero only for °C and °F), and whether it is an angle.
struct Atom {
    names: &'static [&'static str],
    dimension: Dimension,
    mantissa: f64,
    exp10: i32,
    offset: f64,
    angle: bool,
}

/// A unit that is a power of ten of the SI one (`mm` is 10^-3 m).
const fn ten(names: &'static [&'static str], dimension: Dimension, exp10: i32) -> Atom {
    Atom { names, dimension, mantissa: 1.0, exp10, offset: 0.0, angle: false }
}

/// A unit with another factor (`in` is 0.0254 m).
const fn times(names: &'static [&'static str], dimension: Dimension, mantissa: f64) -> Atom {
    Atom { names, dimension, mantissa, exp10: 0, offset: 0.0, angle: false }
}

const fn angle(names: &'static [&'static str], dimension: Dimension, mantissa: f64) -> Atom {
    Atom { names, dimension, mantissa, exp10: 0, offset: 0.0, angle: true }
}

const ATOMS: &[Atom] = &[
    // Length.
    ten(&["m"], L, 0),
    ten(&["km"], L, 3),
    ten(&["cm"], L, -2),
    ten(&["mm"], L, -3),
    ten(&["µm", "μm", "um"], L, -6),
    ten(&["nm"], L, -9),
    times(&["in", "inch"], L, INCH),
    times(&["ft"], L, FOOT),
    // Mass.
    ten(&["kg"], M, 0),
    ten(&["g"], M, -3),
    ten(&["mg"], M, -6),
    ten(&["t", "tonne"], M, 3),
    times(&["lb", "lbm"], M, POUND_MASS),
    // Time.
    ten(&["s"], T, 0),
    ten(&["ms"], T, -3),
    ten(&["µs", "μs", "us"], T, -6),
    times(&["min"], T, 60.0),
    times(&["h"], T, 3600.0),
    ten(&["Hz"], FREQUENCY, 0),
    ten(&["kHz"], FREQUENCY, 3),
    // Temperature: as an absolute temperature °C and °F carry their offset; inside a compound unit
    // (W/(m·°C)) or as a temperature difference they are only a step size.
    ten(&["K"], TH, 0),
    Atom { names: &["°C", "degC", "℃"], dimension: TH, mantissa: 1.0, exp10: 0, offset: CELSIUS_ZERO, angle: false },
    Atom { names: &["°F", "degF", "℉"], dimension: TH, mantissa: FAHRENHEIT_STEP, exp10: 0, offset: FAHRENHEIT_ZERO, angle: false },
    // Force.
    ten(&["N"], FORCE, 0),
    ten(&["mN"], FORCE, -3),
    ten(&["kN"], FORCE, 3),
    ten(&["MN"], FORCE, 6),
    times(&["lbf"], FORCE, POUND_FORCE),
    Atom { names: &["kip"], dimension: FORCE, mantissa: POUND_FORCE, exp10: 3, offset: 0.0, angle: false },
    // Pressure and stress.
    ten(&["Pa"], PRESSURE, 0),
    ten(&["kPa"], PRESSURE, 3),
    ten(&["MPa"], PRESSURE, 6),
    ten(&["GPa"], PRESSURE, 9),
    ten(&["bar"], PRESSURE, 5),
    times(&["psi"], PRESSURE, POUND_FORCE / (INCH * INCH)),
    Atom { names: &["ksi"], dimension: PRESSURE, mantissa: POUND_FORCE / (INCH * INCH), exp10: 3, offset: 0.0, angle: false },
    // Energy and power.
    ten(&["J"], ENERGY, 0),
    ten(&["mJ"], ENERGY, -3),
    ten(&["kJ"], ENERGY, 3),
    ten(&["W"], POWER, 0),
    ten(&["mW"], POWER, -3),
    ten(&["kW"], POWER, 3),
    // Dimensionless: angles and ratios.
    angle(&["rad"], ONE, 1.0),
    angle(&["deg", "°"], ONE, std::f64::consts::PI / 180.0),
    // Revolutions per minute: an angular velocity.
    angle(&["rpm"], FREQUENCY, std::f64::consts::TAU / 60.0),
    ten(&["%"], ONE, -2),
];

/// Longest accepted unit or quantity text; longer input is refused unparsed.
pub(super) const MAX_TEXT: usize = 64;
/// Largest accepted exponent in a unit expression (`m^6` is the largest anything needs).
const MAX_POWER: i8 = 6;
/// Deepest accepted parenthesis nesting.
const MAX_DEPTH: usize = 4;
/// Largest power of ten a unit's factor may have (`(km)^6` is 10^18).
const MAX_EXP10: u32 = 60;

/// A unit, possibly compound, reduced to SI: `value_si = value × mantissa × 10^exp10 + offset`.
#[derive(Debug, Clone, PartialEq)]
pub struct Unit {
    symbol: String,
    dimension: Dimension,
    /// The factor to SI is `mantissa × 10^exp10`: the power of ten of the prefixes is kept apart,
    /// so it scales exactly in decimal ([`super::decimal`]).
    mantissa: f64,
    exp10: i32,
    /// The SI value of the unit's zero: non-zero only for a lone °C or °F.
    offset: f64,
    /// Whether the unit names an angle (`rad`, `deg`, `°`): only angles accept it.
    angle: bool,
}

impl Unit {
    /// Parses a unit expression: named units joined by `*`, `·`, a space or `/`, with integer
    /// powers (`^2`, `²`, `^-1`, `⁻¹`, or trailing digits as in `mm2`) and parentheses; a leading
    /// `1` is a numerator (`1/K`). `1` alone is the unit of a dimensionless quantity. `.` is not a
    /// multiplication sign (it reads as a decimal point).
    pub fn parse(text: &str) -> Result<Unit, UnitError> {
        let text = text.trim();
        if text.is_empty() {
            return Err(UnitError::Empty);
        }
        if text.len() > MAX_TEXT {
            return Err(UnitError::TooLong);
        }
        if text == "1" {
            return Ok(Unit { symbol: text.to_string(), dimension: Dimension::NONE, mantissa: 1.0, exp10: 0, offset: 0.0, angle: false });
        }
        let mut p = Parser { chars: text.chars().collect(), pos: 0, text, angle: false };
        let Parsed { dimension, mantissa, exp10, atoms, offset } = p.expression(0)?;
        if p.pos != p.chars.len() {
            return Err(p.syntax());
        }
        // An offset only means something for a lone temperature unit.
        let offset = if atoms == 1 { offset } else { 0.0 };
        let factor = super::decimal::shift(mantissa, exp10);
        if !(mantissa.is_finite() && mantissa > 0.0 && factor.is_finite() && factor > 0.0 && exp10.unsigned_abs() <= MAX_EXP10) {
            return Err(UnitError::OutOfRange(text.to_string()));
        }
        Ok(Unit { symbol: text.to_string(), dimension, mantissa, exp10, offset, angle: p.angle })
    }

    /// The unit of a ratio, `1`.
    pub fn one() -> Unit {
        Unit { symbol: "1".to_string(), dimension: Dimension::NONE, mantissa: 1.0, exp10: 0, offset: 0.0, angle: false }
    }

    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    pub fn dimension(&self) -> Dimension {
        self.dimension
    }

    /// The size of one unit in SI (for °C and °F, the size of one degree).
    pub fn factor(&self) -> f64 {
        super::decimal::shift(self.mantissa, self.exp10)
    }

    /// Whether this unit is SI's for its dimension (factor exactly one, no offset): converting
    /// through it changes nothing.
    pub fn is_si(&self) -> bool {
        self.mantissa == 1.0 && self.exp10 == 0 && self.offset == 0.0
    }

    /// Whether the unit names an angle (`rad`, `deg`, `°`).
    pub fn is_angle(&self) -> bool {
        self.angle
    }

    /// Converts a value in this unit to SI. `absolute` applies the zero offset (absolute
    /// temperatures); temperature differences and every other quantity pass `false`.
    pub(super) fn to_si(&self, value: f64, absolute: bool) -> f64 {
        // °F on its decimal: (F + 459.67) × 5/9 as one rounding of exact integers, and back to
        // the decimal it came from, so whole and typed °C and °F convert both ways exactly (68 °F
        // is 293.15 K, 80 °C is 176 °F).
        if self.is_fahrenheit() {
            return if absolute { super::decimal::fahrenheit_to_kelvin(value) } else { value * 5.0 / 9.0 };
        }
        let scaled = if self.mantissa == 1.0 { value } else { value * self.mantissa };
        let scaled = super::decimal::shift(scaled, self.exp10);
        if absolute && self.offset != 0.0 { super::decimal::add(scaled, self.offset) } else { scaled }
    }

    /// Whether it is the degree Fahrenheit alone (whose step 5/9 K is no decimal).
    fn is_fahrenheit(&self) -> bool {
        self.mantissa == FAHRENHEIT_STEP && self.exp10 == 0 && self.offset == FAHRENHEIT_ZERO
    }

    /// Converts an SI value to this unit; the inverse of [`Unit::to_si`].
    pub(super) fn express(&self, si: f64, absolute: bool) -> f64 {
        if self.is_fahrenheit() {
            return if absolute { super::decimal::kelvin_to_fahrenheit(si) } else { si * 9.0 / 5.0 };
        }
        let v = if absolute && self.offset != 0.0 { super::decimal::add(si, -self.offset) } else { si };
        let v = super::decimal::shift(v, -self.exp10);
        if self.mantissa == 1.0 { v } else { v / self.mantissa }
    }
}

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    text: &'a str,
    angle: bool,
}

/// A parsed unit expression: its dimension, factor (`mantissa × 10^exp10`), number of named
/// units and the offset of the last named unit.
#[derive(Clone, Copy)]
struct Parsed {
    dimension: Dimension,
    mantissa: f64,
    exp10: i32,
    atoms: usize,
    offset: f64,
}

impl Parser<'_> {
    fn syntax(&self) -> UnitError {
        UnitError::Syntax(self.text.to_string())
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn skip_spaces(&mut self) -> bool {
        let start = self.pos;
        while self.peek().is_some_and(char::is_whitespace) {
            self.pos += 1;
        }
        self.pos > start
    }

    fn expression(&mut self, depth: usize) -> Result<Parsed, UnitError> {
        self.skip_spaces();
        // A leading `1` is a numerator: `1/K`, `1/s`.
        let leading_one = self.peek() == Some('1') && !self.chars.get(self.pos + 1).is_some_and(char::is_ascii_digit);
        let mut acc = if leading_one {
            self.pos += 1;
            self.skip_spaces();
            if self.peek() != Some('/') {
                return Err(self.syntax());
            }
            Parsed { dimension: Dimension::NONE, mantissa: 1.0, exp10: 0, atoms: 0, offset: 0.0 }
        } else {
            self.term(depth)?
        };
        loop {
            let spaced = self.skip_spaces();
            let divide = match self.peek() {
                Some('*' | '·') => false,
                Some('/') => true,
                // A space between two units multiplies: `N mm`, `W/(m K)`.
                Some(c) if spaced && (c == '(' || c.is_alphabetic() || matches!(c, '°' | '℃' | '℉' | '%' | 'µ' | 'μ')) => {
                    self.pos -= 1;
                    false
                }
                _ => break,
            };
            self.pos += 1;
            self.skip_spaces();
            let t = self.term(depth)?;
            let overflow = || UnitError::OutOfRange(self.text.to_string());
            if divide {
                acc.dimension = acc.dimension.checked_div(t.dimension).ok_or_else(overflow)?;
                acc.mantissa /= t.mantissa;
                acc.exp10 = acc.exp10.checked_sub(t.exp10).ok_or_else(overflow)?;
            } else {
                acc.dimension = acc.dimension.checked_mul(t.dimension).ok_or_else(overflow)?;
                acc.mantissa *= t.mantissa;
                acc.exp10 = acc.exp10.checked_add(t.exp10).ok_or_else(overflow)?;
            }
            acc.atoms += t.atoms;
            acc.offset = t.offset;
        }
        Ok(acc)
    }

    fn term(&mut self, depth: usize) -> Result<Parsed, UnitError> {
        let base = if self.peek() == Some('(') {
            if depth >= MAX_DEPTH {
                return Err(self.syntax());
            }
            self.pos += 1;
            let inner = self.expression(depth + 1)?;
            self.skip_spaces();
            if self.peek() != Some(')') {
                return Err(self.syntax());
            }
            self.pos += 1;
            inner
        } else {
            self.atom()?
        };
        let power = self.power()?;
        if power == 1 {
            return Ok(base);
        }
        let overflow = || UnitError::OutOfRange(self.text.to_string());
        let dimension = base.dimension.checked_pow(power).ok_or_else(overflow)?;
        let exp10 = base.exp10.checked_mul(i32::from(power)).ok_or_else(overflow)?;
        // A powered unit is never a lone absolute temperature.
        Ok(Parsed { dimension, mantissa: base.mantissa.powi(i32::from(power)), exp10, atoms: base.atoms + 1, offset: 0.0 })
    }

    fn atom(&mut self) -> Result<Parsed, UnitError> {
        let start = self.pos;
        while let Some(c) = self.peek() {
            if c.is_alphabetic() || matches!(c, '°' | '℃' | '℉' | '%' | 'µ' | 'μ') {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return Err(self.syntax());
        }
        let name: String = self.chars.get(start..self.pos).unwrap_or_default().iter().collect();
        let atom = ATOMS.iter().find(|a| a.names.contains(&name.as_str())).ok_or_else(|| UnitError::UnknownUnit(name.clone()))?;
        self.angle |= atom.angle;
        Ok(Parsed { dimension: atom.dimension, mantissa: atom.mantissa, exp10: atom.exp10, atoms: 1, offset: atom.offset })
    }

    /// An optional integer exponent; 1 when absent.
    fn power(&mut self) -> Result<i8, UnitError> {
        let mut digits = String::new();
        if self.peek() == Some('^') {
            self.pos += 1;
            if self.peek() == Some('-') {
                digits.push('-');
                self.pos += 1;
            }
            while let Some(c) = self.peek().filter(char::is_ascii_digit) {
                digits.push(c);
                self.pos += 1;
            }
        } else {
            while let Some(c) = self.peek() {
                let mapped = match c {
                    '⁻' => '-',
                    '¹' => '1',
                    '²' => '2',
                    '³' => '3',
                    '⁴' => '4',
                    '⁵' => '5',
                    '⁶' => '6',
                    d if d.is_ascii_digit() => d,
                    _ => break,
                };
                digits.push(mapped);
                self.pos += 1;
            }
            if digits.is_empty() {
                return Ok(1);
            }
        }
        match digits.parse::<i8>() {
            Ok(p) if p != 0 && p.unsigned_abs() <= MAX_POWER.unsigned_abs() => Ok(p),
            _ => Err(self.syntax()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
    }

    #[test]
    fn compound_units_reduce_to_si() {
        let cases = [
            ("t/mm³", Dimension::new(-3, 1, 0, 0), 1e12),
            ("t/mm^3", Dimension::new(-3, 1, 0, 0), 1e12),
            ("kg/m3", Dimension::new(-3, 1, 0, 0), 1.0),
            ("N/mm²", PRESSURE, 1e6),
            ("W/(m·K)", Dimension::new(1, 1, -3, -1), 1.0),
            ("mW/(mm*K)", Dimension::new(1, 1, -3, -1), 1.0),
            ("W/m/K", Dimension::new(1, 1, -3, -1), 1.0),
            ("K^-1", Dimension::new(0, 0, 0, -1), 1.0),
            ("lbf·s²/in", M, POUND_FORCE / INCH),
            ("mm/s", Dimension::new(1, 0, -1, 0), 1e-3),
            ("MPa·s", Dimension::new(-1, 1, -1, 0), 1e6),
        ];
        for (text, dim, factor) in cases {
            let u = Unit::parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
            assert_eq!(u.dimension(), dim, "{text}");
            assert!(close(u.factor(), factor), "{text}: {} vs {factor}", u.factor());
        }
    }

    #[test]
    fn numerators_and_spaces() {
        let k = Unit::parse("1/K").unwrap();
        assert_eq!(k.dimension(), Dimension::new(0, 0, 0, -1));
        assert_eq!(Unit::parse("1/s").unwrap().dimension(), FREQUENCY);
        assert_eq!(Unit::parse("1 / (m·K)").unwrap().dimension(), Dimension::new(-1, 0, 0, -1));
        for (spaced, plain) in
            [("N · mm", "N·mm"), ("W/(m K)", "W/(m·K)"), ("kg m^-3", "kg/m^3"), ("W / ( m * K )", "W/(m·K)"), ("lbf s² / in", "lbf·s²/in")]
        {
            let (a, b) = (Unit::parse(spaced).unwrap(), Unit::parse(plain).unwrap());
            assert_eq!(a.dimension(), b.dimension(), "{spaced}");
            assert!(close(a.factor(), b.factor()), "{spaced}");
        }
        assert!(Unit::parse("deg").unwrap().is_angle());
        assert!(Unit::parse("°").unwrap().is_angle());
        assert!(!Unit::parse("%").unwrap().is_angle());
    }

    #[test]
    fn offsets_only_for_lone_temperatures() {
        assert!(close(Unit::parse("°C").unwrap().to_si(20.0, true), 293.15));
        assert!(close(Unit::parse("°F").unwrap().to_si(32.0, true), 273.15));
        assert!(close(Unit::parse("°F").unwrap().to_si(212.0, true), 373.15));
        // In a compound unit the degree is a step.
        let k = Unit::parse("W/(m·°C)").unwrap();
        assert!(close(k.to_si(1.0, true), 1.0));
    }

    #[test]
    fn hostile_units_are_refused() {
        for text in [
            "",
            "   ",
            "m^",
            "m^0",
            "m^99",
            "m^-128",
            "((((((m))))))",
            "(m",
            "m)",
            "m//s",
            "mm2x",
            "foo",
            "m^7",
            "BTU",
            "m.s",
            "1 K",
            "11/K",
            "2/K",
            "1/",
        ] {
            assert!(Unit::parse(text).is_err(), "{text:?}");
        }
        assert!(matches!(Unit::parse(&"m".repeat(100)), Err(UnitError::TooLong)));
        // Exponent overflow across many factors.
        assert!(Unit::parse("(m^6)^6·(m^6)^6·(m^6)^6·(m^6)^6").is_err());
    }
}
