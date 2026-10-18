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
const CELSIUS_ZERO: f64 = 273.15;
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
/// unit in SI, and the SI value of its zero (non-zero only for °C and °F).
struct Atom {
    names: &'static [&'static str],
    dimension: Dimension,
    factor: f64,
    offset: f64,
}

const fn atom(names: &'static [&'static str], dimension: Dimension, factor: f64) -> Atom {
    Atom { names, dimension, factor, offset: 0.0 }
}

const ATOMS: &[Atom] = &[
    // Length.
    atom(&["m"], L, 1.0),
    atom(&["km"], L, 1e3),
    atom(&["cm"], L, 1e-2),
    atom(&["mm"], L, 1e-3),
    atom(&["µm", "μm", "um"], L, 1e-6),
    atom(&["nm"], L, 1e-9),
    atom(&["in", "inch"], L, INCH),
    atom(&["ft"], L, FOOT),
    // Mass.
    atom(&["kg"], M, 1.0),
    atom(&["g"], M, 1e-3),
    atom(&["mg"], M, 1e-6),
    atom(&["t", "tonne"], M, 1e3),
    atom(&["lb", "lbm"], M, POUND_MASS),
    // Time.
    atom(&["s"], T, 1.0),
    atom(&["ms"], T, 1e-3),
    atom(&["µs", "μs", "us"], T, 1e-6),
    atom(&["min"], T, 60.0),
    atom(&["h"], T, 3600.0),
    atom(&["Hz"], FREQUENCY, 1.0),
    atom(&["kHz"], FREQUENCY, 1e3),
    // Temperature: as an absolute temperature °C and °F carry their offset; inside a compound unit
    // (W/(m·°C)) or as a temperature difference they are only a step size.
    atom(&["K"], TH, 1.0),
    Atom { names: &["°C", "degC", "℃"], dimension: TH, factor: 1.0, offset: CELSIUS_ZERO },
    Atom { names: &["°F", "degF", "℉"], dimension: TH, factor: FAHRENHEIT_STEP, offset: FAHRENHEIT_ZERO },
    // Force.
    atom(&["N"], FORCE, 1.0),
    atom(&["mN"], FORCE, 1e-3),
    atom(&["kN"], FORCE, 1e3),
    atom(&["MN"], FORCE, 1e6),
    atom(&["lbf"], FORCE, POUND_FORCE),
    atom(&["kip"], FORCE, 1e3 * POUND_FORCE),
    // Pressure and stress.
    atom(&["Pa"], PRESSURE, 1.0),
    atom(&["kPa"], PRESSURE, 1e3),
    atom(&["MPa"], PRESSURE, 1e6),
    atom(&["GPa"], PRESSURE, 1e9),
    atom(&["bar"], PRESSURE, 1e5),
    atom(&["psi"], PRESSURE, POUND_FORCE / (INCH * INCH)),
    atom(&["ksi"], PRESSURE, 1e3 * POUND_FORCE / (INCH * INCH)),
    // Energy and power.
    atom(&["J"], ENERGY, 1.0),
    atom(&["mJ"], ENERGY, 1e-3),
    atom(&["kJ"], ENERGY, 1e3),
    atom(&["W"], POWER, 1.0),
    atom(&["mW"], POWER, 1e-3),
    atom(&["kW"], POWER, 1e3),
    // Dimensionless: angles and ratios.
    atom(&["rad"], ONE, 1.0),
    atom(&["deg", "°"], ONE, std::f64::consts::PI / 180.0),
    atom(&["%"], ONE, 1e-2),
];

/// Longest accepted unit or quantity text; longer input is refused unparsed.
pub(super) const MAX_TEXT: usize = 64;
/// Largest accepted exponent in a unit expression (`m^6` is the largest anything needs).
const MAX_POWER: i8 = 6;
/// Deepest accepted parenthesis nesting.
const MAX_DEPTH: usize = 4;

/// A unit, possibly compound, reduced to SI: `value_si = value * factor + offset`.
#[derive(Debug, Clone, PartialEq)]
pub struct Unit {
    symbol: String,
    dimension: Dimension,
    factor: f64,
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
            return Ok(Unit { symbol: text.to_string(), dimension: Dimension::NONE, factor: 1.0, offset: 0.0, angle: false });
        }
        let mut p = Parser { chars: text.chars().collect(), pos: 0, text, angle: false };
        let (dimension, factor, atoms, offset) = p.expression(0)?;
        if p.pos != p.chars.len() {
            return Err(p.syntax());
        }
        // An offset only means something for a lone temperature unit.
        let offset = if atoms == 1 { offset } else { 0.0 };
        if !(factor.is_finite() && factor > 0.0) {
            return Err(UnitError::OutOfRange(text.to_string()));
        }
        Ok(Unit { symbol: text.to_string(), dimension, factor, offset, angle: p.angle })
    }

    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    pub fn dimension(&self) -> Dimension {
        self.dimension
    }

    /// The size of one unit in SI (for °C and °F, the size of one degree).
    pub fn factor(&self) -> f64 {
        self.factor
    }

    /// Whether the unit names an angle (`rad`, `deg`, `°`).
    pub fn is_angle(&self) -> bool {
        self.angle
    }

    /// Converts a value in this unit to SI. `absolute` applies the zero offset (absolute
    /// temperatures); temperature differences and every other quantity pass `false`.
    pub(super) fn to_si(&self, value: f64, absolute: bool) -> f64 {
        value * self.factor + if absolute { self.offset } else { 0.0 }
    }

    /// Converts an SI value to this unit; the inverse of [`Unit::to_si`].
    pub(super) fn express(&self, si: f64, absolute: bool) -> f64 {
        (si - if absolute { self.offset } else { 0.0 }) / self.factor
    }
}

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    text: &'a str,
    angle: bool,
}

/// Dimension, factor, number of named units, offset of the last named unit.
type Parsed = (Dimension, f64, usize, f64);

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
        let (mut dim, mut factor, mut atoms, mut offset) = if leading_one {
            self.pos += 1;
            self.skip_spaces();
            if self.peek() != Some('/') {
                return Err(self.syntax());
            }
            (Dimension::NONE, 1.0, 0, 0.0)
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
            let (d, f, a, o) = self.term(depth)?;
            let overflow = || UnitError::OutOfRange(self.text.to_string());
            if divide {
                dim = dim.checked_div(d).ok_or_else(overflow)?;
                factor /= f;
            } else {
                dim = dim.checked_mul(d).ok_or_else(overflow)?;
                factor *= f;
            }
            atoms += a;
            offset = o;
        }
        Ok((dim, factor, atoms, offset))
    }

    fn term(&mut self, depth: usize) -> Result<Parsed, UnitError> {
        let (dim, factor, atoms, offset) = if self.peek() == Some('(') {
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
            return Ok((dim, factor, atoms, offset));
        }
        let dim = dim.checked_pow(power).ok_or_else(|| UnitError::OutOfRange(self.text.to_string()))?;
        // A powered unit is never a lone absolute temperature.
        Ok((dim, factor.powi(i32::from(power)), atoms + 1, 0.0))
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
        self.angle |= matches!(atom.names.first(), Some(&("rad" | "deg")));
        Ok((atom.dimension, atom.factor, 1, atom.offset))
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
