//! Reading an OpenFOAM field file's internal field (ASCII format), hostile-input safe.
//!
//! The format comes from OpenFOAM's public user guide ("Basic input/output file format"): a
//! `FoamFile { … }` header dictionary, then entries ending in `;`, with `//` and `/* … */`
//! comments. Only the `internalField` entry is read; its values are counted and checked against
//! the bound the caller passes, before any value is read. No allocation is sized from a count
//! read from the file, and the text is scanned in place, never copied.

use crate::{Error, Result};

const FORMAT: &str = "OpenFOAM";

/// How many values the caller accepts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Count {
    /// Exactly this many (a field on a mesh whose cell count is known).
    Exactly(usize),
    /// At most this many (the cell centres, which set the count).
    AtMost(usize),
}

/// The values of an internal field.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldValues {
    UniformScalar(f64),
    UniformVector([f64; 3]),
    Scalars(Vec<f64>),
    Vectors(Vec<[f64; 3]>),
}

impl FieldValues {
    /// The scalar at cell `i` (a uniform field has one value for every cell).
    pub fn scalar(&self, i: usize) -> Option<f64> {
        match self {
            Self::UniformScalar(value) => Some(*value),
            Self::Scalars(values) => values.get(i).copied(),
            Self::UniformVector(_) | Self::Vectors(_) => None,
        }
    }

    /// The vector at cell `i` (a uniform field has one value for every cell).
    pub fn vector(&self, i: usize) -> Option<[f64; 3]> {
        match self {
            Self::UniformVector(value) => Some(*value),
            Self::Vectors(values) => values.get(i).copied(),
            Self::UniformScalar(_) | Self::Scalars(_) => None,
        }
    }

    /// The number of values in a nonuniform field; `None` for a uniform one.
    pub fn count(&self) -> Option<usize> {
        match self {
            Self::UniformScalar(_) | Self::UniformVector(_) => None,
            Self::Scalars(values) => Some(values.len()),
            Self::Vectors(values) => Some(values.len()),
        }
    }
}

/// Components of a vector value.
const VECTOR_LEN: usize = 3;
/// How many characters of file text an error message quotes.
const QUOTE_CHARS: usize = 32;
/// The `List<…>` type names this reader takes.
const SCALAR_LIST: &str = "List<scalar>";
const VECTOR_LIST: &str = "List<vector>";
const BINARY_REFUSAL: &str = "a binary field file (the case is written with writeFormat ascii)";

/// The internal field of an ASCII field file. The text must tokenize to its end (a comment or
/// string left open is an error); only the `internalField` entry is read into values.
pub fn read_field(text: &str, count: Count) -> Result<FieldValues> {
    let mut lexer = Lexer { text, pos: 0 };
    let mut depth = 0_usize;
    while let Some(token) = lexer.next_token()? {
        match token {
            Token::BraceOpen => depth += 1,
            Token::BraceClose => {
                depth = depth.checked_sub(1).ok_or_else(|| read_error("a closing } with no opening {"))?;
            }
            Token::Word("FoamFile") if depth == 0 => read_header(&mut lexer)?,
            Token::Word("internalField") if depth == 0 => {
                let values = read_internal_field(&mut lexer, count)?;
                while lexer.next_token()?.is_some() {}
                return Ok(values);
            }
            Token::Word(_) | Token::Open | Token::Close | Token::Semi => {}
        }
    }
    Err(read_error("no internalField"))
}

/// One lexical token of a field file. Quoted strings are words (their quotes included).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Token<'a> {
    Open,
    Close,
    BraceOpen,
    BraceClose,
    Semi,
    Word(&'a str),
}

impl<'a> Token<'a> {
    fn text(self) -> &'a str {
        match self {
            Self::Open => "(",
            Self::Close => ")",
            Self::BraceOpen => "{",
            Self::BraceClose => "}",
            Self::Semi => ";",
            Self::Word(word) => word,
        }
    }
}

/// Splits the text into tokens, borrowing from it (no copy).
struct Lexer<'a> {
    text: &'a str,
    /// The byte offset of the next unread character; always a character boundary.
    pos: usize,
}

impl<'a> Lexer<'a> {
    fn next_token(&mut self) -> Result<Option<Token<'a>>> {
        loop {
            let rest = self.text.get(self.pos..).unwrap_or_default();
            let Some(first) = rest.chars().next() else {
                return Ok(None);
            };
            if first.is_whitespace() {
                self.pos += first.len_utf8();
            } else if rest.starts_with("//") {
                self.pos += rest.find('\n').unwrap_or(rest.len());
            } else if rest.starts_with("/*") {
                let body = rest.get(2..).unwrap_or_default();
                let close = body.find("*/").ok_or_else(|| read_error("an unterminated /* comment"))?;
                self.pos += 2 + close + 2;
            } else {
                return self.token(rest, first).map(Some);
            }
        }
    }

    /// The token that starts at `rest`, whose first character is `first`.
    fn token(&mut self, rest: &'a str, first: char) -> Result<Token<'a>> {
        let (token, len) = match first {
            '(' => (Token::Open, 1),
            ')' => (Token::Close, 1),
            '{' => (Token::BraceOpen, 1),
            '}' => (Token::BraceClose, 1),
            ';' => (Token::Semi, 1),
            '"' => {
                let body = rest.get(1..).unwrap_or_default();
                let close = body.find('"').ok_or_else(|| read_error("an unterminated string"))?;
                let len = 2 + close;
                (Token::Word(rest.get(..len).unwrap_or_default()), len)
            }
            _ => {
                let len = word_len(rest);
                (Token::Word(rest.get(..len).unwrap_or_default()), len)
            }
        };
        self.pos += len;
        Ok(token)
    }
}

/// The byte length of the bare word at the start of `rest`. Its first character is neither
/// whitespace, a delimiter nor a comment's start (the caller has dispatched those), so it is at
/// least one byte long.
fn word_len(rest: &str) -> usize {
    rest.char_indices()
        .find(|&(i, c)| c.is_whitespace() || is_delimiter(c) || rest.get(i..).is_some_and(|tail| tail.starts_with("//") || tail.starts_with("/*")))
        .map_or(rest.len(), |(i, _)| i)
}

fn is_delimiter(c: char) -> bool {
    matches!(c, '(' | ')' | '{' | '}' | ';' | '"')
}

fn read_error(message: impl Into<String>) -> Error {
    Error::Read { format: FORMAT, message: message.into() }
}

/// The first characters of file text, for a message.
fn quoted(text: &str) -> String {
    text.chars().take(QUOTE_CHARS).collect()
}

/// The next token; the end of the text before it is an error.
fn next<'a>(lexer: &mut Lexer<'a>, expected: &str) -> Result<Token<'a>> {
    lexer.next_token()?.ok_or_else(|| read_error(format!("unexpected end of file, expected {expected}")))
}

fn expect(lexer: &mut Lexer<'_>, wanted: Token<'_>, expected: &str) -> Result<()> {
    let token = next(lexer, expected)?;
    if token == wanted { Ok(()) } else { Err(unexpected(token, expected)) }
}

fn unexpected(token: Token<'_>, expected: &str) -> Error {
    read_error(format!("expected {expected}, found \"{}\"", quoted(token.text())))
}

/// A word that must be a finite decimal number.
fn number(word: &str) -> Result<f64> {
    match word.parse::<f64>() {
        Ok(value) if value.is_finite() => Ok(value),
        _ => Err(read_error(format!("\"{}\" is not a finite number", quoted(word)))),
    }
}

/// Reads the `FoamFile` dictionary (its braces included) and refuses a binary `format`.
fn read_header(lexer: &mut Lexer<'_>) -> Result<()> {
    expect(lexer, Token::BraceOpen, "{ after FoamFile")?;
    let mut depth = 1_usize;
    while depth > 0 {
        match next(lexer, "}")? {
            Token::BraceOpen => depth += 1,
            Token::BraceClose => depth -= 1,
            Token::Word("format") => {
                if next(lexer, "ascii or binary")? == Token::Word("binary") {
                    return Err(Error::Unsupported { format: FORMAT, what: BINARY_REFUSAL.into() });
                }
            }
            Token::Word(_) | Token::Open | Token::Close | Token::Semi => {}
        }
    }
    Ok(())
}

/// Reads the value of `internalField` (after its keyword), up to and including its `;`.
fn read_internal_field(lexer: &mut Lexer<'_>, count: Count) -> Result<FieldValues> {
    let values = match next(lexer, "uniform or nonuniform")? {
        Token::Word("uniform") => read_uniform(lexer)?,
        Token::Word("nonuniform") => read_nonuniform(lexer, count)?,
        other => return Err(unexpected(other, "uniform or nonuniform")),
    };
    expect(lexer, Token::Semi, ";")?;
    Ok(values)
}

fn read_uniform(lexer: &mut Lexer<'_>) -> Result<FieldValues> {
    match next(lexer, "a number or (x y z)")? {
        Token::Word(word) => Ok(FieldValues::UniformScalar(number(word)?)),
        Token::Open => Ok(FieldValues::UniformVector(read_vector_body(lexer)?)),
        other => Err(unexpected(other, "a number or (x y z)")),
    }
}

/// Reads a nonuniform list: its type, its count (checked against `count` before any value),
/// and the values between the parentheses.
fn read_nonuniform(lexer: &mut Lexer<'_>, count: Count) -> Result<FieldValues> {
    let kind = next(lexer, "List<scalar> or List<vector>")?;
    let vector = match kind {
        Token::Word(SCALAR_LIST) => false,
        Token::Word(VECTOR_LIST) => true,
        Token::Word(other) if other.starts_with("List<") => {
            return Err(Error::Unsupported { format: FORMAT, what: format!("a field of type {}", quoted(other)) });
        }
        other => return Err(unexpected(other, "List<scalar> or List<vector>")),
    };
    let n_word = next(lexer, "the number of values")?.text();
    let n = n_word.parse::<usize>().map_err(|_| read_error(format!("the number of values is \"{}\", not a whole number", quoted(n_word))))?;
    check_count(n, count)?;
    expect(lexer, Token::Open, "(")?;
    if vector { Ok(FieldValues::Vectors(read_vectors(lexer, n)?)) } else { Ok(FieldValues::Scalars(read_scalars(lexer, n)?)) }
}

fn check_count(n: usize, count: Count) -> Result<()> {
    match count {
        Count::Exactly(m) if n != m => Err(read_error(format!("a field of {n} values on a mesh of {m} cells"))),
        Count::AtMost(m) if n > m => Err(Error::Limit(format!("{n} values in a field (at most {m})"))),
        Count::Exactly(_) | Count::AtMost(_) => Ok(()),
    }
}

/// Reads scalars up to the closing parenthesis. Values are pushed one by one, and the list
/// stops at the first value past `n`, so the vector never holds more than `n` (which the
/// caller has bounded).
fn read_scalars(lexer: &mut Lexer<'_>, n: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    loop {
        match next(lexer, "a value or )")? {
            Token::Close => break,
            Token::Word(word) => {
                if values.len() == n {
                    return Err(read_error(format!("more than {n} values in the list")));
                }
                values.push(number(word)?);
            }
            other => return Err(unexpected(other, "a value or )")),
        }
    }
    if values.len() != n {
        return Err(read_error(format!("a list of {} values, not the {n} its count gives", values.len())));
    }
    Ok(values)
}

/// Reads vectors up to the closing parenthesis, with the same bound as [`read_scalars`].
fn read_vectors(lexer: &mut Lexer<'_>, n: usize) -> Result<Vec<[f64; VECTOR_LEN]>> {
    let mut values = Vec::new();
    loop {
        match next(lexer, "a vector or )")? {
            Token::Close => break,
            Token::Open => {
                if values.len() == n {
                    return Err(read_error(format!("more than {n} values in the list")));
                }
                values.push(read_vector_body(lexer)?);
            }
            other => return Err(unexpected(other, "a vector or )")),
        }
    }
    if values.len() != n {
        return Err(read_error(format!("a list of {} values, not the {n} its count gives", values.len())));
    }
    Ok(values)
}

/// The components of a vector after its `(`, through its `)`.
fn read_vector_body(lexer: &mut Lexer<'_>) -> Result<[f64; VECTOR_LEN]> {
    let mut vector = [0.0; VECTOR_LEN];
    for component in &mut vector {
        match next(lexer, "a vector component")? {
            Token::Word(word) => *component = number(word)?,
            other => return Err(unexpected(other, "a vector component")),
        }
    }
    expect(lexer, Token::Close, ")")?;
    Ok(vector)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(extra: &str) -> String {
        format!(
            "/*--------------------------------*- C++ -*----------------------------------*\\\n\
             \\*---------------------------------------------------------------------------*/\n\
             FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       volVectorField;\n    object      U;\n}}\n\
             // a comment\n{extra}\n"
        )
    }

    #[test]
    fn reads_uniform_and_nonuniform_fields() {
        let vectors = header(
            "dimensions      [0 1 -1 0 0 0 0];\n\n\
             internalField   nonuniform List<vector> 3\n(\n(1 2 3) (4 5 6)\n(-1e-3 0 7.5)\n);\n\n\
             boundaryField\n{\n    inlet { type fixedValue; value uniform (9 9 9); }\n}\n",
        );
        let field = read_field(&vectors, Count::Exactly(3)).expect("vectors read");
        assert_eq!(field, FieldValues::Vectors(vec![[1., 2., 3.], [4., 5., 6.], [-1e-3, 0., 7.5]]));

        let scalars = header("internalField nonuniform List<scalar> 2(101325 1.5e5);");
        let field = read_field(&scalars, Count::AtMost(10)).expect("scalars read");
        assert_eq!(field, FieldValues::Scalars(vec![101325., 1.5e5]));
        assert_eq!(field.count(), Some(2));

        let uniform_vector = header("internalField uniform (0 0 1);");
        let field = read_field(&uniform_vector, Count::Exactly(7)).expect("uniform vector read");
        assert_eq!(field, FieldValues::UniformVector([0., 0., 1.]));
        assert_eq!(field.vector(6), Some([0., 0., 1.]));
        assert_eq!(field.count(), None);

        let uniform_scalar = header("internalField uniform 293.15;");
        let field = read_field(&uniform_scalar, Count::Exactly(7)).expect("uniform scalar read");
        assert_eq!(field.scalar(6), Some(293.15));
        assert_eq!(field.vector(0), None);
    }

    #[test]
    fn refuses_malformed_fields() {
        let cases: [(&str, Count); 9] = [
            ("internalField nonuniform List<scalar> 3(1 2);", Count::AtMost(10)),
            ("internalField nonuniform List<scalar> 1(1 2);", Count::AtMost(10)),
            ("internalField nonuniform List<vector> 1((1 2));", Count::AtMost(10)),
            ("internalField nonuniform List<tensor> 1((1 2 3 4 5 6 7 8 9));", Count::AtMost(10)),
            ("internalField nonuniform List<scalar> 1(NaN);", Count::AtMost(10)),
            ("internalField uniform inf;", Count::AtMost(10)),
            ("internalField uniform 1; /* unterminated", Count::AtMost(10)),
            ("internalField nonuniform List<scalar> 2(1 2);", Count::Exactly(3)),
            ("dimensions [0 0 0 0 0 0 0];", Count::AtMost(10)),
        ];
        for (internal, count) in cases {
            let text = header(internal);
            assert!(read_field(&text, count).is_err(), "accepted: {internal}");
        }

        let binary = "FoamFile\n{\n    format binary;\n}\ninternalField uniform 1;\n";
        let err = read_field(binary, Count::AtMost(10)).expect_err("binary refused");
        assert!(err.to_string().contains("binary"), "{err}");

        let two = header("internalField nonuniform List<scalar> 3(1 2 3);");
        let err = read_field(&two, Count::Exactly(2)).expect_err("count mismatch");
        assert!(err.to_string().contains("mesh of 2 cells"), "{err}");
    }

    #[test]
    fn a_huge_count_is_refused_before_any_value() {
        let text = "internalField nonuniform List<vector> 50000000(";
        let limited = read_field(text, Count::AtMost(1000)).expect_err("over the bound");
        assert!(matches!(limited, Error::Limit(_)), "{limited}");
        let short = read_field(text, Count::Exactly(50_000_000)).expect_err("list ends early");
        assert!(matches!(short, Error::Read { .. }), "{short}");
    }

    #[test]
    fn deep_or_long_text_never_panics() {
        let deep = format!("{}\ninternalField uniform 1;", "{".repeat(100_000));
        assert!(read_field(&deep, Count::AtMost(10)).is_err());
        let slashes = "/".repeat(1_000_000);
        assert!(read_field(&slashes, Count::AtMost(10)).is_err());
    }
}
