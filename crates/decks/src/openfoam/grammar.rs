//! Our own reader of the OpenFOAM dictionary grammar, used by the tests to check every dictionary
//! the case writer produces. It is test-only: nothing in the product reads OpenFOAM files.
//!
//! The grammar is the one the OpenFOAM user guide describes for its input files ("Basic
//! input/output file format"), restated here in our own words: a file is a sequence of entries;
//! `keyword { entries }` is a sub-dictionary; `keyword value… ;` is a primitive entry; `//` and
//! `/* */` are comments; a word may contain balanced parentheses (`div(phi,U)`).

#[derive(Debug, Clone, PartialEq)]
pub(super) enum Token {
    Word(String),
    Number(f64),
    Str(String),
    Punct(char),
}

#[derive(Debug, Clone, PartialEq)]
pub(super) enum Entry {
    Dict(Dict),
    Value(Vec<Token>),
}

#[derive(Debug, Clone, PartialEq, Default)]
pub(super) struct Dict {
    pub(super) entries: Vec<(String, Entry)>,
    pub(super) directives: Vec<String>,
}

impl Dict {
    /// The entry at a `/`-separated path of keywords (`"castellatedMeshControls/maxLocalCells"`);
    /// the first match at each level.
    pub(super) fn get(&self, path: &str) -> Option<&Entry> {
        let mut dict = self;
        let mut parts = path.split('/').peekable();
        loop {
            let key = parts.next()?;
            let (_, entry) = dict.entries.iter().find(|(k, _)| k == key)?;
            if parts.peek().is_none() {
                return Some(entry);
            }
            match entry {
                Entry::Dict(inner) => dict = inner,
                Entry::Value(_) => return None,
            }
        }
    }

    /// The value at `path` as text: its tokens joined with single spaces, numbers via `{}` of f64,
    /// strings with their quotes, punctuation as itself. `None` if absent or a sub-dictionary.
    pub(super) fn text(&self, path: &str) -> Option<String> {
        match self.get(path)? {
            Entry::Value(tokens) => Some(tokens.iter().map(render).collect::<Vec<_>>().join(" ")),
            Entry::Dict(_) => None,
        }
    }

    /// The value at `path` if it is exactly one number.
    pub(super) fn number(&self, path: &str) -> Option<f64> {
        match self.get(path)? {
            Entry::Value(tokens) => match tokens.as_slice() {
                [Token::Number(n)] => Some(*n),
                _ => None,
            },
            Entry::Dict(_) => None,
        }
    }

    /// Every directive anywhere in the tree (including nested dictionaries), as written.
    pub(super) fn all_directives(&self) -> Vec<String> {
        let mut out = self.directives.clone();
        for (_, entry) in &self.entries {
            if let Entry::Dict(inner) = entry {
                out.extend(inner.all_directives());
            }
        }
        out
    }
}

/// Nesting limit for `{ }` blocks and for brackets inside a value; a hostile file must not
/// overflow the stack.
const MAX_DEPTH: usize = 64;

/// Entry keys a dictionary or field file's `FoamFile` header must carry.
const HEADER_KEYS: [&str; 3] = ["version", "class", "object"];

/// A token with the line it starts on and its byte offset in the source.
#[derive(Debug, Clone)]
struct Lexed {
    tok: Token,
    line: usize,
    start: usize,
}

/// Parse a whole file. Errors are a `String` naming the line (1-based) and what is wrong.
pub(super) fn parse(text: &str) -> Result<Dict, String> {
    let mut parser = Parser { toks: lex(text)?, pos: 0, text };
    parser.entries(0, false)
}

/// Parse a dictionary or field file and require its `FoamFile` header: a sub-dictionary with
/// `version`, `format` (`ascii`), `class` and `object` entries. Returns the whole tree.
pub(super) fn check_file(text: &str) -> Result<Dict, String> {
    let tree = parse(text)?;
    let Some(Entry::Dict(header)) = tree.get("FoamFile") else {
        return Err("no FoamFile header".into());
    };
    for key in HEADER_KEYS {
        if header.get(key).is_none() {
            return Err(format!("FoamFile header has no `{key}` entry"));
        }
    }
    if header.text("format").as_deref() != Some("ascii") {
        return Err("FoamFile header does not say `format ascii`".into());
    }
    Ok(tree)
}

/// Split the text into tokens: whitespace and comments are dropped, a word keeps its balanced
/// parentheses, a number must be finite.
fn lex(text: &str) -> Result<Vec<Lexed>, String> {
    let chars: Vec<(usize, char)> = text.char_indices().collect();
    let char_at = |i: usize| chars.get(i).map(|&(_, c)| c);
    let mut out = Vec::new();
    let mut line = 1;
    let mut i = 0;
    while let Some(&(start, c)) = chars.get(i) {
        if c == '\n' {
            line += 1;
            i += 1;
        } else if c.is_whitespace() {
            i += 1;
        } else if c == '/' && char_at(i + 1) == Some('/') {
            while char_at(i).is_some_and(|c| c != '\n') {
                i += 1;
            }
        } else if c == '/' && char_at(i + 1) == Some('*') {
            let open_line = line;
            i += 2;
            loop {
                match char_at(i) {
                    None => return Err(format!("line {open_line}: unterminated `/*` comment")),
                    Some('*') if char_at(i + 1) == Some('/') => {
                        i += 2;
                        break;
                    }
                    Some(ch) => {
                        if ch == '\n' {
                            line += 1;
                        }
                        i += 1;
                    }
                }
            }
        } else if "{}()[];".contains(c) {
            out.push(Lexed { tok: Token::Punct(c), line, start });
            i += 1;
        } else if c == '"' {
            let open_line = line;
            let mut s = String::new();
            i += 1;
            loop {
                match char_at(i) {
                    None => return Err(format!("line {open_line}: unterminated string")),
                    Some('"') => break,
                    Some(ch) => {
                        if ch == '\n' {
                            line += 1;
                        }
                        s.push(ch);
                        i += 1;
                    }
                }
            }
            out.push(Lexed { tok: Token::Str(s), line: open_line, start });
            i += 1;
        } else if c.is_ascii_alphabetic() || matches!(c, '_' | '#' | '$') {
            let mut s = String::new();
            let mut depth = 0usize;
            while let Some(ch) = char_at(i) {
                if ch.is_whitespace() || matches!(ch, ';' | '{' | '}' | '"' | '[' | ']') {
                    break;
                }
                if ch == '(' {
                    depth += 1;
                } else if ch == ')' {
                    if depth == 0 {
                        break;
                    }
                    depth -= 1;
                }
                s.push(ch);
                i += 1;
            }
            out.push(Lexed { tok: Token::Word(s), line, start });
        } else if c.is_ascii_digit() || matches!(c, '-' | '+' | '.') {
            let mut s = String::new();
            while let Some(ch) = char_at(i) {
                if ch.is_whitespace() || "{}()[];\"".contains(ch) {
                    break;
                }
                s.push(ch);
                i += 1;
            }
            let n: f64 = s.parse().map_err(|_| format!("line {line}: `{s}` is not a number"))?;
            if !n.is_finite() {
                return Err(format!("line {line}: `{s}` is not a finite number"));
            }
            out.push(Lexed { tok: Token::Number(n), line, start });
        } else {
            return Err(format!("line {line}: unexpected character `{c}`"));
        }
    }
    Ok(out)
}

struct Parser<'a> {
    toks: Vec<Lexed>,
    pos: usize,
    text: &'a str,
}

impl Parser<'_> {
    /// The line of the next token (or of the last one at end of file).
    fn line_here(&self) -> usize {
        self.toks.get(self.pos).or(self.toks.last()).map_or(1, |t| t.line)
    }

    fn error(&self, what: &str) -> String {
        format!("line {}: {what}", self.line_here())
    }

    /// The entries of a file (`in_block` false) or of a `{ }` block (`in_block` true).
    fn entries(&mut self, depth: usize, in_block: bool) -> Result<Dict, String> {
        let mut dict = Dict::default();
        loop {
            let Some(next) = self.toks.get(self.pos).cloned() else {
                return if in_block { Err(self.error("missing `}` at end of file")) } else { Ok(dict) };
            };
            match next.tok {
                Token::Punct('}') => {
                    if !in_block {
                        return Err(self.error("`}` closes nothing"));
                    }
                    self.pos += 1;
                    return Ok(dict);
                }
                Token::Word(ref w) if w.starts_with('#') => {
                    // A directive runs to the end of its line and has no `;`.
                    let rest = self.text.get(next.start..).unwrap_or_default();
                    let end = rest.find('\n').unwrap_or(rest.len());
                    dict.directives.push(rest.get(..end).unwrap_or_default().trim_end().to_string());
                    while self.toks.get(self.pos).is_some_and(|t| t.line == next.line) {
                        self.pos += 1;
                    }
                    continue;
                }
                _ => {}
            }
            let key = match next.tok {
                Token::Word(w) | Token::Str(w) => w,
                _ => return Err(self.error("expected a keyword")),
            };
            self.pos += 1;
            let entry = if matches!(self.toks.get(self.pos).map(|t| &t.tok), Some(Token::Punct('{'))) {
                if depth >= MAX_DEPTH {
                    return Err(self.error("dictionaries nested too deeply"));
                }
                self.pos += 1;
                Entry::Dict(self.entries(depth + 1, true)?)
            } else {
                Entry::Value(self.value()?)
            };
            dict.entries.push((key, entry));
        }
    }

    /// The tokens of a primitive entry, up to the `;` that is outside every bracket (possibly none).
    fn value(&mut self) -> Result<Vec<Token>, String> {
        let mut out = Vec::new();
        // The closing character each open bracket still waits for.
        let mut closers: Vec<char> = Vec::new();
        loop {
            let Some(next) = self.toks.get(self.pos).cloned() else {
                return Err(self.error("missing `;` at end of file"));
            };
            self.pos += 1;
            match &next.tok {
                // A bare `$name;` (macro substitution) has no value tokens at all.
                Token::Punct(';') if closers.is_empty() => return Ok(out),
                Token::Punct(c) if matches!(*c, '(' | '[' | '{') => {
                    if closers.len() >= MAX_DEPTH {
                        return Err(self.error("brackets nested too deeply"));
                    }
                    closers.push(match *c {
                        '(' => ')',
                        '[' => ']',
                        _ => '}',
                    });
                }
                Token::Punct(c) if matches!(*c, ')' | ']' | '}') => match closers.pop() {
                    Some(open) if open == *c => {}
                    _ => return Err(self.error(&format!("`{c}` does not close the open bracket (missing `;`?)"))),
                },
                _ => {}
            }
            out.push(next.tok);
        }
    }
}

fn render(token: &Token) -> String {
    match token {
        Token::Word(w) => w.clone(),
        Token::Number(n) => format!("{n}"),
        Token::Str(s) => format!("\"{s}\""),
        Token::Punct(c) => c.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nested_dictionaries_lists_and_comments() {
        let text = r#"
FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
// a comment
/* a block
   comment */
divSchemes { default none; div(phi,U) Gauss linearUpwind grad(U); div((nuEff*dev2(T(grad(U))))) Gauss linear; }
solvers { "(U|UFinal)" { solver smoothSolver; tolerance 1e-6; } }
features ( { file "a.eMesh"; level 0; } );
dimensions [0 1 -1 0 0 0 0];
internalField uniform (2e0 0 0);
"#;
        let tree = check_file(text).unwrap();
        assert_eq!(tree.text("divSchemes/div(phi,U)").as_deref(), Some("Gauss linearUpwind grad(U)"));
        assert_eq!(tree.number("solvers/(U|UFinal)/tolerance"), Some(1e-6));
        assert_eq!(tree.text("internalField").as_deref(), Some("uniform ( 2 0 0 )"));
        assert!(matches!(tree.get("divSchemes/div((nuEff*dev2(T(grad(U)))))"), Some(Entry::Value(_))));
        assert!(tree.all_directives().is_empty());
    }

    #[test]
    fn records_include_directives() {
        let text = "FoamFile { version 2.0; format ascii; class dictionary; object a; }\n\
                    m { #includeEtc \"caseDicts/meshQualityDict\"\n nSmoothScale 4; }\n";
        let tree = parse(text).unwrap();
        assert_eq!(tree.all_directives(), vec!["#includeEtc \"caseDicts/meshQualityDict\"".to_string()]);
        assert_eq!(tree.number("m/nSmoothScale"), Some(4.0));
    }

    #[test]
    fn refuses_broken_files() {
        let deep = "a { ".repeat(100);
        for text in ["a 1", "a { b 1; ", "}", "a (1 2;", "/* open", "a \"open;", "a 1e999;", deep.as_str()] {
            assert!(parse(text).is_err(), "accepted: {text}");
        }
    }

    #[test]
    fn check_file_needs_the_header() {
        let no_header = check_file("a 1;").unwrap_err();
        assert!(no_header.contains("FoamFile"), "{no_header}");
        let no_object = check_file("FoamFile { version 2.0; format ascii; class dictionary; }").unwrap_err();
        assert!(no_object.contains("object"), "{no_object}");
    }
}
