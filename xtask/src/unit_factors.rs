//! `cargo xtask unit-factors`: literal unit-conversion factors (1e-3, 1e6, 273.15, 25.4, …) written
//! outside the unit table (`crates/model/src/units/`), a check of `cargo xtask ci` (docs/units.md).
//!
//! Every conversion goes through the table (#579), so a float literal equal to a common factor is
//! taken for a hand-written conversion. A number that only looks like one, a tolerance or a cap,
//! is a named constant whose line carries the marker `// unit-factors: <why it is not a
//! conversion>`. The sources are read as syntax trees: comments, strings (help texts, WGSL) and
//! test code (`#[cfg(test)]` items and the module files they declare) are not code that converts.
//!
//! One part of the tree is not counted, listed here with its reason: the frozen migration step
//! ([`FROZEN`]).

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use proc_macro2::{TokenStream, TokenTree};
use syn::visit::{self, Visit};
use syn::{Attribute, ImplItem, Item, Lit, Meta, TraitItem};

/// The factors flagged, with what each usually converts.
const FACTORS: &[(f64, &str)] = &[
    (1e-3, "milli (mm↔m, ms↔s, g↔kg)"),
    (1e3, "kilo (t↔kg, m↔mm)"),
    (1e6, "mega (MPa↔Pa, mm²↔m²)"),
    (1e9, "giga (GPa↔Pa, mm³↔m³)"),
    (1e-6, "micro (Pa↔MPa)"),
    (1e-9, "nano (Pa↔GPa, m³↔mm³)"),
    (1e12, "t/mm³↔kg/m³"),
    (1e-12, "kg/m³↔t/mm³"),
    (273.15, "°C↔K"),
    (459.67, "°F↔°R"),
    (25.4, "in↔mm"),
    (0.0254, "in↔m"),
    (9.80665, "standard gravity"),
    (9.81, "gravity"),
    (9810.0, "gravity in mm/s²"),
];

/// The one place factors belong.
const TABLE: &str = "crates/model/src/units/";

/// Files that keep their literal factors, and why. Nothing else is exempt as a whole.
const FROZEN: &[(&str, &str)] = &[(
    "crates/io/src/migrate/v2_to_v3.rs",
    "the format 2 → 3 migration step is frozen: it converts old files with its own copy of the \
     factors, so a later change to the table cannot change how an old file migrates",
)];

/// The per-line exemption: `// unit-factors: <reason>` on the line of the literal.
const MARKER: &str = "// unit-factors:";

const SKIP_DIRS: &[&str] = &["target", ".git", "plan", "dist", "xtask", "tests", "benches", "node_modules"];

pub fn run(root: &Path) -> Result<(), String> {
    let mut files = Vec::new();
    for dir in ["crates", "apps"] {
        walk(&root.join(dir), &mut files);
    }
    files.sort();

    let mut parsed = Vec::new();
    let mut tests = TestPaths::default();
    let mut errors = Vec::new();
    for path in &files {
        let rel = relative(root, path);
        if rel.starts_with(TABLE) {
            continue;
        }
        let text = std::fs::read_to_string(path).map_err(|e| format!("{rel}: {e}"))?;
        let ast = match syn::parse_file(&text) {
            Ok(ast) => ast,
            Err(e) => {
                let at = e.span().start();
                errors.push(format!("{rel}:{}:{}: cannot be parsed: {e}", at.line, at.column + 1));
                continue;
            }
        };
        let mut scan = Scan { file: path, inline: Vec::new(), found: Vec::new(), tests: &mut tests };
        scan.visit_file(&ast);
        let found = scan.found;
        parsed.push((path.clone(), rel, text, found));
    }

    let mut checked = Vec::new();
    let mut frozen = 0usize;
    for (path, rel, text, found) in &parsed {
        if tests.contains(path) {
            continue;
        }
        if FROZEN.iter().any(|(f, _)| f == rel) {
            frozen += found.iter().filter(|(_, v)| factor(*v).is_some()).count();
            continue;
        }
        let lines: Vec<&str> = text.lines().collect();
        let mut flagged = BTreeSet::new();
        for (line, value) in found {
            let Some(what) = factor(*value) else { continue };
            flagged.insert(*line);
            let source = lines.get(line.wrapping_sub(1)).copied().unwrap_or("");
            match marker(source) {
                Some("") => errors.push(format!("{rel}:{line}: `{MARKER}` needs a reason")),
                Some(_) => {}
                None => {
                    let hit = format!("{rel}:{line}: {value:e} ({what}): {}", source.trim());
                    checked.push(hit);
                }
            }
        }
        for (n, source) in lines.iter().enumerate() {
            if marker(source).is_some() && !flagged.contains(&(n + 1)) {
                errors.push(format!("{rel}:{}: `{MARKER}` on a line without a factor; remove it", n + 1));
            }
        }
    }

    for h in &checked {
        println!("{h}");
    }
    for e in &errors {
        println!("{e}");
    }
    for (file, why) in FROZEN {
        println!("unit-factors: {file} keeps its factors: {why}");
    }
    println!("unit-factors: {} literal factor(s) outside {TABLE}; {frozen} in the frozen step", checked.len());
    if checked.is_empty() && errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "unit-factors: {} literal factor(s) and {} marker error(s): convert through {TABLE}, or name the constant and \
             mark its line `{MARKER} <why it is not a conversion>` (docs/units.md)",
            checked.len(),
            errors.len()
        ))
    }
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root).unwrap_or(path).to_string_lossy().replace('\\', "/")
}

fn factor(value: f64) -> Option<&'static str> {
    FACTORS.iter().find(|(f, _)| *f == value).map(|(_, what)| *what)
}

/// The reason of a line's marker (empty when it gives none), if it has one.
fn marker(line: &str) -> Option<&str> {
    line.find(MARKER).map(|i| line.get(i + MARKER.len()..).unwrap_or("").trim())
}

fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = std::fs::read_dir(dir) else { return };
    for e in rd.flatten() {
        let p = e.path();
        let name = e.file_name().to_string_lossy().to_string();
        if p.is_dir() {
            if !SKIP_DIRS.contains(&name.as_str()) && !name.starts_with('.') {
                walk(&p, out);
            }
        } else if p.extension().is_some_and(|x| x == "rs") {
            out.push(p);
        }
    }
}

/// The module files and directories that are compiled only for tests.
#[derive(Default)]
struct TestPaths {
    files: BTreeSet<PathBuf>,
    dirs: BTreeSet<PathBuf>,
}

impl TestPaths {
    fn contains(&self, path: &Path) -> bool {
        self.files.contains(path) || self.dirs.iter().any(|d| path.starts_with(d))
    }
}

/// One file's float literals (line, value) outside test code.
struct Scan<'a> {
    file: &'a Path,
    /// The inline modules (`mod a { … }`) the visitor is in.
    inline: Vec<String>,
    found: Vec<(usize, f64)>,
    tests: &'a mut TestPaths,
}

impl Scan<'_> {
    /// The directory of the modules `mod x;` declares here.
    fn module_dir(&self) -> PathBuf {
        let parent = self.file.parent().unwrap_or(Path::new("")).to_path_buf();
        let mod_rs = self.file.file_name().is_some_and(|n| n == "lib.rs" || n == "main.rs" || n == "mod.rs");
        let mut dir = if mod_rs { parent } else { parent.join(self.file.file_stem().unwrap_or_default()) };
        for m in &self.inline {
            dir.push(m);
        }
        dir
    }

    /// A test-only module: its file (or inline body) and the directory of its own submodules.
    fn test_module(&mut self, m: &syn::ItemMod) {
        let name = m.ident.to_string();
        let dir = self.module_dir();
        if m.content.is_none() {
            let file = match path_attr(&m.attrs) {
                // Outside inline modules `#[path]` is relative to the declaring file's directory.
                Some(p) if self.inline.is_empty() => self.file.parent().unwrap_or(Path::new("")).join(p),
                Some(p) => dir.join(p),
                None => dir.join(format!("{name}.rs")),
            };
            self.tests.files.insert(file);
        }
        self.tests.dirs.insert(dir.join(name));
    }

    /// Conversions are floating-point arithmetic: float literals, and integers typed as floats
    /// (`3f64`). Plain integer literals are counts and caps.
    fn literal(&mut self, lit: &Lit) {
        let value = match lit {
            Lit::Float(f) => f.base10_parse::<f64>().ok(),
            Lit::Int(i) if i.suffix().starts_with('f') => i.base10_digits().parse::<f64>().ok(),
            _ => None,
        };
        if let Some(v) = value {
            self.found.push((lit.span().start().line, v));
        }
    }

    fn tokens(&mut self, tokens: &TokenStream) {
        for t in tokens.clone() {
            match t {
                TokenTree::Literal(l) => self.literal(&Lit::new(l)),
                TokenTree::Group(g) => self.tokens(&g.stream()),
                TokenTree::Ident(_) | TokenTree::Punct(_) => {}
            }
        }
    }
}

impl<'ast> Visit<'ast> for Scan<'_> {
    fn visit_item(&mut self, item: &'ast Item) {
        if test_only(item_attrs(item)) {
            if let Item::Mod(m) = item {
                self.test_module(m);
            }
            return;
        }
        if let Item::Mod(m) = item
            && m.content.is_some()
        {
            self.inline.push(m.ident.to_string());
            visit::visit_item(self, item);
            self.inline.pop();
            return;
        }
        visit::visit_item(self, item);
    }

    fn visit_impl_item(&mut self, item: &'ast ImplItem) {
        let attrs = match item {
            ImplItem::Const(i) => &i.attrs,
            ImplItem::Fn(i) => &i.attrs,
            ImplItem::Type(i) => &i.attrs,
            ImplItem::Macro(i) => &i.attrs,
            _ => return visit::visit_impl_item(self, item),
        };
        if !test_only(attrs) {
            visit::visit_impl_item(self, item);
        }
    }

    fn visit_trait_item(&mut self, item: &'ast TraitItem) {
        let attrs = match item {
            TraitItem::Const(i) => &i.attrs,
            TraitItem::Fn(i) => &i.attrs,
            TraitItem::Type(i) => &i.attrs,
            TraitItem::Macro(i) => &i.attrs,
            _ => return visit::visit_trait_item(self, item),
        };
        if !test_only(attrs) {
            visit::visit_trait_item(self, item);
        }
    }

    fn visit_lit(&mut self, lit: &'ast Lit) {
        self.literal(lit);
    }

    /// Macro bodies (`json!`, `format!`, `vec!`) are tokens to syn: their literals are read there.
    fn visit_macro(&mut self, mac: &'ast syn::Macro) {
        self.tokens(&mac.tokens);
    }
}

fn item_attrs(item: &Item) -> &[Attribute] {
    match item {
        Item::Const(i) => &i.attrs,
        Item::Enum(i) => &i.attrs,
        Item::ExternCrate(i) => &i.attrs,
        Item::Fn(i) => &i.attrs,
        Item::ForeignMod(i) => &i.attrs,
        Item::Impl(i) => &i.attrs,
        Item::Macro(i) => &i.attrs,
        Item::Mod(i) => &i.attrs,
        Item::Static(i) => &i.attrs,
        Item::Struct(i) => &i.attrs,
        Item::Trait(i) => &i.attrs,
        Item::TraitAlias(i) => &i.attrs,
        Item::Type(i) => &i.attrs,
        Item::Union(i) => &i.attrs,
        Item::Use(i) => &i.attrs,
        _ => &[],
    }
}

/// Whether the attributes compile an item only for tests: `#[test]`, `#[cfg(test)]` or
/// `#[cfg(all(test, …))]` (`any(test, …)` is not test-only).
fn test_only(attrs: &[Attribute]) -> bool {
    attrs.iter().any(|a| a.path().is_ident("test") || (a.path().is_ident("cfg") && a.parse_args::<Meta>().is_ok_and(|m| requires_test(&m))))
}

fn requires_test(m: &Meta) -> bool {
    match m {
        Meta::Path(p) => p.is_ident("test"),
        Meta::List(l) if l.path.is_ident("all") => {
            l.parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated).is_ok_and(|all| all.iter().any(requires_test))
        }
        _ => false,
    }
}

/// A module's `#[path = "…"]`.
fn path_attr(attrs: &[Attribute]) -> Option<String> {
    attrs.iter().find_map(|a| match &a.meta {
        Meta::NameValue(nv) if nv.path.is_ident("path") => match &nv.value {
            syn::Expr::Lit(syn::ExprLit { lit: Lit::Str(s), .. }) => Some(s.value()),
            _ => None,
        },
        _ => None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scan(file: &str, code: &str) -> (Vec<(usize, f64)>, TestPaths) {
        let mut tests = TestPaths::default();
        let ast = syn::parse_file(code).unwrap();
        let mut s = Scan { file: Path::new(file), inline: Vec::new(), found: Vec::new(), tests: &mut tests };
        s.visit_file(&ast);
        let found = s.found;
        (found, tests)
    }

    #[test]
    fn literals_in_code_and_macros_not_in_strings_or_comments() {
        let code = "fn f(y: f64) -> f64 {\n    let x = 1e-3 * y + 273.15_f64; // 25.4\n    let s = \"1e6\";\n    \
                    println!(\"{}\", 9810.0 + 3f64);\n    x\n}\n";
        let (found, _) = scan("a/src/lib.rs", code);
        assert_eq!(found, vec![(2, 1e-3), (2, 273.15), (4, 9810.0), (4, 3.0)]);
    }

    #[test]
    fn test_code_is_skipped_and_its_module_files_found() {
        let code = "#[cfg(test)]\nmod tests;\n#[cfg(all(test, unix))]\nfn t() { 1e-3; }\n\
                    #[cfg(any(test, feature = \"x\"))]\nfn kept() { 1e6; }\n#[cfg(test)]\n#[path = \"x_tests.rs\"]\nmod x;\n\
                    mod inner {\n    #[cfg(test)]\n    mod deep;\n    #[test]\n    fn t() { 1e9; }\n}\n";
        let (found, tests) = scan("a/src/b.rs", code);
        assert_eq!(found, vec![(6, 1e6)]);
        assert!(tests.contains(Path::new("a/src/b/tests.rs")));
        assert!(tests.contains(Path::new("a/src/b/tests/more.rs")));
        assert!(tests.contains(Path::new("a/src/x_tests.rs")));
        assert!(tests.contains(Path::new("a/src/b/inner/deep.rs")));
        assert!(!tests.contains(Path::new("a/src/b/inner.rs")));
    }

    #[test]
    fn markers() {
        assert_eq!(marker("const X: f64 = 1e-9; // unit-factors: relative to the size"), Some("relative to the size"));
        assert_eq!(marker("let x = 1e-9; // unit-factors:"), Some(""));
        assert_eq!(marker("let x = 1e-9;"), None);
    }
}
