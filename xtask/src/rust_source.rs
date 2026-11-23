//! The workspace's Rust sources read as syntax trees, for the source checks: `cargo xtask
//! unit-factors` (literal conversion factors) and `cargo xtask libm` (std transcendental methods).
//!
//! A check sees each file's float literals and its calls to a float maths method, with their
//! lines. Comments and strings are not code, and test code (`#[cfg(test)]` items, `#[test]`
//! functions, the module files they declare, and the `tests`/`benches` directories) is not read at
//! all: test code never feeds a result the solvers produce.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use proc_macro2::{Delimiter, TokenStream, TokenTree};
use syn::visit::{self, Visit};
use syn::{Attribute, Expr, ExprMethodCall, ExprPath, ImplItem, Item, Lit, Meta, TraitItem};

/// Directories never read: build output, VCS, planning notes, and code that is not shipped.
const SKIP_DIRS: &[&str] = &["target", ".git", "plan", "dist", "xtask", "tests", "benches", "node_modules"];

/// Float methods that evaluate a transcendental function through the platform's maths library, so
/// their last bits depend on it. Not listed, because they never call that library: `sqrt`, `abs`
/// and `mul_add` are correctly rounded by IEEE 754, and `powi` is a sequence of IEEE multiplies
/// (LLVM's expansion of a constant exponent, and compiler-builtins' `__powidf2` otherwise).
pub const FLOAT_METHODS: &[&str] = &[
    "exp", "exp_m1", "ln", "ln_1p", "log", "log10", "log2", "powf", "sin", "cos", "sin_cos", "tan", "asin", "acos", "atan", "atan2", "sinh", "cosh",
    "tanh", "asinh", "acosh", "atanh", "cbrt", "hypot",
];

/// One source file: its path relative to the root, its text, and what it holds outside test code.
pub struct Source {
    pub rel: String,
    pub text: String,
    /// Float literals (line, value) in code and macro bodies. Integers that are counts are not.
    pub literals: Vec<(usize, f64)>,
    /// Calls to a [`FLOAT_METHODS`] entry: `x.exp()` as `exp`, `f64::exp` as `f64::exp`.
    pub calls: Vec<(usize, String)>,
}

/// The reason of a line's marker (empty when it gives none), if it has one.
pub fn marker_reason<'a>(line: &'a str, marker: &str) -> Option<&'a str> {
    line.find(marker).map(|i| line.get(i + marker.len()..).unwrap_or("").trim())
}

/// Reads every non-test Rust file under `root/<dir>` for each `dir`, sorted by path. Returns the
/// files and the files that do not parse (one message each).
pub fn read_tree(root: &Path, dirs: &[&str]) -> Result<(Vec<Source>, Vec<String>), String> {
    let mut files = Vec::new();
    for dir in dirs {
        walk(&root.join(dir), &mut files);
    }
    files.sort();

    let mut parsed = Vec::new();
    let mut tests = TestPaths::default();
    let mut errors = Vec::new();
    for path in &files {
        let rel = relative(root, path);
        let text = std::fs::read_to_string(path).map_err(|e| format!("{rel}: {e}"))?;
        let ast = match syn::parse_file(&text) {
            Ok(ast) => ast,
            Err(e) => {
                let at = e.span().start();
                errors.push(format!("{rel}:{}:{}: cannot be parsed: {e}", at.line, at.column + 1));
                continue;
            }
        };
        let mut scan = Scan { file: path, inline: Vec::new(), literals: Vec::new(), calls: Vec::new(), tests: &mut tests };
        scan.visit_file(&ast);
        parsed.push((path.clone(), Source { rel, text, literals: scan.literals, calls: scan.calls }));
    }

    let sources = parsed.into_iter().filter(|(path, _)| !tests.contains(path)).map(|(_, source)| source).collect();
    Ok((sources, errors))
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root).unwrap_or(path).to_string_lossy().replace('\\', "/")
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

/// One file's literals and float-method calls outside test code.
struct Scan<'a> {
    file: &'a Path,
    /// The inline modules (`mod a { … }`) the visitor is in.
    inline: Vec<String>,
    literals: Vec<(usize, f64)>,
    calls: Vec<(usize, String)>,
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
            self.literals.push((lit.span().start().line, v));
        }
    }

    /// Macro bodies are tokens to syn: literals are read there, and so is `.name(…)`.
    fn tokens(&mut self, tokens: &TokenStream) {
        let trees: Vec<TokenTree> = tokens.clone().into_iter().collect();
        for (i, tree) in trees.iter().enumerate() {
            match tree {
                TokenTree::Literal(l) => self.literal(&Lit::new(l.clone())),
                TokenTree::Group(g) => self.tokens(&g.stream()),
                TokenTree::Ident(name) => {
                    let after_dot =
                        i.checked_sub(1).and_then(|j| trees.get(j)).is_some_and(|t| matches!(t, TokenTree::Punct(p) if p.as_char() == '.'));
                    let called = matches!(trees.get(i + 1), Some(TokenTree::Group(g)) if g.delimiter() == Delimiter::Parenthesis);
                    let name = name.to_string();
                    if after_dot && called && FLOAT_METHODS.contains(&name.as_str()) {
                        self.calls.push((tree.span().start().line, name));
                    }
                }
                TokenTree::Punct(_) => {}
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

    fn visit_expr_method_call(&mut self, node: &'ast ExprMethodCall) {
        let name = node.method.to_string();
        if FLOAT_METHODS.contains(&name.as_str()) {
            self.calls.push((node.method.span().start().line, name));
        }
        visit::visit_expr_method_call(self, node);
    }

    /// `f64::exp` as a path, called or passed as a function value (`map(f64::sin)`).
    fn visit_expr_path(&mut self, node: &'ast ExprPath) {
        let segments: Vec<&syn::Ident> = node.path.segments.iter().map(|s| &s.ident).collect();
        if let [.., ty, name] = segments.as_slice()
            && (*ty == "f64" || *ty == "f32")
            && FLOAT_METHODS.contains(&name.to_string().as_str())
        {
            self.calls.push((name.span().start().line, format!("{ty}::{name}")));
        }
        visit::visit_expr_path(self, node);
    }

    /// Macro bodies (`json!`, `format!`, `vec!`) are tokens to syn: their literals and calls are read there.
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
            Expr::Lit(syn::ExprLit { lit: Lit::Str(s), .. }) => Some(s.value()),
            _ => None,
        },
        _ => None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scan(file: &str, code: &str) -> (Source, TestPaths) {
        let mut tests = TestPaths::default();
        let ast = syn::parse_file(code).unwrap();
        let mut s = Scan { file: Path::new(file), inline: Vec::new(), literals: Vec::new(), calls: Vec::new(), tests: &mut tests };
        s.visit_file(&ast);
        let source = Source { rel: file.to_string(), text: code.to_string(), literals: s.literals, calls: s.calls };
        (source, tests)
    }

    #[test]
    fn literals_in_code_and_macros_not_in_strings_or_comments() {
        let code = "fn f(y: f64) -> f64 {\n    let x = 1e-3 * y + 273.15_f64; // 25.4\n    let s = \"1e6\";\n    \
                    println!(\"{}\", 9810.0 + 3f64);\n    x\n}\n";
        let (source, _) = scan("a/src/lib.rs", code);
        assert_eq!(source.literals, vec![(2, 1e-3), (2, 273.15), (4, 9810.0), (4, 3.0)]);
    }

    #[test]
    fn test_code_is_skipped_and_its_module_files_found() {
        let code = "#[cfg(test)]\nmod tests;\n#[cfg(all(test, unix))]\nfn t() { 1e-3; }\n\
                    #[cfg(any(test, feature = \"x\"))]\nfn kept() { 1e6; }\n#[cfg(test)]\n#[path = \"x_tests.rs\"]\nmod x;\n\
                    mod inner {\n    #[cfg(test)]\n    mod deep;\n    #[test]\n    fn t() { 1e9; }\n}\n";
        let (source, tests) = scan("a/src/b.rs", code);
        assert_eq!(source.literals, vec![(6, 1e6)]);
        assert!(tests.contains(Path::new("a/src/b/tests.rs")));
        assert!(tests.contains(Path::new("a/src/b/tests/more.rs")));
        assert!(tests.contains(Path::new("a/src/x_tests.rs")));
        assert!(tests.contains(Path::new("a/src/b/inner/deep.rs")));
        assert!(!tests.contains(Path::new("a/src/b/inner.rs")));
    }

    #[test]
    fn float_method_calls_in_code_macros_and_paths() {
        let code = "fn f(x: f64, n: usize) -> f64 {\n    let a = x.exp() + x.sqrt() + f64::sin(x);\n    let v: Vec<f64> = [x].iter().copied().map(f64::cos).collect();\n    \
                    println!(\"{}\", x.ln());\n    let n = n.pow(2);\n    a + v[0]\n}\n";
        let (source, _) = scan("a/src/lib.rs", code);
        let names: Vec<(usize, &str)> = source.calls.iter().map(|(l, n)| (*l, n.as_str())).collect();
        assert_eq!(names, vec![(2, "exp"), (2, "f64::sin"), (3, "f64::cos"), (4, "ln")]);
    }
}
