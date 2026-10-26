//! `cargo xtask book`: build the documentation site from `docs/` with mdBook.
//!
//! `docs/` is the book's source as it stands (`docs/book.toml`, `docs/SUMMARY.md`), so every page
//! still reads on GitHub. Pages link to files outside `docs/` (`../tests/oracle/…`,
//! `../crates/…`) that the site does not carry: the build copies `docs/` to `target/book-src`,
//! points those links at the repository on GitHub, and runs `mdbook` (a separate tool, not a
//! dependency) into `target/book`. After SolidCraft's `xtask book` (MIT OR Apache-2.0; see
//! ATTRIBUTION.md).

use std::path::{Component, Path, PathBuf};

/// Where links that leave `docs/` go: the trunk on GitHub (a directory's `blob` URL shows its tree).
const REPOSITORY: &str = "https://github.com/tangency-labs/fieldcraft/blob/develop/";

pub fn run(root: &Path, build: bool) -> Result<(), String> {
    let docs = root.join("docs");
    let staged = root.join("target/book-src");
    if staged.exists() {
        std::fs::remove_dir_all(&staged).map_err(|e| format!("{}: {e}", staged.display()))?;
    }
    copy_tree(&docs, &staged, &docs)?;
    if build {
        let mut c = std::process::Command::new("mdbook");
        c.arg("build").arg(&staged).arg("--dest-dir").arg(root.join("target/book"));
        crate::run(c, "mdbook build").map_err(|e| format!("{e}\n(install mdbook with `cargo install mdbook --version 0.5.4`)"))?;
        println!("book: built into target/book");
    }
    Ok(())
}

/// Copies `from` into `to`, rewriting the links of each Markdown file; `docs` is the source root.
fn copy_tree(from: &Path, to: &Path, docs: &Path) -> Result<(), String> {
    std::fs::create_dir_all(to).map_err(|e| format!("{}: {e}", to.display()))?;
    let entries = std::fs::read_dir(from).map_err(|e| format!("{}: {e}", from.display()))?;
    for entry in entries.flatten() {
        let (src, dst) = (entry.path(), to.join(entry.file_name()));
        if src.is_dir() {
            copy_tree(&src, &dst, docs)?;
        } else if src.extension().is_some_and(|x| x == "md") {
            let text = std::fs::read_to_string(&src).map_err(|e| format!("{}: {e}", src.display()))?;
            let dir = src.parent().and_then(|d| d.strip_prefix(docs).ok()).unwrap_or(Path::new(""));
            std::fs::write(&dst, rewrite(&text, dir)).map_err(|e| format!("{}: {e}", dst.display()))?;
        } else {
            std::fs::copy(&src, &dst).map_err(|e| format!("{}: {e}", src.display()))?;
        }
    }
    Ok(())
}

/// The page `text` in `dir` (relative to `docs/`) with every Markdown link target `](…)` that
/// leaves `docs/` pointed at the repository. Fenced code blocks are left as they are.
fn rewrite(text: &str, dir: &Path) -> String {
    let mut out = String::with_capacity(text.len());
    let mut fenced = false;
    for line in text.split_inclusive('\n') {
        if line.trim_start().starts_with("```") {
            fenced = !fenced;
        }
        if fenced {
            out.push_str(line);
            continue;
        }
        let mut rest = line;
        while let Some(i) = rest.find("](") {
            let (head, tail) = rest.split_at(i + 2);
            out.push_str(head);
            let end = tail.find(')').unwrap_or(tail.len());
            let (target, after) = tail.split_at(end);
            out.push_str(&outside(target, dir).unwrap_or_else(|| target.to_string()));
            rest = after;
        }
        out.push_str(rest);
    }
    out
}

/// The repository URL of a relative link `target` from `dir`, when it leaves `docs/`.
fn outside(target: &str, dir: &Path) -> Option<String> {
    if target.is_empty() || target.starts_with('#') || target.starts_with('/') || target.contains("://") || target.starts_with("mailto:") {
        return None;
    }
    let (path, anchor) = target.split_once('#').map_or((target, None), |(p, a)| (p, Some(a)));
    // The target's path from the repository root; one that climbs above the root is left alone.
    let mut repo: Vec<String> = vec!["docs".into()];
    for c in dir.join(path).components() {
        match c {
            Component::ParentDir => {
                repo.pop()?;
            }
            Component::CurDir => {}
            Component::Normal(s) => repo.push(s.to_string_lossy().into_owned()),
            _ => return None,
        }
    }
    if repo.first().is_some_and(|s| s == "docs") {
        return None;
    }
    let repo_path: PathBuf = repo.iter().collect();
    let mut url = format!("{REPOSITORY}{}", repo_path.to_string_lossy());
    if let Some(a) = anchor {
        url.push('#');
        url.push_str(a);
    }
    Some(url)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links_inside_docs_stay_and_links_outside_point_at_the_repository() {
        let page = "See [the policy](verification.md), [cases](../tests/oracle/), \
                    [units](../crates/engine/src/units.rs#L3) and [the web](https://example.org).\n";
        assert_eq!(
            rewrite(page, Path::new("")),
            "See [the policy](verification.md), \
             [cases](https://github.com/tangency-labs/fieldcraft/blob/develop/tests/oracle), \
             [units](https://github.com/tangency-labs/fieldcraft/blob/develop/crates/engine/src/units.rs#L3) \
             and [the web](https://example.org).\n"
        );
    }

    #[test]
    fn links_from_a_subfolder_resolve_from_it() {
        assert_eq!(
            rewrite("[skill](../../skills/fieldcraft-meshing/SKILL.md) [pic](skill/a.png) [up](../units.md)", Path::new("meshing")),
            "[skill](https://github.com/tangency-labs/fieldcraft/blob/develop/skills/fieldcraft-meshing/SKILL.md) [pic](skill/a.png) [up](../units.md)"
        );
    }

    #[test]
    fn fenced_code_is_left_alone() {
        let page = "```\n[x](../a.md)\n```\n[y](../b.md)\n";
        assert_eq!(rewrite(page, Path::new("")), "```\n[x](../a.md)\n```\n[y](https://github.com/tangency-labs/fieldcraft/blob/develop/b.md)\n");
    }
}
