//! Bounded INCLUDE expansion with an immutable original file bundle.
//! Abaqus Input Syntax Rules: relative INPUT paths use the job submission directory.
//! https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-inputsyntax.htm
use super::{MAX_INPUT_BYTES, MAX_RECORDS, lex, source::SourceDocument};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

const MAX_FILES: usize = 1024;
const MAX_DEPTH: usize = 32;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Location {
    pub file: String,
    pub line: usize,
}

/// Original files plus a derived stream for interpretation. No file system access in decks.
#[derive(Debug, Clone)]
pub struct SourceBundle {
    pub root: String,
    pub files: BTreeMap<String, SourceDocument>,
    pub expanded: String,
    pub locations: Vec<Location>,
}

/// Portable relative paths only, so a bundle can also be exported without escaping its folder.
pub fn relative_path(value: &str, line: usize) -> Result<String> {
    let value = lex::unquote(value, line)?;
    if value.is_empty() || value.starts_with('/') || value.contains(['\\', ':']) || value.split('/').any(|s| s.is_empty() || s == ".." || s == ".") {
        return Err(lex::error(line, "INCLUDE requires a relative path without '.', '..', drive or backslash components"));
    }
    Ok(value)
}

impl SourceBundle {
    /// The loader must enforce access roots and resolve symbolic links before reading.
    pub fn load(root: &str, text: &str, mut loader: impl FnMut(&str) -> Result<String>) -> Result<Self> {
        let root = relative_path(&format!("\"{root}\""), 1)?;
        let mut out = Self { root: root.clone(), files: BTreeMap::new(), expanded: String::new(), locations: Vec::new() };
        let mut bytes = 0usize;
        out.expand(&root, text, &mut loader, &mut BTreeSet::new(), &mut bytes, &mut 0)?;
        Ok(out)
    }

    fn expand(
        &mut self,
        name: &str,
        text: &str,
        loader: &mut impl FnMut(&str) -> Result<String>,
        active: &mut BTreeSet<String>,
        bytes: &mut usize,
        work: &mut usize,
    ) -> Result<()> {
        *work = work.saturating_add(1);
        super::check_work(*work)?;
        if *work > MAX_RECORDS {
            return Err(lex::error(1, "INCLUDE traversal work limit reached"));
        }
        if active.len() >= MAX_DEPTH || !active.insert(name.to_string()) {
            return Err(lex::error(1, format!("{name}: INCLUDE cycle or depth limit ({MAX_DEPTH})")));
        }
        let file_count = self.files.len();
        if let std::collections::btree_map::Entry::Vacant(entry) = self.files.entry(name.to_string()) {
            *bytes = bytes.saturating_add(text.len());
            if file_count >= MAX_FILES || *bytes > MAX_INPUT_BYTES {
                return Err(lex::error(1, "INCLUDE bundle exceeds file count or 64 MiB source budget"));
            }
            entry.insert(SourceDocument::parse(text).map_err(|e| match e {
                Error::Parse { line, message } => lex::error(line, format!("{name}: {message}")),
                other => other,
            })?);
        }
        // Prefix permits data-only INCLUDE fragments. The combined stream is validated later.
        let prefixed = format!("*SOURCE\n{text}");
        let blocks = lex::blocks(&prefixed, MAX_RECORDS).map_err(|e| match e {
            Error::Parse { line, message } => lex::error(line.saturating_sub(1).max(1), format!("{name}: {message}")),
            other => other,
        })?;
        let mut includes = BTreeMap::new();
        let physical_lines: Vec<&str> = text.lines().collect();
        for b in blocks.iter().filter(|b| b.keyword == "INCLUDE") {
            let line = b.line.saturating_sub(1);
            let location = |e| match e {
                Error::Parse { line, message } => lex::error(line.saturating_sub(1).max(1), format!("{name}: {message}")),
                other => other,
            };
            b.only(&["INPUT"]).map_err(location)?;
            if !b.rows.is_empty() {
                return Err(lex::error(line, format!("{name}: INCLUDE has no data lines")));
            }
            let target = relative_path(b.required("INPUT").map_err(location)?, line)?;
            // Keyword continuations may contain comments and blank lines.
            let mut end = line;
            let mut continued = false;
            for (index, physical) in physical_lines.iter().enumerate().skip(line.saturating_sub(1)) {
                let trimmed = physical.trim();
                end = index + 1;
                if trimmed.is_empty() || trimmed.starts_with("**") {
                    continue;
                }
                continued = trimmed.ends_with(',');
                if !continued {
                    break;
                }
            }
            if continued {
                return Err(lex::error(line, "unfinished INCLUDE continuation"));
            }
            includes.insert(line, (end, target));
        }
        let mut skip_until = 0;
        for (index, physical) in text.split_inclusive('\n').enumerate() {
            super::check_work(index)?;
            let line = index + 1;
            if line <= skip_until {
                continue;
            }
            if let Some((end, target)) = includes.get(&line) {
                skip_until = *end;
                let content = match self.files.get(target) {
                    Some(doc) => doc.original().to_string(),
                    None => loader(target).map_err(|e| lex::error(line, format!("{name}: INCLUDE {target}: {e}")))?,
                };
                self.expand(target, &content, loader, active, bytes, work).map_err(|e| lex::error(line, format!("{name}: INCLUDE {target}: {e}")))?;
            } else {
                if self.expanded.len().saturating_add(physical.len()).saturating_add(1) > MAX_INPUT_BYTES || self.locations.len() >= MAX_RECORDS {
                    return Err(lex::error(line, format!("{name}: expanded INCLUDE stream limit reached")));
                }
                self.expanded.push_str(physical);
                if !physical.ends_with('\n') {
                    self.expanded.push('\n');
                }
                self.locations.push(Location { file: name.to_string(), line });
            }
        }
        active.remove(name);
        Ok(())
    }

    pub fn location(&self, line: usize) -> Option<&Location> {
        self.locations.get(line.saturating_sub(1))
    }

    pub fn source_error(&self, error: Error) -> Error {
        match error {
            Error::Parse { line, message } => match self.location(line) {
                Some(at) => lex::error(at.line, format!("{}: {message}", at.file)),
                None => Error::Parse { line, message },
            },
            other => other,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn include_fragments_keep_bytes_and_use_root_relative_paths() {
        let root = "** root\r\n*NODE\r\n*INCLUDE,\r\n INPUT=\"sub/nodes.inc\"\r\n*STEP\r\n*STATIC\r\n1,1\r\n*END STEP";
        let nodes = "1,0,0,0\n*INCLUDE,INPUT=other.inc\n";
        let bundle = SourceBundle::load("model.inp", root, |name| match name {
            "sub/nodes.inc" => Ok(nodes.into()),
            "other.inc" => Ok("2,1,0,0".into()),
            _ => Err(lex::error(1, "unexpected path")),
        })
        .unwrap();
        assert_eq!(bundle.files["model.inp"].original(), root);
        assert_eq!(bundle.files["sub/nodes.inc"].original(), nodes);
        assert_eq!(bundle.files["other.inc"].original(), "2,1,0,0");
        let read = crate::abaqus::read(&bundle.expanded).unwrap();
        assert_eq!(read.model.nodes().len(), 2);
        assert!(read.model.steps().is_empty());
        assert!(bundle.locations.iter().any(|at| at.file == "other.inc" && at.line == 1));
    }
    #[test]
    fn refuses_cycle_escape_and_encrypted_include() {
        for input in ["*INCLUDE,INPUT=model.inp", "*INCLUDE,INPUT=../x", "*INCLUDE,INPUT=/x", "*INCLUDE,INPUT=C:\\x", "*INCLUDE,INPUT=x,PASSWORD=p"] {
            assert!(SourceBundle::load("model.inp", input, |_| Ok(input.to_string())).is_err());
        }
    }
    #[test]
    fn reports_error_in_dependency_at_its_own_line() {
        let bundle = SourceBundle::load("root.inp", "** comment\n*INCLUDE,INPUT=bad.inc\n", |_| Ok("*NODE\n1,NaN,0,0\n".into())).unwrap();
        let error = bundle.source_error(crate::abaqus::read(&bundle.expanded).unwrap_err()).to_string();
        assert!(error.contains("line 2") && error.contains("bad.inc"), "{error}");
    }
}
