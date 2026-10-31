//! Bounded option blocks from the public Abaqus input syntax reference.
use crate::{Error, Result};
use std::collections::BTreeMap;

pub(super) fn error(line: usize, message: impl Into<String>) -> Error {
    Error::Parse { line, message: message.into() }
}

#[derive(Debug)]
pub(super) struct Row<'a> {
    pub line: usize,
    pub text: &'a str,
}

#[derive(Debug)]
pub(super) struct Block<'a> {
    pub line: usize,
    pub keyword: String,
    pub params: BTreeMap<String, Option<String>>,
    pub rows: Vec<Row<'a>>,
}

pub(super) fn canonical(s: &str) -> String {
    s.chars().filter(|c| !c.is_ascii_whitespace()).map(|c| c.to_ascii_uppercase()).collect()
}

/// CSV tokens, retaining empty fields and not splitting inside a quoted name.
pub(super) fn fields(text: &str, line: usize) -> Result<Vec<String>> {
    let mut out = Vec::new();
    let (mut start, mut quoted) = (0, false);
    for (at, c) in text.char_indices() {
        if c == '"' {
            quoted = !quoted;
        }
        if c == ',' && !quoted {
            out.push(text.get(start..at).unwrap_or_default().trim().to_string());
            start = at + 1;
        }
    }
    if quoted {
        return Err(error(line, "unterminated quoted name"));
    }
    out.push(text.get(start..).unwrap_or_default().trim().to_string());
    Ok(out)
}

pub(super) fn unquote(value: &str, line: usize) -> Result<String> {
    let value = value.trim();
    if value.starts_with('"') && value.ends_with('"') && value.len() >= 2 {
        let inner = value.get(1..value.len() - 1).unwrap_or_default();
        if inner.contains('"') {
            return Err(error(line, "a name contains an embedded quote"));
        }
        Ok(inner.to_string())
    } else if value.contains('"') {
        Err(error(line, "quotes must surround the complete name"))
    } else {
        Ok(value.chars().filter(|c| !c.is_ascii_whitespace()).collect())
    }
}

impl Block<'_> {
    pub fn value(&self, key: &str) -> Option<&str> {
        self.params.get(key).and_then(|p| p.as_deref())
    }
    pub fn required(&self, key: &str) -> Result<&str> {
        self.value(key).filter(|s| !s.is_empty()).ok_or_else(|| error(self.line, format!("*{} needs {key}=", self.keyword)))
    }
    pub fn only(&self, allowed: &[&str]) -> Result<()> {
        if let Some(key) = self.params.keys().find(|k| !allowed.contains(&k.as_str())) {
            return Err(error(self.line, format!("*{} parameter {key} is not supported; no data imported from this block", self.keyword)));
        }
        for (key, value) in &self.params {
            let valid = if key == "GENERATE" { value.is_none() } else { value.as_ref().is_some_and(|s| !s.is_empty()) };
            if !valid {
                return Err(error(
                    self.line,
                    if key == "GENERATE" {
                        "GENERATE is a flag and does not take a value".to_string()
                    } else {
                        format!("*{} needs {key}=value, not a bare parameter", self.keyword)
                    },
                ));
            }
        }
        Ok(())
    }
}

fn keyword(text: &str, line: usize) -> Result<Block<'static>> {
    let tokens = fields(text, line)?;
    let name = tokens.first().map(|s| canonical(s)).unwrap_or_default();
    if name.is_empty() {
        return Err(error(line, "empty keyword"));
    }
    let mut params = BTreeMap::new();
    for token in tokens.iter().skip(1).filter(|s| !s.is_empty()) {
        let (key, value) = match token.split_once('=') {
            Some((k, v)) => (canonical(k), Some(v.trim().to_string())),
            None => (canonical(token), None),
        };
        if key.is_empty() || params.insert(key.clone(), value).is_some() {
            return Err(error(line, format!("empty or repeated keyword parameter {key}")));
        }
    }
    Ok(Block { line, keyword: name, params, rows: Vec::new() })
}

pub(super) fn blocks(text: &str, max_records: usize) -> Result<Vec<Block<'_>>> {
    let mut out: Vec<Block<'_>> = Vec::new();
    let mut pending: Option<(usize, String)> = None;
    let mut records = 0usize;
    super::check_work(0)?;
    for (index, line) in text.lines().enumerate() {
        super::check_work(index)?;
        let number = index + 1;
        if line.len() > super::super::MAX_LINE {
            return Err(error(number, "line exceeds 256 bytes"));
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("**") {
            continue;
        }
        records = records.saturating_add(1);
        if records > max_records {
            return Err(error(number, "input record limit reached"));
        }
        if !trimmed.is_ascii() {
            return Err(error(number, "keyword and data records must be ASCII"));
        }
        if let Some((start, mut value)) = pending.take() {
            if trimmed.starts_with('*') {
                return Err(error(start, "unfinished keyword continuation"));
            }
            value.push_str(trimmed);
            if value.len() > 4096 {
                return Err(error(start, "continued keyword exceeds 4096 bytes"));
            }
            if trimmed.ends_with(',') {
                pending = Some((start, value));
            } else {
                out.push(keyword(&value, start)?);
            }
        } else if let Some(name) = trimmed.strip_prefix('*') {
            if name.ends_with(',') {
                pending = Some((number, name.to_string()));
            } else {
                out.push(keyword(name, number)?);
            }
        } else if let Some(block) = out.last_mut() {
            block.rows.push(Row { line: number, text: trimmed });
        } else {
            return Err(error(number, "data before the first keyword"));
        }
    }
    if let Some((line, _)) = pending {
        return Err(error(line, "unfinished keyword continuation"));
    }
    Ok(out)
}
