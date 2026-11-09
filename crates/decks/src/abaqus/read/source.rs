//! Original source storage, separate from the interpreted finite element model.
//! Keeping the bytes avoids reconstructing comments, spelling or line endings from tokens.
use super::{MAX_INPUT_BYTES, MAX_RECORDS, lex};
use crate::Result;

/// A keyword's original location. The text includes its data and intervening comments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceBlock {
    pub keyword: String,
    pub line: usize,
    start: usize,
    end: usize,
}

/// Immutable input text. Replaying this document does not export edits to a model and does
/// not imply that its keywords can be solved. Includes are retained, not opened here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceDocument {
    text: String,
    blocks: Vec<SourceBlock>,
}

impl SourceDocument {
    pub fn parse(text: &str) -> Result<Self> {
        if text.len() > MAX_INPUT_BYTES {
            return Err(lex::error(1, "input exceeds 64 MiB"));
        }
        let prefixed = format!("*SOURCE\n{text}");
        let parsed = lex::blocks(&prefixed, MAX_RECORDS).map_err(|e| match e {
            crate::Error::Parse { line, message } => lex::error(line.saturating_sub(1).max(1), message),
            other => other,
        })?;
        let mut offsets = Vec::new();
        let mut offset = 0;
        for line in text.split_inclusive('\n') {
            super::check_work(offsets.len())?;
            offsets.push(offset);
            offset += line.len();
        }
        let mut blocks: Vec<SourceBlock> = Vec::new();
        for block in parsed.into_iter().skip(1) {
            let line = block.line.saturating_sub(1);
            let start = offsets.get(line.saturating_sub(1)).copied().ok_or_else(|| lex::error(block.line, "source location is outside input"))?;
            if let Some(previous) = blocks.last_mut() {
                previous.end = start;
            }
            blocks.push(SourceBlock { keyword: block.keyword, line, start, end: text.len() });
        }
        Ok(Self { text: text.to_string(), blocks })
    }

    /// Exact source, including preamble, blank lines, CRLF and absence of a final newline.
    pub fn original(&self) -> &str {
        &self.text
    }

    pub fn blocks(&self) -> &[SourceBlock] {
        &self.blocks
    }

    /// Indexed access prevents using another document's offsets with this document.
    pub fn block_text(&self, index: usize) -> Option<&str> {
        let block = self.blocks.get(index)?;
        self.text.get(block.start..block.end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_unknown_blocks_and_physical_lines() {
        let text = "** original\r\n\r\n*Future Option,\r\n flag=YES\r\n1, 2\r\n** note\r\n*NODE\r\n1,0,0,0";
        let document = SourceDocument::parse(text).unwrap();
        assert_eq!(document.original(), text);
        assert_eq!(document.blocks()[0].line, 3);
        assert_eq!(document.blocks()[0].keyword, "FUTUREOPTION");
        assert_eq!(document.block_text(0), Some("*Future Option,\r\n flag=YES\r\n1, 2\r\n** note\r\n"));
        assert_eq!(SourceDocument::parse(document.original()).unwrap(), document);
        assert_eq!(document.block_text(2), None);
    }
}
