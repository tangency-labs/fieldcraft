//! Clean-room ASCII Universal File Format reader. Public UC-SDRL dataset definitions:
//! <https://www.ceas3.uc.edu/sdrluff/files/15.asc> (1987-08-30),
//! <https://www.ceas3.uc.edu/sdrluff/files/55.asc> (1997-03-07),
//! <https://www.ceas3.uc.edu/sdrluff/files/58.asc> (1999-06-30).
//! Fixed-width numeric records are decoded without depending on vendor implementations.
//! Values retain file units. Coordinate systems, axis metadata and complex values are never
//! discarded or guessed. Conversion into solver observations is a separate checked boundary.
use crate::{Error, Result};
use serde::Serialize;
use std::collections::BTreeSet;

pub const MAX_BYTES: u64 = 16 << 20;
const MAX_LINES: usize = 500_000;
const MAX_VALUES: usize = 1_000_000;

#[derive(Debug, Serialize)]
pub struct Node {
    pub label: u32,
    pub definition_system: i32,
    pub displacement_system: i32,
    pub xyz: [f64; 3],
}
#[derive(Debug, Serialize)]
pub struct Shape {
    pub node: u32,
    /// Pairs [real, imaginary], in global translation then rotation order.
    pub components: Vec<[f64; 2]>,
}
#[derive(Debug, Serialize)]
pub struct Mode {
    pub source_line: usize,
    pub label: u32,
    pub load_case: i32,
    pub frequency_hz: f64,
    pub modal_mass: f64,
    pub viscous_damping_ratio: f64,
    pub hysteretic_damping_ratio: f64,
    pub specific_data_type: i32,
    pub complex: bool,
    pub shapes: Vec<Shape>,
}
#[derive(Debug, Serialize)]
pub struct Axis {
    pub data_type: i32,
    pub exponents: [i32; 3],
    pub label: String,
    pub unit: String,
}
#[derive(Debug, Serialize)]
pub struct Frf {
    pub response_entity: String,
    pub response_node: u32,
    pub response_direction: i32,
    pub reference_entity: String,
    pub reference_node: u32,
    pub reference_direction: i32,
    pub axes: Vec<Axis>,
    /// [abscissa, real, imaginary] in file units; no assumed Fourier convention.
    pub samples: Vec<[f64; 3]>,
}
#[derive(Debug, Default, Serialize)]
pub struct Imported {
    pub nodes: Vec<Node>,
    pub modes: Vec<Mode>,
    pub frfs: Vec<Frf>,
    pub skipped_datasets: Vec<i32>,
}
fn invalid(message: impl Into<String>) -> Error {
    Error::Invalid(format!("UFF: {}", message.into()))
}
fn field(line: &str, start: usize, width: usize) -> Result<&str> {
    line.get(start..start.saturating_add(width)).map(str::trim).ok_or_else(|| invalid("truncated fixed-width record"))
}
fn int(s: &str) -> Result<i32> {
    s.parse().map_err(|_| invalid(format!("invalid integer `{s}`")))
}
fn label(s: &str) -> Result<u32> {
    let n = int(s)?;
    if n <= 0 {
        return Err(invalid("node/mode label must be positive"));
    }
    Ok(n as u32)
}
fn real(s: &str) -> Result<f64> {
    let x: f64 = s.replace(['D', 'd'], "E").parse().map_err(|_| invalid(format!("invalid real `{s}`")))?;
    if !x.is_finite() {
        return Err(invalid("nonfinite value"));
    }
    Ok(x)
}
fn integers(line: &str, count: usize) -> Result<Vec<i32>> {
    (0..count).map(|i| int(field(line, i * 10, 10)?)).collect()
}
fn reals(line: &str, width: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    for chunk in line.as_bytes().chunks(width) {
        let text = std::str::from_utf8(chunk).map_err(|_| invalid("non-ASCII numeric record"))?.trim();
        if !text.is_empty() {
            values.push(real(text)?);
        }
    }
    Ok(values)
}
fn row<'a>(rows: &[&'a str], i: usize) -> Result<&'a str> {
    rows.get(i).copied().ok_or_else(|| invalid("missing record"))
}
fn value(values: &[f64], i: usize) -> Result<f64> {
    values.get(i).copied().ok_or_else(|| invalid("missing value"))
}

/// Unknown datasets are skipped and listed. 55 supports normal-mode global 3/6-DOF vectors
/// (real or complex); 58 supports single-point complex FRFs in single/double precision,
/// evenly/unevenly spaced abscissa. Unsupported variants fail by name, including 58b.
pub fn read(text: &str) -> Result<Imported> {
    if text.len() as u64 > MAX_BYTES || !text.is_ascii() {
        return Err(invalid("expected ASCII UFF of at most 16 MiB"));
    }
    let mut lines = Vec::new();
    for l in text.lines() {
        if l.len() > 4096 || lines.len() >= MAX_LINES {
            return Err(invalid("line limit exceeded"));
        }
        lines.push(l);
    }
    let mut out = Imported::default();
    let mut cursor = 0;
    let mut labels = BTreeSet::new();
    let mut modes = BTreeSet::new();
    while cursor < lines.len() {
        if row(&lines, cursor)?.trim().is_empty() {
            cursor += 1;
            continue;
        }
        if row(&lines, cursor)?.trim() != "-1" {
            return Err(invalid(format!("expected delimiter at line {}", cursor + 1)));
        }
        cursor += 1;
        let id_text = row(&lines, cursor)?.trim();
        if id_text.to_ascii_lowercase().starts_with("58b") {
            return Err(invalid("binary dataset 58b is not supported"));
        }
        let id = int(id_text)?;
        cursor += 1;
        let begin = cursor;
        while cursor < lines.len() && row(&lines, cursor)?.trim() != "-1" {
            cursor += 1;
        }
        if cursor == lines.len() {
            return Err(invalid(format!("unterminated dataset {id}")));
        }
        let rows = lines.get(begin..cursor).ok_or_else(|| invalid("invalid block"))?;
        match id {
            15 => {
                for l in rows {
                    let label = label(field(l, 0, 10)?)?;
                    if !labels.insert(label) {
                        return Err(invalid("duplicate node label"));
                    }
                    let definition_system = int(field(l, 10, 10)?)?;
                    let displacement_system = int(field(l, 20, 10)?)?;
                    let _color = int(field(l, 30, 10)?)?;
                    let xyz = [real(field(l, 40, 13)?)?, real(field(l, 53, 13)?)?, real(field(l, 66, 13)?)?];
                    out.nodes.push(Node { label, definition_system, displacement_system, xyz });
                }
            }
            55 => {
                let m = mode(rows, begin + 1)?;
                if !modes.insert((m.load_case, m.label)) {
                    return Err(invalid("duplicate modal label within load case"));
                }
                out.modes.push(m);
            }
            58 => out.frfs.push(frf(rows)?),
            _ => out.skipped_datasets.push(id),
        }
        cursor += 1;
    }
    if out.nodes.is_empty() && out.modes.is_empty() && out.frfs.is_empty() && out.skipped_datasets.is_empty() {
        return Err(invalid("empty file"));
    }
    Ok(out)
}
fn mode(rows: &[&str], source_line: usize) -> Result<Mode> {
    let def = integers(row(rows, 5)?, 6)?;
    let [model, analysis, characteristic, specific, data_type, ndv] = def.as_slice() else {
        return Err(invalid("55 data definition"));
    };
    if *model != 1 || *analysis != 2 || !matches!((*characteristic, *ndv), (2, 3) | (3, 6)) || !matches!(*data_type, 2 | 5) {
        return Err(invalid("55 supports structural normal modes with global 3/6-DOF real/complex vectors only"));
    }
    let ints = integers(row(rows, 6)?, 4)?;
    let [2, 4, load_case, mode_label] = ints.as_slice() else {
        return Err(invalid("55 normal-mode parameter counts must be 2,4"));
    };
    if *mode_label <= 0 {
        return Err(invalid("55 mode label must be positive"));
    }
    let parameters = reals(row(rows, 7)?, 13)?;
    if parameters.len() != 4 {
        return Err(invalid("55 needs four real modal parameters"));
    }
    let frequency_hz = value(&parameters, 0)?;
    if frequency_hz < 0.0 {
        return Err(invalid("55 negative frequency"));
    }
    let complex = *data_type == 5;
    let count = *ndv as usize * if complex { 2 } else { 1 };
    let mut cursor = 8;
    let mut seen = BTreeSet::new();
    let mut shapes = Vec::new();
    while cursor < rows.len() {
        let node = label(field(row(rows, cursor)?, 0, 10)?)?;
        if !seen.insert(node) {
            return Err(invalid("55 duplicate shape node"));
        }
        cursor += 1;
        let mut v = Vec::new();
        while v.len() < count {
            v.extend(reals(row(rows, cursor)?, 13)?);
            cursor += 1;
        }
        if v.len() != count {
            return Err(invalid("55 wrong component count"));
        }
        let mut components = Vec::new();
        for i in 0..*ndv as usize {
            components.push(if complex { [value(&v, 2 * i)?, value(&v, 2 * i + 1)?] } else { [value(&v, i)?, 0.0] });
        }
        shapes.push(Shape { node, components });
    }
    if shapes.is_empty() {
        return Err(invalid("55 has no shape data"));
    }
    Ok(Mode {
        source_line,
        label: *mode_label as u32,
        load_case: *load_case,
        frequency_hz,
        modal_mass: value(&parameters, 1)?,
        viscous_damping_ratio: value(&parameters, 2)?,
        hysteretic_damping_ratio: value(&parameters, 3)?,
        specific_data_type: *specific,
        complex,
        shapes,
    })
}
fn axis(line: &str) -> Result<Axis> {
    Ok(Axis {
        data_type: int(field(line, 0, 10)?)?,
        exponents: [int(field(line, 10, 5)?)?, int(field(line, 15, 5)?)?, int(field(line, 20, 5)?)?],
        label: field(line, 26, 20)?.into(),
        unit: field(line, 47, 20)?.into(),
    })
}
fn frf(rows: &[&str]) -> Result<Frf> {
    let dofs = row(rows, 5)?;
    if int(field(dofs, 0, 5)?)? != 4 || int(field(dofs, 20, 10)?)? != 0 {
        return Err(invalid("58 supports single-point frequency response functions only"));
    }
    let response_direction = int(field(dofs, 51, 4)?)?;
    let reference_direction = int(field(dofs, 76, 4)?)?;
    if response_direction.unsigned_abs() > 6 || reference_direction.unsigned_abs() > 6 {
        return Err(invalid("58 direction must be -6..6"));
    }
    let form = row(rows, 6)?;
    let dtype = int(field(form, 0, 10)?)?;
    let n = int(field(form, 10, 10)?)?;
    let spacing = int(field(form, 20, 10)?)?;
    if !matches!(dtype, 5 | 6) || !matches!(spacing, 0 | 1) || n <= 0 || n as usize > MAX_VALUES / 3 {
        return Err(invalid("58 unsupported data type/spacing or invalid sample count"));
    }
    let start = real(field(form, 30, 13)?)?;
    let increment = real(field(form, 43, 13)?)?;
    if spacing == 1 && (increment <= 0.0 || start < 0.0) {
        return Err(invalid("58 even spacing must be positive with nonnegative start"));
    }
    let axes = (7..11).map(|i| axis(row(rows, i)?)).collect::<Result<Vec<_>>>()?;
    let mut values = Vec::new();
    for l in rows.get(11..).ok_or_else(|| invalid("58 missing samples"))? {
        if dtype == 6 && spacing == 0 {
            // Complex double uneven: exactly E13.5,2E20.12 per sample record.
            values.extend([real(field(l, 0, 13)?)?, real(field(l, 13, 20)?)?, real(field(l, 33, 20)?)?]);
            if l.get(53..).is_some_and(|s| !s.trim().is_empty()) {
                return Err(invalid("58 excess values"));
            }
        } else {
            values.extend(reals(l, if dtype == 6 { 20 } else { 13 })?);
        }
        if values.len() > MAX_VALUES {
            return Err(invalid("58 data cap exceeded"));
        }
    }
    let stride = if spacing == 0 { 3 } else { 2 };
    if values.len() != n as usize * stride {
        return Err(invalid("58 sample count mismatch"));
    }
    let mut samples = Vec::new();
    for (i, c) in values.chunks(stride).enumerate() {
        let sample =
            if spacing == 0 { [value(c, 0)?, value(c, 1)?, value(c, 2)?] } else { [start + i as f64 * increment, value(c, 0)?, value(c, 1)?] };
        if !sample[0].is_finite() || sample[0] < 0.0 || samples.last().is_some_and(|old: &[f64; 3]| old[0] >= sample[0]) {
            return Err(invalid("58 abscissa must be finite, nonnegative, strictly increasing"));
        }
        samples.push(sample);
    }
    Ok(Frf {
        response_entity: field(dofs, 31, 10)?.into(),
        response_node: label(field(dofs, 41, 10)?)?,
        response_direction,
        reference_entity: field(dofs, 56, 10)?.into(),
        reference_node: label(field(dofs, 66, 10)?)?,
        reference_direction,
        axes,
        samples,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn block(id: i32, rows: &str) -> String {
        format!("    -1\n{id:6}\n{rows}\n    -1\n")
    }
    fn nums(v: &[f64], width: usize) -> String {
        v.iter().map(|x| format!("{x:>width$.5E}")).collect()
    }
    fn node() -> String {
        format!("{:10}{:10}{:10}{:10}{}", 4, 0, 0, 1, nums(&[1.0, -2.0, 3.0], 13))
    }
    fn modal() -> String {
        let mut s = "NONE\n".repeat(5);
        s += &format!("{:10}{:10}{:10}{:10}{:10}{:10}\n", 1, 2, 2, 8, 5, 3);
        s += &format!("{:10}{:10}{:10}{:10}\n", 2, 4, 1, 9);
        s += &nums(&[63.039, 1.0, 0.02, 0.0], 13);
        s += &format!("\n{:10}\n", 4);
        s += &nums(&[1.0, -2.0, 3.0, -4.0, 5.0, -6.0], 13);
        s
    }
    fn response(double: bool, uneven: bool) -> String {
        let mut s = "NONE\n".repeat(5);
        s += &format!("{:5}{:10}{:5}{:10} {:10}{:10}{:4} {:10}{:10}{:4}\n", 4, 1, 0, 0, "response", 4, -3, "force", 7, 3);
        s += &format!("{:10}{:10}{:10}{}\n", if double { 6 } else { 5 }, 2, if uneven { 0 } else { 1 }, nums(&[1.0, 2.0, 0.0], 13));
        for (t, l, u) in [(18, "frequency", "Hz"), (12, "acceleration", "m/s^2"), (13, "force", "N"), (0, "NONE", "NONE")] {
            s += &format!("{t:10}{:5}{:5}{:5} {l:20} {u:20}\n", 0, 0, 0);
        }
        if uneven && double {
            s += &format!("{}{}\n{}{}", nums(&[1.0], 13), nums(&[2.0, -3.0], 20), nums(&[3.0], 13), nums(&[4.0, -5.0], 20));
        } else {
            s += &nums(if uneven { &[1.0, 2.0, -3.0, 3.0, 4.0, -5.0] } else { &[2.0, -3.0, 4.0, -5.0] }, if double { 20 } else { 13 });
        }
        s
    }
    #[test]
    fn fixed_width_nodes_complex_shapes_and_skips() {
        let text = block(15, &node()) + &block(55, &modal()) + &block(82, "arbitrary connectivity");
        let data = read(&text).unwrap();
        assert_eq!(data.nodes[0].xyz, [1.0, -2.0, 3.0]);
        assert_eq!(data.modes[0].shapes[0].components[1], [3.0, -4.0]);
        assert_eq!(data.modes[0].label, 9);
        assert_eq!(data.skipped_datasets, vec![82]);
    }
    #[test]
    fn all_complex_frf_layouts_retain_reference_and_phase() {
        for double in [false, true] {
            for uneven in [false, true] {
                let data = read(&block(58, &response(double, uneven))).unwrap();
                let f = &data.frfs[0];
                assert_eq!(f.samples, vec![[1.0, 2.0, -3.0], [3.0, 4.0, -5.0]]);
                assert_eq!(f.response_direction, -3);
                assert_eq!(f.reference_node, 7);
                assert_eq!(f.axes[1].unit, "m/s^2");
            }
        }
    }
    #[test]
    fn damaged_files_and_hostile_directions_fail() {
        for text in [
            String::new(),
            "    -1\n58b\n".into(),
            block(15, "short"),
            block(55, "NONE"),
            block(58, "NONE"),
            block(15, &node()) + &block(15, &node()),
        ] {
            assert!(read(&text).is_err());
        }
        let frf = response(false, false).replace("  -3", "xxxx");
        assert!(read(&block(58, &frf)).is_err());
        let complete = block(55, &modal());
        assert!(read(complete.trim_end_matches("    -1\n")).is_err());
    }
}
