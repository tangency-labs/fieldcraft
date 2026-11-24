//! External result frames carry the visualization grid on which the solver sampled its fields.
//! Point and cell indices are local to each frame. Even coincident points may carry different
//! values across a discontinuity. No correspondence to the analysis mesh is implied.

use fieldcraft_model::units::{Quantity, Unit};
use std::collections::BTreeMap;

use crate::{Error, Result, invalid};

mod binary;
pub(crate) use binary::{decode, encode};

pub(crate) const PART: &str = "external.bin";
const MAGIC: &[u8; 8] = b"FCEXT001";

/// Aggregate storage budget for a result series, including frame geometry and field values.
pub const MAX_BYTES: usize = crate::MAX_GEOMETRY_PART_BYTES as usize;
/// File and decoded storage budget for any single visualization frame.
pub const MAX_FRAME_BYTES: usize = 256 << 20;
pub const MAX_POINTS: usize = 2_000_000;
pub const MAX_CELLS: usize = 2_000_000;
pub const MAX_FIELDS: usize = 256;
pub const MAX_COMPONENTS: usize = 9;
pub const MAX_FRAMES: usize = 4096;
pub(crate) const FRAME_BYTES: usize = 256;
pub(crate) const FIELD_BYTES: usize = 192;
pub(crate) const CELL_BYTES: usize = 64;

#[derive(Debug, Clone, PartialEq)]
pub struct Cell {
    pub vtk_type: u8,
    pub points: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Grid {
    /// Undeformed point coordinates in meters, in file order.
    pub points: Vec<[f64; 3]>,
    pub cells: Vec<Cell>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Field {
    pub components: u8,
    /// Tuples in VTK point or cell order, components adjacent.
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    pub time: f64,
    /// Source unit of arrays selected for displacement. The grid itself is always in meters.
    /// A missing unit denotes meters, including for older embedded frames.
    pub displacement_unit: Option<Unit>,
    pub grid: Grid,
    pub point_data: BTreeMap<String, Field>,
    pub cell_data: BTreeMap<String, Field>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Series {
    pub frames: Vec<Frame>,
}

pub(crate) fn valid_name(name: &str) -> bool {
    !name.trim().is_empty() && name.len() <= 256 && !name.chars().any(char::is_control)
}

/// Point count of a supported fixed-arity cell in the VTK XML unstructured grid format.
pub fn cell_arity(vtk_type: u8) -> Option<usize> {
    Some(match vtk_type {
        1 => 1,       // vertex
        3 => 2,       // line
        5 | 21 => 3,  // triangle, quadratic edge
        9 | 10 => 4,  // quad, tetrahedron
        14 => 5,      // pyramid
        13 | 22 => 6, // wedge, quadratic triangle
        12 | 23 => 8, // hexahedron, quadratic quad
        24 => 10,     // quadratic tetrahedron
        25 => 20,     // quadratic hexahedron
        _ => return None,
    })
}

impl Frame {
    /// Budget for live and encoded data. Counts include collection entries and each tuple's
    /// full width. The same conservative charges apply on 32-bit and 64-bit hosts.
    pub fn storage_bytes(&self) -> usize {
        let mut size = FRAME_BYTES + self.displacement_unit.as_ref().map_or(0, |u| u.symbol().len());
        size = size.saturating_add(self.grid.points.len().saturating_mul(24));
        for cell in &self.grid.cells {
            size = size.saturating_add(CELL_BYTES).saturating_add(cell.points.len().saturating_mul(4));
        }
        for (name, field) in self.point_data.iter().chain(&self.cell_data) {
            size = size.saturating_add(FIELD_BYTES).saturating_add(name.len()).saturating_add(field.values.len().saturating_mul(8));
        }
        size
    }

    pub fn check(&self) -> Result<()> {
        self.checked_storage(0).map(|_| ())
    }

    /// Check a frame and return the aggregate storage through this frame.
    fn checked_storage(&self, previous: usize) -> Result<usize> {
        if let Some(unit) = &self.displacement_unit {
            Quantity::Length.to_si(0.0, unit).map_err(|e| Error::Invalid(e.to_string()))?;
        }
        if !self.time.is_finite() {
            return invalid("external result frame time must be finite");
        }
        if !(1..=MAX_POINTS).contains(&self.grid.points.len()) || !(1..=MAX_CELLS).contains(&self.grid.cells.len()) {
            return invalid(format!("an external result grid needs 1 to {MAX_POINTS} points and 1 to {MAX_CELLS} cells"));
        }
        let frame_size = self.storage_bytes();
        if frame_size > MAX_FRAME_BYTES {
            return Err(Error::Limit(format!("external result frame exceeds the {} MiB storage limit", MAX_FRAME_BYTES >> 20)));
        }
        let size = previous.saturating_add(frame_size);
        if size > MAX_BYTES {
            return Err(storage_limit());
        }
        if self.grid.points.iter().flatten().any(|v| !v.is_finite() || v.abs() > fieldcraft_model::MAX_COORD) {
            return invalid("external result grid coordinates must be finite and within the model coordinate limit");
        }
        for cell in &self.grid.cells {
            if cell_arity(cell.vtk_type) != Some(cell.points.len()) {
                return invalid(format!("unsupported VTK cell type {} or wrong number of points for that type", cell.vtk_type));
            }
            if cell.points.iter().any(|&index| index as usize >= self.grid.points.len()) {
                return invalid("external result cell refers to a point outside its grid");
            }
        }
        let fields = self.point_data.len().saturating_add(self.cell_data.len());
        if !(1..=MAX_FIELDS).contains(&fields) {
            return invalid(format!("an external result frame needs 1 to {MAX_FIELDS} fields"));
        }
        for (data, tuples) in [(&self.point_data, self.grid.points.len()), (&self.cell_data, self.grid.cells.len())] {
            for (name, field) in data {
                if !valid_name(name) {
                    return invalid("external result field names must be nonempty, at most 256 bytes, and free of control characters");
                }
                let components = usize::from(field.components);
                if !(1..=MAX_COMPONENTS).contains(&components) || tuples.checked_mul(components) != Some(field.values.len()) {
                    return invalid(format!("external field `{name}` has an invalid component or tuple count"));
                }
                if !field.values.iter().all(|value| value.is_finite()) {
                    return invalid(format!("external field `{name}` contains a nonfinite value"));
                }
            }
        }
        Ok(size)
    }
}

/// Stored display settings use the same flat scalar parameters as the display command.
pub(crate) fn check_view(view: Option<&serde_json::Value>, has_frames: bool) -> Result<()> {
    let Some(view) = view else { return Ok(()) };
    let Some(params) = view.as_object() else { return invalid("external_view must be an object") };
    if !has_frames {
        return invalid("external_view requires external result frames");
    }
    if params.len() > 16
        || params.iter().any(|(name, v)| name.len() > 64 || v.is_array() || v.is_object() || v.as_str().is_some_and(|s| s.len() > 256))
    {
        return invalid("external_view must contain bounded scalar display parameters");
    }
    Ok(())
}

pub(crate) fn charge_storage(remaining: &mut usize, cost: usize) -> Result<()> {
    *remaining = remaining.checked_sub(cost).ok_or_else(storage_limit)?;
    Ok(())
}

fn storage_limit() -> Error {
    Error::Limit(format!("external results exceed the storage budget ({} MiB per frame, {} MiB per series)", MAX_FRAME_BYTES >> 20, MAX_BYTES >> 20))
}

impl Series {
    fn check_timeline(&self) -> Result<()> {
        if !(1..=MAX_FRAMES).contains(&self.frames.len()) {
            return invalid(format!("external results must have 1 to {MAX_FRAMES} frames"));
        }
        if self.frames.iter().zip(self.frames.iter().skip(1)).any(|(a, b)| a.time >= b.time) {
            return invalid("external result times must increase strictly");
        }
        Ok(())
    }

    pub fn check(&self) -> Result<()> {
        self.check_timeline()?;
        self.frames.iter().try_fold(0usize, |size, frame| frame.checked_storage(size)).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(time: f64) -> Frame {
        Frame {
            displacement_unit: None,
            time,
            grid: Grid { points: vec![[0.0; 3], [1.0, 0.0, 0.0]], cells: vec![Cell { vtk_type: 3, points: vec![1, 0] }] },
            point_data: BTreeMap::from([("u".into(), Field { components: 2, values: vec![0.0, -0.0, 1.0, -0.25] })]),
            cell_data: BTreeMap::from([("u".into(), Field { components: 1, values: vec![3.25] })]),
        }
    }

    #[test]
    fn binary_round_trip_keeps_cell_order_and_frame_specific_grids() {
        let first = frame(0.0);
        let mut second = frame(2.5);
        second.grid.points[0] = [0.5, -0.5, 1.25];
        second.grid.cells[0].points.reverse();
        let original = Series { frames: vec![first, second] };
        let bytes = encode(&original).unwrap();
        let restored = decode(&bytes).unwrap();
        assert_eq!(restored, original);
        assert_eq!(encode(&restored).unwrap(), bytes);
        assert_eq!(restored.frames[0].point_data["u"].values[1].to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn documented_vertex_record_reads_and_writes_exactly() {
        let bytes = [
            b"FCEXT001".as_slice(),
            &1u32.to_le_bytes(),   // frames
            &0.5f64.to_le_bytes(), // time
            &1u32.to_le_bytes(),   // points
            &[0; 24],              // origin
            &1u32.to_le_bytes(),   // cells
            &[1],                  // VTK_VERTEX
            &0u32.to_le_bytes(),   // point index
            &1u16.to_le_bytes(),   // point fields
            &1u16.to_le_bytes(),   // name length
            b"s",
            &[1], // components
            &1.5f64.to_le_bytes(),
            &0u16.to_le_bytes(), // cell fields
        ]
        .concat();
        let restored = decode(&bytes).unwrap();
        assert_eq!(restored.frames[0].time, 0.5);
        assert_eq!(restored.frames[0].point_data["s"].values, [1.5]);
        assert_eq!(encode(&restored).unwrap(), bytes);
        let mut excessive = bytes;
        excessive[20..24].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(decode(&excessive).is_err());
    }

    #[test]
    fn invalid_frames_and_damaged_lengths_are_refused() {
        let valid = frame(1.0);
        let mut bad = valid.clone();
        bad.grid.cells[0].points[0] = 2;
        assert!(bad.check().is_err());
        bad = valid.clone();
        bad.cell_data.get_mut("u").unwrap().values.clear();
        assert!(bad.check().is_err());
        bad = valid.clone();
        bad.point_data.get_mut("u").unwrap().values[1] = f64::INFINITY;
        assert!(bad.check().is_err());
        assert!(Series { frames: vec![valid.clone(), valid.clone()] }.check().is_err());
        let bytes = encode(&Series { frames: vec![valid] }).unwrap();
        for end in 0..bytes.len() {
            assert!(decode(&bytes[..end]).is_err(), "{end}");
        }
        let mut hostile = bytes;
        hostile[8..16].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(decode(&hostile).is_err());
    }
}
