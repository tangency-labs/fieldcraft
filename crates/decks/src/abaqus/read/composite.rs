//! Constant material frames and laminate rows from the public Abaqus 2026 keyword reference.
use super::{Block, Result, canonical, error, number, quantity, row_fields};
use fieldcraft_model::units::{Quantity, Unit};
use fieldcraft_model::{DisplaySystem, MaterialOrientation, Ply, SectionOrientation};

pub(super) struct Orientation {
    frame: MaterialOrientation,
    rotation_axis: u8,
}

impl Orientation {
    pub(super) fn section(&self, shell: bool) -> SectionOrientation {
        SectionOrientation { frame: self.frame, shell_reference_axis: if shell { self.rotation_axis % 3 + 1 } else { 1 } }
    }
}

pub(super) fn orientation(block: &Block<'_>, units: &DisplaySystem) -> Result<Orientation> {
    block.only(&["NAME", "SYSTEM", "DEFINITION"])?;
    if block.value("DEFINITION").is_some_and(|v| canonical(v) != "COORDINATES") {
        return Err(error(block.line, "*ORIENTATION supports constant DEFINITION=COORDINATES only"));
    }
    if block.rows.is_empty() || block.rows.len() > 2 {
        return Err(error(block.line, "*ORIENTATION needs one coordinate row and at most one rotation row"));
    }
    let row = block.rows.first().ok_or_else(|| error(block.line, "missing orientation coordinates"))?;
    let values = row_fields(row.text, row.line)?.iter().map(|v| number(v, row.line)).collect::<Result<Vec<_>>>()?;
    let at = |i| values.get(i).copied().unwrap_or_default();
    let frame = match canonical(block.value("SYSTEM").unwrap_or("RECTANGULAR")).as_str() {
        "RECTANGULAR" if values.len() == 6 || values.len() == 9 => MaterialOrientation::Rectangular {
            axis1: std::array::from_fn(|i| at(i) - at(i + 6)),
            axis2: std::array::from_fn(|i| at(i + 3) - at(i + 6)),
        },
        "CYLINDRICAL" if values.len() == 6 => {
            let mut origin = [0.0; 3];
            for (i, v) in origin.iter_mut().enumerate() {
                *v = quantity(units, Quantity::Length, at(i), row.line)?;
            }
            MaterialOrientation::Cylindrical { origin, axis: std::array::from_fn(|i| at(i + 3) - at(i)) }
        }
        _ => {
            return Err(error(
                row.line,
                "*ORIENTATION requires 6 (or 9 rectangular) coordinates; spherical, user and distributed frames are not imported",
            ));
        }
    };
    frame.validate().map_err(|e| error(row.line, e.to_string()))?;
    let mut rotation_axis = 1;
    if let Some(row) = block.rows.get(1) {
        let f = row_fields(row.text, row.line)?;
        if f.len() > 2 {
            return Err(error(row.line, "orientation rotation has only axis and angle"));
        }
        if let Some(v) = f.first().filter(|v| !v.is_empty()) {
            rotation_axis =
                v.parse::<u8>().ok().filter(|v| (1..=3).contains(v)).ok_or_else(|| error(row.line, "orientation rotation axis must be 1, 2 or 3"))?;
        }
        if let Some(v) = f.get(1).filter(|v| !v.is_empty())
            && number(v, row.line)? != 0.0
        {
            return Err(error(row.line, "nonzero additional *ORIENTATION rotation is not imported; supply the rotated constant frame explicitly"));
        }
    }
    Ok(Orientation { frame, rotation_axis })
}

pub(super) fn plies(block: &Block<'_>, units: &DisplaySystem) -> Result<Vec<Ply>> {
    if block.rows.is_empty() || block.rows.len() > fieldcraft_model::composite::MAX_PLIES {
        return Err(error(block.line, "a composite shell needs 1 to 256 ply rows"));
    }
    let angle_unit = Unit::parse("deg").map_err(|e| error(block.line, e.to_string()))?;
    let mut plies = Vec::with_capacity(block.rows.len());
    for row in &block.rows {
        let f = row_fields(row.text, row.line)?;
        if !(3..=4).contains(&f.len()) {
            return Err(error(row.line, "a ply needs thickness, point count, material and optional angle; ply names are not retained"));
        }
        let points = f
            .get(1)
            .filter(|v| !v.is_empty())
            .map(|v| v.parse::<u32>().map_err(|_| error(row.line, "invalid ply point count")))
            .transpose()?
            .unwrap_or(3);
        if points < 3 || points.is_multiple_of(2) {
            return Err(error(row.line, "a ply requires odd Simpson integration count at least 3 to preserve its elastic stiffness"));
        }
        let material = super::key(f.get(2).map(String::as_str).unwrap_or_default(), row.line)?;
        let angle = f.get(3).filter(|v| !v.is_empty()).map(|v| number(v, row.line)).transpose()?.unwrap_or(0.0);
        let ply = Ply {
            material,
            points,
            thickness: quantity(units, Quantity::Length, number(f.first().map(String::as_str).unwrap_or_default(), row.line)?, row.line)?,
            angle: Quantity::Angle.to_si(angle, &angle_unit).map_err(|e| error(row.line, e.to_string()))?,
        };
        ply.validate().map_err(|e| error(row.line, e.to_string()))?;
        plies.push(ply);
    }
    Ok(plies)
}

pub(super) fn offset(value: Option<&str>, line: usize) -> Result<f64> {
    match value.map(canonical).as_deref() {
        None => Ok(0.0),
        Some("SPOS") => Ok(0.5),
        Some("SNEG") => Ok(-0.5),
        Some(v) => number(v, line),
    }
}
