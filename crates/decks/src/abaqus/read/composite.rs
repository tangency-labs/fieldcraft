//! Constant material frames and laminate rows from the public Abaqus 2026 keyword reference.
use super::{Block, Result, canonical, error, number, quantity, row_fields};
use fieldcraft_model::units::{Quantity, Unit};
use fieldcraft_model::{DisplaySystem, MaterialOrientation, Ply, SectionOrientation};

/// Public *ELASTIC, TYPE=ORTHOTROPIC stiffness order. Invert the normal block
/// to the native reciprocal engineering constants; no plane-stress completion.
pub(super) fn orthotropic(v: &[f64], units: &DisplaySystem, line: usize) -> Result<fieldcraft_model::Elastic> {
    let scale = v.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(error(line, "*ELASTIC TYPE=ORTHOTROPIC requires a positive definite stiffness"));
    }
    let at = |i| v.get(i).copied().unwrap_or_default() / scale;
    let (a, b, c, d, e, f) = (at(0), at(1), at(2), at(3), at(4), at(5));
    let minor = a * c - b * b;
    let determinant = a * (c * f - e * e) - b * (b * f - d * e) + d * (b * e - c * d);
    if a <= 0.0 || minor <= 0.0 || determinant <= 0.0 || (6..9).any(|i| at(i) <= 0.0) {
        return Err(error(line, "*ELASTIC TYPE=ORTHOTROPIC requires a positive definite stiffness"));
    }
    let diagonal = [(c * f - e * e) / determinant, (a * f - d * d) / determinant, minor / determinant];
    let mut young = [0.0; 3];
    let mut shear = [0.0; 3];
    for i in 0..3 {
        young[i] = quantity(units, Quantity::Pressure, scale / diagonal[i], line)?;
        shear[i] = quantity(units, Quantity::Pressure, v.get(i + 6).copied().unwrap_or_default(), line)?;
    }
    let poisson =
        [(b * f - d * e) / (determinant * diagonal[0]), (c * d - b * e) / (determinant * diagonal[0]), (a * e - b * d) / (determinant * diagonal[1])];
    Ok(fieldcraft_model::Elastic::EngineeringConstants { young, poisson, shear })
}

pub(super) struct Orientation {
    frame: MaterialOrientation,
    rotation_axis: u8,
    rotation: f64,
    line: usize,
}

impl Orientation {
    pub(super) fn section(&self, shell: bool) -> Result<SectionOrientation> {
        let mut frame = self.frame;
        if self.rotation != 0.0 {
            if shell || !matches!(frame, MaterialOrientation::Rectangular { .. }) {
                return Err(error(
                    self.line,
                    "nonzero additional *ORIENTATION rotation requires a rectangular solid frame; shell projection and cylindrical twist are not mapped",
                ));
            }
            let axes = frame.axes_at([0.0; 3]).map_err(|e| error(self.line, e.to_string()))?;
            let axis = axes.get(usize::from(self.rotation_axis - 1)).copied().ok_or_else(|| error(self.line, "invalid rotation axis"))?;
            let (sin, cos) = self.rotation.sin_cos();
            let rotate = |v: [f64; 3]| {
                let dot = axis.iter().zip(v).map(|(a, b)| a * b).sum::<f64>();
                let cross = [axis[1] * v[2] - axis[2] * v[1], axis[2] * v[0] - axis[0] * v[2], axis[0] * v[1] - axis[1] * v[0]];
                std::array::from_fn(|i| cos * v[i] + sin * cross[i] + (1.0 - cos) * dot * axis[i])
            };
            frame = MaterialOrientation::Rectangular { axis1: rotate(axes[0]), axis2: rotate(axes[1]) };
        }
        Ok(SectionOrientation { frame, shell_reference_axis: if shell { self.rotation_axis % 3 + 1 } else { 1 } })
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
    let mut rotation = 0.0;
    if let Some(row) = block.rows.get(1) {
        let f = row_fields(row.text, row.line)?;
        if f.len() > 2 {
            return Err(error(row.line, "orientation rotation has only axis and angle"));
        }
        if let Some(v) = f.first().filter(|v| !v.is_empty()) {
            rotation_axis =
                v.parse::<u8>().ok().filter(|v| (1..=3).contains(v)).ok_or_else(|| error(row.line, "orientation rotation axis must be 1, 2 or 3"))?;
        }
        if let Some(v) = f.get(1).filter(|v| !v.is_empty()) {
            let unit = Unit::parse("deg").map_err(|e| error(row.line, e.to_string()))?;
            rotation = Quantity::Angle.to_si(number(v, row.line)?, &unit).map_err(|e| error(row.line, e.to_string()))?;
        }
    }
    Ok(Orientation { frame, rotation_axis, rotation, line: block.line })
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
