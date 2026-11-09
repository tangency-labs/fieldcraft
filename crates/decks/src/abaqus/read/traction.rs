//! Fixed global traction on linear triangular faces.
//! Public syntax: https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dsload.htm
//! The native consistent face integral is exact for constant traction on a Tet4 face.
use super::{Block, Model, Result, Set, canonical, error, number, unsupported, value};
use fieldcraft_model::ElementKind;

pub(super) fn vector(b: &Block<'_>, fields: &[String], set: &Set, model: &Model, line: usize) -> Result<[f64; 3]> {
    if canonical(value(b, "FOLLOWER", "YES")?) != "NO" {
        return Err(unsupported(line, "DSLOAD TRVEC requires explicit FOLLOWER=NO"));
    }
    if canonical(value(b, "CONSTANTRESULTANT", "NO")?) != "YES" {
        return Err(unsupported(line, "DSLOAD TRVEC requires explicit CONSTANT RESULTANT=YES"));
    }
    for face in &set.faces {
        let element = model.element(face.element).ok_or_else(|| error(line, "missing TRVEC surface element"))?;
        if element.kind != ElementKind::Tet4 {
            return Err(unsupported(
                line,
                format!("DSLOAD TRVEC element {} ({:?}): only Tet4 triangular faces are mapped", face.element, element.kind),
            ));
        }
    }
    if fields.len() != 6 {
        return Err(error(line, "DSLOAD TRVEC needs surface, magnitude and three direction components"));
    }
    let magnitude = number(fields.get(2).map(String::as_str).unwrap_or_default(), line)?;
    let mut vector = [0.0; 3];
    for (axis, component) in vector.iter_mut().enumerate() {
        *component = number(fields.get(axis + 3).map(String::as_str).unwrap_or_default(), line)?;
    }
    let norm = vector.iter().fold(0.0_f64, |n, v| n.hypot(*v));
    if norm == 0.0 || !norm.is_finite() {
        return Err(error(line, "DSLOAD TRVEC direction must be finite and nonzero"));
    }
    Ok(vector.map(|v| magnitude * (v / norm)))
}
