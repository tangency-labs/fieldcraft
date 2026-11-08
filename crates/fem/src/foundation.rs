//! Consistent Winkler surface foundation on a fixed reference geometry.
//!
//! Virtual work gives K_ab = ∫ N_a N_b [k_n nnᵀ + k_t(I-nnᵀ)] dA
//! (Hetényi, *Beams on Elastic Foundation*, 1946; Zienkiewicz, Taylor & Zhu,
//! *The Finite Element Method*, 7th ed., 2013, virtual-work formulation).
//! The same kernel integrates viscous coefficients. Values are per unit area;
//! the returned matrix acts on three translations per node, never rotations.
//! Shape functions and quadrature are shared with surface tractions. Flat tri6
//! products have degree four; quad8 products have degree at most four per axis.

use crate::{Error, Result};

/// Row-major 3n × 3n surface matrix, with isotropic in-plane tangential support.
/// Zero coefficients are allowed, including for parameter derivatives at zero.
pub fn matrix(xyz: &[[f64; 3]], normal: f64, tangential: f64) -> Result<Vec<f64>> {
    if [normal, tangential].iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(Error::Invalid("foundation coefficients must be finite and non-negative".into()));
    }
    // Validate the node count before allocating or multiplying input-derived sizes.
    fieldcraft_model::nodal::face_points(xyz.len()).map_err(|e| Error::Invalid(e.to_string()))?;
    let width = 3 * xyz.len();
    let mut matrix = vec![0.0; width * width];
    let mut invalid = false;
    fieldcraft_model::nodal::integrate_face(xyz, |shape, weighted_normal| {
        let area = weighted_normal[0].hypot(weighted_normal[1]).hypot(weighted_normal[2]);
        if !area.is_finite() || area <= 0.0 {
            invalid = true;
            return;
        }
        let n = weighted_normal.map(|v| v / area);
        for (a, na) in shape.iter().enumerate() {
            for (b, nb) in shape.iter().enumerate() {
                for r in 0..3 {
                    for c in 0..3 {
                        let projection = (normal - tangential) * n[r] * n[c] + if r == c { tangential } else { 0.0 };
                        if let Some(entry) = matrix.get_mut((3 * a + r) * width + 3 * b + c) {
                            *entry += na * nb * area * projection;
                        }
                    }
                }
            }
        }
    })
    .map_err(|e| Error::Invalid(e.to_string()))?;
    if invalid || matrix.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("foundation face is degenerate or its matrix overflows".into()));
    }
    Ok(matrix)
}
