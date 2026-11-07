//! Plane-stress first-ply failure measures, without stiffness degradation.
//!
//! References: S. W. Tsai & E. M. Wu, *A general theory of strength for anisotropic
//! materials*, J. Composite Materials 5:58–80 (1971), doi:10.1177/002199837100500106;
//! Z. Hashin, *Failure Criteria for Unidirectional Fiber Composites*, J. Applied Mechanics
//! 47:329–334 (1980), doi:10.1115/1.3153664. Public equations are also given in the Abaqus
//! theory/material manuals, "Plane stress orthotropic failure measures" and "Damage
//! initiation for fiber-reinforced composites" (documentation only, no vendor code or decks).
//!
//! Tensile stress is positive, strengths are positive magnitudes. Tsai–Hill uses signed
//! σ11σ22 (not Azzi–Tsai–Hill's absolute product), with X/Y selected by each normal stress
//! sign. All quadratic measures are the original polynomials, not their square roots or
//! inverse reserve factors. Tsai–Wu may be negative. Hashin is the 1980 plane-stress form,
//! α=1, with the squared fibre-compression convention; only the active sign branch is stored.

use fieldcraft_model::PlyStrength;

use crate::{Error, Result};

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct HashinIndices {
    pub fiber_tension: Option<f64>,
    pub fiber_compression: Option<f64>,
    pub matrix_tension: Option<f64>,
    pub matrix_compression: Option<f64>,
    pub maximum: f64,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FailureIndices {
    pub maximum_stress: f64,
    pub tsai_hill: f64,
    pub tsai_wu: f64,
    pub hashin: HashinIndices,
}

/// [σ11, σ22, τ12] in the ply's material axes, in Pa. No interlaminar or 3D failure claim.
pub fn plane_stress_indices(stress: [f64; 3], strengths: &PlyStrength) -> Result<FailureIndices> {
    strengths.validate().map_err(|e| Error::Invalid(e.to_string()))?;
    if stress.iter().any(|s| !s.is_finite()) {
        return Err(Error::Invalid("failure indices require finite material-axis stresses".into()));
    }
    let [s1, s2, tau] = stress;
    let [xt, yt] = strengths.tension;
    let [xc, yc] = strengths.compression;
    let x = if s1 >= 0.0 { xt } else { xc };
    let y = if s2 >= 0.0 { yt } else { yc };
    let shear = (tau / strengths.shear12).powi(2);
    let maximum_stress = (s1 / x).abs().max((s2 / y).abs()).max((tau / strengths.shear12).abs());
    let tsai_hill = (s1 / x).powi(2) - (s1 / x) * (s2 / x) + (s2 / y).powi(2) + shear;
    // Normalize before multiplying: products of dimensional strengths can overflow even
    // while the nondimensional index is representable.
    let a = s1 / xt.sqrt() / xc.sqrt();
    let b = s2 / yt.sqrt() / yc.sqrt();
    // Factoring each uniaxial polynomial preserves FI=1 at either strength even when
    // tensile and compressive strengths have very different magnitudes.
    let normal = |s: f64, tension: f64, compression: f64| {
        let (u, v) = (s / tension, s / compression);
        if s >= 0.0 { u + v * (u - 1.0) } else { -v + u * (1.0 + v) }
    };
    let tsai_wu = normal(s1, xt, xc) + normal(s2, yt, yc) + 2.0 * strengths.tsai_wu_interaction * a * b + shear;
    let fiber = if s1 >= 0.0 { (s1 / xt).powi(2) + shear } else { (s1 / xc).powi(2) };
    let matrix = if s2 >= 0.0 {
        (s2 / yt).powi(2) + shear
    } else {
        let r = s2 / yc;
        // σ²/(4S23²) + (Yc²/(4S23²) - 1) σ/Yc, factored at σ=-Yc.
        let coupled = (s2 / strengths.shear23 / 2.0) * (1.0 + r);
        coupled * (yc / strengths.shear23 / 2.0) - r + shear
    };
    if [maximum_stress, tsai_hill, tsai_wu, fiber, matrix].iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("failure index overflow; check the stresses and strengths".into()));
    }
    Ok(FailureIndices {
        maximum_stress,
        tsai_hill,
        tsai_wu,
        hashin: HashinIndices {
            fiber_tension: (s1 >= 0.0).then_some(fiber),
            fiber_compression: (s1 < 0.0).then_some(fiber),
            matrix_tension: (s2 >= 0.0).then_some(matrix),
            matrix_compression: (s2 < 0.0).then_some(matrix),
            maximum: fiber.max(matrix),
        },
    })
}
