//! Solver-neutral laminate data and material directions.
//!
//! Plies are ordered from the negative to the positive shell normal. Angles are radians,
//! right handed about that normal. The reference surface follows [`crate::SectionKind::Shell`]:
//! the laminate centre is at `-offset * thickness`. See A. T. Nettles, NASA RP-1351 (1994),
//! §§II–IV, for material axes and laminate constitutive relations; our stack order is explicit
//! and does not depend on a reference publication's viewing direction.

use serde::{Deserialize, Serialize};

use crate::{Error, Result, check_name, finite, positive};

/// Bound storage and all through-thickness loops before allocating resolved plies.
pub const MAX_PLIES: usize = 256;

/// A section's material frame, including the axis projected into a shell's plane.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SectionOrientation {
    pub frame: MaterialOrientation,
    #[serde(default = "default_reference_axis")]
    pub shell_reference_axis: u8,
}

fn default_reference_axis() -> u8 {
    1
}

impl SectionOrientation {
    pub fn validate(&self) -> Result<()> {
        self.frame.validate()?;
        if !(1..=3).contains(&self.shell_reference_axis) {
            return Err(Error::Invalid("a shell reference axis must be 1, 2 or 3".into()));
        }
        Ok(())
    }
}

/// Material axes in global coordinates. The section's selected reference axis is projected
/// onto a shell's plane before applying a ply angle; solid elements use all three axes at the evaluation position.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum MaterialOrientation {
    /// The first direction and a reference direction in the 1–2 plane (not necessarily unit
    /// or perpendicular). Axis 3 is their right-hand normal.
    Rectangular { axis1: [f64; 3], axis2: [f64; 3] },
    /// Radial, circumferential, axial directions about the line `origin + t * axis`.
    /// `origin` is in metres, `axis` is a dimensionless direction. Points on the axis have
    /// no radial direction and are refused.
    Cylindrical { origin: [f64; 3], axis: [f64; 3] },
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn direction(v: [f64; 3]) -> Result<[f64; 3]> {
    let scale = v.iter().map(|c| c.abs()).fold(0.0, f64::max);
    if v.iter().any(|c| !c.is_finite()) || scale == 0.0 {
        return Err(Error::Invalid("a material direction must be finite and nonzero".into()));
    }
    let v = v.map(|c| c / scale);
    let length = dot(v, v).sqrt();
    Ok(v.map(|c| c / length))
}

fn transverse(v: [f64; 3], normal: [f64; 3]) -> Result<[f64; 3]> {
    let v = direction(v)?;
    let perpendicular = cross(normal, v);
    if dot(perpendicular, perpendicular) <= 1e-24 {
        return Err(Error::Invalid("material reference directions are parallel or numerically indistinguishable".into()));
    }
    // Reconstruct by a second cross product: subtracting nearly parallel vectors once
    // leaves an amplified component along the normal after normalization.
    direction(cross(direction(perpendicular)?, normal))
}

impl MaterialOrientation {
    pub fn validate(&self) -> Result<()> {
        match *self {
            Self::Rectangular { axis1, axis2 } => {
                transverse(axis2, direction(axis1)?)?;
            }
            Self::Cylindrical { origin, axis } => {
                crate::check_coordinates(&origin)?;
                direction(axis)?;
            }
        }
        Ok(())
    }

    /// Right-handed orthonormal axes, as rows in global coordinates. No arbitrary fallback
    /// is chosen for a degenerate reference or a point on a cylindrical axis.
    pub fn axes_at(&self, position: [f64; 3]) -> Result<[[f64; 3]; 3]> {
        self.validate()?;
        crate::check_coordinates(&position)?;
        match *self {
            Self::Rectangular { axis1, axis2 } => {
                let e1 = direction(axis1)?;
                let e2 = transverse(axis2, e1)?;
                Ok([e1, e2, direction(cross(e1, e2))?])
            }
            Self::Cylindrical { origin, axis } => {
                let e3 = direction(axis)?;
                let displacement = direction(std::array::from_fn(|i| position[i] - origin[i]))?;
                // Unlike two reference directions, a point far along the cylinder axis
                // can have a perfectly valid small radius. Do not use an angular cutoff.
                let hoop = direction(cross(e3, displacement))?;
                let radial = direction(cross(hoop, e3))?;
                Ok([radial, hoop, e3])
            }
        }
    }

    /// Shell material axes before the ply rotation. Positive axis 3 is the element normal.
    /// `reference_axis` is 1, 2 or 3; choose 2 (hoop) or 3 (axial) on a cylindrical shell,
    /// where the radial axis cannot define an in-plane fibre direction.
    pub fn shell_axes_at(&self, position: [f64; 3], normal: [f64; 3], reference_axis: u8) -> Result<[[f64; 3]; 3]> {
        let normal = direction(normal)?;
        let axes = self.axes_at(position)?;
        let reference = reference_axis
            .checked_sub(1)
            .and_then(|i| axes.get(usize::from(i)))
            .ok_or_else(|| Error::Invalid("a shell reference axis must be 1, 2 or 3".into()))?;
        let e1 = transverse(*reference, normal)?;
        Ok([e1, direction(cross(normal, e1))?, normal])
    }
}

/// One homogeneous ply. `points` controls sampling through its thickness; the elastic
/// stiffness and mass moments can be integrated exactly even with one sampling point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Ply {
    pub material: String,
    pub thickness: f64,
    #[serde(default)]
    pub angle: f64,
    #[serde(default = "default_ply_points")]
    pub points: u32,
}

fn default_ply_points() -> u32 {
    3
}

impl Ply {
    pub fn validate(&self) -> Result<()> {
        check_name("ply material", &self.material)?;
        positive("ply thickness", self.thickness)?;
        finite("ply angle", self.angle)?;
        check_ply_points(self.points)
    }
}

pub fn check_ply_points(points: u32) -> Result<()> {
    if !(1..=crate::sections::MAX_SHELL_POINTS).contains(&points) {
        return Err(Error::Invalid(format!("a ply needs 1 to {} section points, got {points}", crate::sections::MAX_SHELL_POINTS)));
    }
    Ok(())
}

/// A laminate's ordered plies and offset. Total thickness is derived, never a second input.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Layup {
    pub plies: Vec<Ply>,
    #[serde(default)]
    pub offset: f64,
}

impl Layup {
    pub fn thickness(&self) -> f64 {
        self.plies.iter().map(|p| p.thickness).sum()
    }

    pub fn validate(&self) -> Result<()> {
        if self.plies.is_empty() || self.plies.len() > MAX_PLIES {
            return Err(Error::Invalid(format!("a laminate needs 1 to {MAX_PLIES} plies")));
        }
        if !(self.offset.abs() <= 0.5) {
            return Err(Error::Invalid("a laminate offset must be in -0.5..=0.5".into()));
        }
        for (i, ply) in self.plies.iter().enumerate() {
            ply.validate().map_err(|e| Error::Invalid(format!("ply {}: {e}", i + 1)))?;
        }
        positive("total laminate thickness", self.thickness())?;
        Ok(())
    }
}

/// Plane-stress strengths, all positive magnitudes in Pa. S23 is independent of S12: the
/// 1980 Hashin matrix-compression branch needs it even when the applied stress is plane stress.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlyStrength {
    pub tension: [f64; 2],
    pub compression: [f64; 2],
    pub shear12: f64,
    pub shear23: f64,
    /// F12 / sqrt(F11 F22), dimensionless, strictly between -1 and 1 for a closed Tsai–Wu
    /// ellipse. The default -0.5 is an assumption, not a measured biaxial strength.
    #[serde(default = "default_interaction")]
    pub tsai_wu_interaction: f64,
}

fn default_interaction() -> f64 {
    -0.5
}

impl PlyStrength {
    pub fn validate(&self) -> Result<()> {
        for (name, v) in [
            ("Xt", self.tension[0]),
            ("Yt", self.tension[1]),
            ("Xc", self.compression[0]),
            ("Yc", self.compression[1]),
            ("S12", self.shear12),
            ("S23", self.shear23),
        ] {
            positive(name, v)?;
        }
        if !(self.tsai_wu_interaction.abs() < 1.0) {
            return Err(Error::Invalid("Tsai–Wu normalized interaction must be strictly between -1 and 1".into()));
        }
        Ok(())
    }
}
