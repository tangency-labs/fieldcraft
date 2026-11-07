//! Classical laminate constitutive integration, shared by the MITC shell kernels.
//!
//! Written from A. T. Nettles, *Basic Mechanics of Laminated Composite Plates*, NASA RP-1351
//! (1994), §§II–IV, <https://ntrs.nasa.gov/citations/19950009349>. The shell uses engineering
//! shear and ε(z) = ε₀ + zκ. Plies run bottom to top along the right-hand element normal;
//! positive angles turn the shell x axis towards y. This explicitly fixes the conventions
//! that differ between texts. A, B, D are ∫Q̄ dz, ∫zQ̄ dz, ∫z²Q̄ dz about the nodal surface.
//!
//! Each homogeneous ply is integrated analytically, including its own t³/12 term. The section
//! point count controls sampling, not stiffness accuracy. Reissner–Mindlin shear uses the
//! existing shell's 5/6 correction in both directions; this equivalent-single-layer model
//! does not recover equilibrated interlaminar stresses or a traction-free shear profile.

use fieldcraft_model::{Elastic, Layup, Model, PlyStrength};

use crate::failure::{FailureIndices, plane_stress_indices};
use crate::{Error, Result};

pub type InPlane = [[f64; 3]; 3];

/// A ply resolved from the model, with the material angle expressed in the shell's axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlyProperties {
    pub elastic: Elastic,
    pub thickness: f64,
    pub angle: f64,
    pub points: u32,
    pub density: Option<f64>,
    pub strength: Option<PlyStrength>,
}

/// Resolve the materials of an ordered model layup without assigning an element frame.
/// Angles remain relative to the section reference direction; [`crate::shell::Shell::laminated`]
/// composes that direction with the element frame. No missing density or strength is invented.
pub fn resolve_plies(layup: &Layup, model: &Model) -> Result<Vec<PlyProperties>> {
    layup.validate().map_err(|e| Error::Invalid(e.to_string()))?;
    layup
        .plies
        .iter()
        .map(|ply| {
            let material = model.material(&ply.material).ok_or_else(|| Error::Invalid(format!("ply material `{}` is missing", ply.material)))?;
            let elastic = material.elastic.ok_or_else(|| Error::Invalid(format!("ply material `{}` has no elastic constants", ply.material)))?;
            Ok(PlyProperties {
                elastic,
                thickness: ply.thickness,
                angle: ply.angle,
                points: ply.points,
                density: material.density,
                strength: material.ply_strength,
            })
        })
        .collect()
}

/// One ply's in-plane constitutive relation in its own material axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lamina {
    q: InPlane,
    shear: [f64; 2],
    /// ε33 = coefficients · [σ11, σ22] under plane stress.
    normal_compliance: [f64; 2],
}

impl Lamina {
    pub fn new(elastic: &Elastic) -> Result<Self> {
        elastic.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        let (e1, e2, nu12, g12, g13, g23, normal_compliance) = match *elastic {
            Elastic::Isotropic { young, poisson } => {
                let g = young / (2.0 * (1.0 + poisson));
                (young, young, poisson, g, g, g, [-poisson / young; 2])
            }
            Elastic::EngineeringConstants { young, poisson, shear } => {
                (young[0], young[1], poisson[0], shear[0], shear[1], shear[2], [-poisson[1] / young[0], -poisson[2] / young[1]])
            }
        };
        let denominator = 1.0 - nu12 * (nu12 * e2 / e1);
        let q11 = e1 / denominator;
        let q22 = e2 / denominator;
        let q12 = nu12 * q22;
        let q = [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, g12]];
        if !(denominator > 0.0)
            || q.iter().flatten().chain(&normal_compliance).any(|v| !v.is_finite())
            || ![g12, g13, g23].iter().all(|v| *v > 0.0 && v.is_finite())
        {
            return Err(Error::Invalid("lamina stiffness is not finite and positive definite".into()));
        }
        Ok(Self { q, shear: [g13, g23], normal_compliance })
    }

    pub fn plane_stress(&self) -> InPlane {
        self.q
    }

    /// Q̄ = Tεᵀ Q Tε: energy conjugacy supplies the stress transformation, avoiding the
    /// common factor-of-two error from transforming engineering shear like stress.
    pub fn rotated(&self, angle: f64) -> Result<InPlane> {
        let t = strain_transform(angle)?;
        let mut q = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let value = (0..3).map(|a| (0..3).map(|b| t[a][i] * self.q[a][b] * t[b][j]).sum::<f64>()).sum();
                q[i][j] = value;
                q[j][i] = value;
            }
        }
        if q.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("rotated lamina stiffness is not finite".into()));
        }
        Ok(q)
    }
}

fn sine_cosine(angle: f64) -> Result<(f64, f64)> {
    if !angle.is_finite() {
        return Err(Error::Invalid("a ply angle must be finite".into()));
    }
    Ok(libm::sincos(angle))
}

/// [ε11, ε22, γ12] in material axes from [εxx, εyy, γxy] in shell axes.
pub fn strain_transform(angle: f64) -> Result<InPlane> {
    let (s, c) = sine_cosine(angle)?;
    Ok([[c * c, s * s, c * s], [s * s, c * c, -c * s], [-2.0 * c * s, 2.0 * c * s, c * c - s * s]])
}

fn product(a: &InPlane, v: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i].iter().zip(v).map(|(x, y)| x * y).sum())
}

#[derive(Debug, Clone, PartialEq)]
struct IntegratedPly {
    input: PlyProperties,
    lamina: Lamina,
    bottom: f64,
    middle: f64,
    top: f64,
}

/// Integrated section about the reference surface. All matrices use shell-local axes.
#[derive(Debug, Clone, PartialEq)]
pub struct Laminate {
    plies: Vec<IntegratedPly>,
    abd: [[f64; 6]; 6],
    shear: [[f64; 2]; 2],
    mass: Option<[f64; 3]>,
    geometry: [f64; 3],
    drilling: f64,
    thickness: f64,
    eccentricity: f64,
}

impl Laminate {
    pub fn new(plies: &[PlyProperties], offset: f64) -> Result<Self> {
        if plies.is_empty() || plies.len() > fieldcraft_model::composite::MAX_PLIES {
            return Err(Error::Invalid(format!("a laminate needs 1 to {} plies", fieldcraft_model::composite::MAX_PLIES)));
        }
        if !(offset.abs() <= 0.5) {
            return Err(Error::Invalid("a laminate offset must be in -0.5..=0.5".into()));
        }
        if !plies.iter().all(|p| p.thickness > 0.0 && p.thickness.is_finite()) {
            return Err(Error::Invalid("ply thicknesses must be positive and finite".into()));
        }
        let thickness: f64 = plies.iter().map(|p| p.thickness).sum();
        if !thickness.is_finite() {
            return Err(Error::Invalid("total laminate thickness is not finite".into()));
        }
        let eccentricity = -offset * thickness;
        let mut out = Self {
            plies: Vec::with_capacity(plies.len()),
            abd: [[0.0; 6]; 6],
            shear: [[0.0; 2]; 2],
            mass: Some([0.0; 3]),
            geometry: [0.0; 3],
            drilling: 0.0,
            thickness,
            eccentricity,
        };
        let mut bottom = eccentricity - thickness / 2.0;
        for p in plies {
            fieldcraft_model::composite::check_ply_points(p.points).map_err(|e| Error::Invalid(e.to_string()))?;
            if let Some(s) = p.strength {
                s.validate().map_err(|e| Error::Invalid(e.to_string()))?;
            }
            if p.density.is_some_and(|d| !(d > 0.0 && d.is_finite())) {
                return Err(Error::Invalid("ply density must be positive and finite".into()));
            }
            let lamina = Lamina::new(&p.elastic)?;
            let q = lamina.rotated(p.angle)?;
            let top = bottom + p.thickness;
            let middle = bottom + p.thickness / 2.0;
            if !(top > bottom && top.is_finite()) {
                return Err(Error::Invalid("ply thickness is not representable at its section position".into()));
            }
            // Centred moments avoid subtracting nearly equal powers at thin ply interfaces.
            let moments = [p.thickness, p.thickness * middle, p.thickness * (middle * middle + p.thickness * p.thickness / 12.0)];
            for i in 0..3 {
                out.geometry[i] += moments[i];
                for j in 0..3 {
                    out.abd[i][j] += q[i][j] * moments[0];
                    out.abd[i][j + 3] += q[i][j] * moments[1];
                    out.abd[i + 3][j] += q[i][j] * moments[1];
                    out.abd[i + 3][j + 3] += q[i][j] * moments[2];
                }
            }
            out.mass = match (out.mass, p.density) {
                (Some(mass), Some(density)) => Some(std::array::from_fn(|i| mass[i] + density * moments[i])),
                _ => None,
            };
            let (s, c) = sine_cosine(p.angle)?;
            let [g13, g23] = lamina.shear;
            let shear = [[g13 * c * c + g23 * s * s, (g13 - g23) * c * s], [(g13 - g23) * c * s, g13 * s * s + g23 * c * c]];
            for (row, values) in out.shear.iter_mut().zip(shear) {
                for (v, value) in row.iter_mut().zip(values) {
                    *v += crate::shell::SHEAR_CORRECTION * p.thickness * value;
                }
            }
            out.drilling += crate::shell::DRILLING * lamina.q[2][2] * p.thickness;
            out.plies.push(IntegratedPly { input: *p, lamina, bottom, middle, top });
            bottom = top;
        }
        if out.abd.iter().flatten().chain(out.shear.iter().flatten()).chain(&out.geometry).chain(out.mass.iter().flatten()).any(|v| !v.is_finite())
            || !(out.drilling > 0.0 && out.drilling.is_finite())
            || (0..6).any(|i| out.abd[i][i] <= 0.0)
            || !positive_shear(out.shear)
            || out.mass.is_some_and(|m| !(m[0] > 0.0 && m[2] > 0.0))
        {
            return Err(Error::Invalid("laminate stiffness or mass moments overflow or underflow".into()));
        }
        Ok(out)
    }

    pub fn abd(&self) -> [[f64; 6]; 6] {
        self.abd
    }
    pub fn shear(&self) -> [[f64; 2]; 2] {
        self.shear
    }
    pub fn thickness(&self) -> f64 {
        self.thickness
    }
    pub fn eccentricity(&self) -> f64 {
        self.eccentricity
    }
    pub fn drilling(&self) -> f64 {
        self.drilling
    }
    /// ∫ρ dz, ∫ρz dz, ∫ρz² dz. Missing density in any ply is an error, not zero mass.
    pub fn mass_moments(&self) -> Result<[f64; 3]> {
        self.mass.ok_or_else(|| Error::Invalid("every laminate ply needs a density for mass or gravity".into()))
    }
    /// ∫1 dz, ∫z dz, ∫z² dz, for uniform force per volume independent of mass density.
    pub fn geometric_moments(&self) -> [f64; 3] {
        self.geometry
    }

    /// Forces per unit length [Nx, Ny, Nxy, Mx, My, Mxy, Qx, Qy] for generalized shell strains.
    pub fn resultants(&self, strain: [f64; 8]) -> Result<[f64; 8]> {
        check_strain(strain)?;
        let mut out = [0.0; 8];
        for (value, row) in out.iter_mut().zip(self.abd) {
            *value = row.iter().zip(strain).map(|(c, e)| c * e).sum();
        }
        for i in 0..2 {
            out[6 + i] = self.shear[i][0] * strain[6] + self.shear[i][1] * strain[7];
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("laminate resultants are not finite".into()));
        }
        Ok(out)
    }

    /// Each interface has two distinct one-sided stress values; never average across plies.
    pub fn recover(&self, strain: [f64; 8]) -> Result<Vec<PlyFields>> {
        check_strain(strain)?;
        self.plies
            .iter()
            .map(|p| {
                Ok(PlyFields {
                    angle: p.input.angle,
                    bottom: p.at(strain, p.bottom)?,
                    middle: p.at(strain, p.middle)?,
                    top: p.at(strain, p.top)?,
                    points: section_abscissae(p.input.points)?
                        .iter()
                        .map(|r| p.at(strain, p.middle + r * p.input.thickness / 2.0))
                        .collect::<Result<_>>()?,
                })
            })
            .collect()
    }

    /// One position in the laminate. At a shared interface this selects the lower ply;
    /// [`Self::recover`] retains both one-sided values instead.
    pub fn at(&self, strain: [f64; 8], z: f64) -> Result<(f64, PlyPoint)> {
        check_strain(strain)?;
        let ply = self
            .plies
            .iter()
            .find(|p| z >= p.bottom && z <= p.top)
            .ok_or_else(|| Error::Invalid("the requested position is outside the laminate".into()))?;
        Ok((ply.input.angle, ply.at(strain, z)?))
    }
}

fn positive_shear(shear: [[f64; 2]; 2]) -> bool {
    shear[0][0] > 0.0 && shear[1][1] > 0.0 && (shear[0][1] / shear[0][0].sqrt() / shear[1][1].sqrt()).abs() < 1.0
}

/// Gauss–Legendre section-point positions, 1..=15. Newton's method on the Legendre
/// recurrence (Abramowitz & Stegun, 1964, §22.7), bounded to 32 updates per root.
fn section_abscissae(points: u32) -> Result<Vec<f64>> {
    fieldcraft_model::composite::check_ply_points(points).map_err(|e| Error::Invalid(e.to_string()))?;
    let n = points as usize;
    let mut out = vec![0.0; n];
    for i in 0..n / 2 {
        let mut x = libm::cos(std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5));
        let mut converged = false;
        for _ in 0..32 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p;
            }
            let derivative = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / derivative;
            x -= dx;
            if dx.abs() <= 4.0 * f64::EPSILON {
                converged = true;
                break;
            }
        }
        if !converged || !x.is_finite() {
            return Err(Error::Invalid("laminate section-point quadrature did not converge".into()));
        }
        out[i] = -x;
        out[n - 1 - i] = x;
    }
    Ok(out)
}

fn check_strain(strain: [f64; 8]) -> Result<()> {
    if strain.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("generalized laminate strains must be finite".into()));
    }
    Ok(())
}

/// Material-axis fields at one physical thickness position. Stress/strain order is
/// 11,22,33,12,13,23, with engineering shear strains. Transverse shear is the FSDT effective
/// stress (5/6 Gγ), not an equilibrated interlaminar stress. Failure is evaluated on 11,22,12.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PlyPoint {
    pub z: f64,
    pub strain: [f64; 6],
    pub stress: [f64; 6],
    pub failure: Option<FailureIndices>,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PlyFields {
    pub angle: f64,
    pub bottom: PlyPoint,
    pub middle: PlyPoint,
    pub top: PlyPoint,
    /// Additional Gauss–Legendre section samples, bottom to top, one to fifteen per ply.
    pub points: Vec<PlyPoint>,
}

impl IntegratedPly {
    fn at(&self, g: [f64; 8], z: f64) -> Result<PlyPoint> {
        let e = product(&strain_transform(self.input.angle)?, [g[0] + z * g[3], g[1] + z * g[4], g[2] + z * g[5]]);
        let s = product(&self.lamina.q, e);
        let (sn, cs) = sine_cosine(self.input.angle)?;
        let transverse = [cs * g[6] + sn * g[7], -sn * g[6] + cs * g[7]];
        let strain =
            [e[0], e[1], self.lamina.normal_compliance[0] * s[0] + self.lamina.normal_compliance[1] * s[1], e[2], transverse[0], transverse[1]];
        let stress = [
            s[0],
            s[1],
            0.0,
            s[2],
            crate::shell::SHEAR_CORRECTION * self.lamina.shear[0] * transverse[0],
            crate::shell::SHEAR_CORRECTION * self.lamina.shear[1] * transverse[1],
        ];
        if strain.iter().chain(&stress).any(|v| !v.is_finite()) {
            return Err(Error::Invalid("recovered ply fields are not finite".into()));
        }
        let failure = self.input.strength.map(|strength| plane_stress_indices(s, &strength)).transpose()?;
        Ok(PlyPoint { z, strain, stress, failure })
    }
}
