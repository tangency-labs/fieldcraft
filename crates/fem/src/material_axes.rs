//! Spatially varying material axes in linear solid elements.
//!
//! The displacement form is K = ∫ Bᵀ D(x) B dV (Hughes, *The Finite Element Method*,
//! 2000, ch. 3), with the tensor transformation of Jones, *Mechanics of Composite
//! Materials*, 2nd ed., 1999, §2.4. A cylindrical frame is evaluated at each stiffness
//! point, not at the element centre. Its components are radial, circumferential, axial;
//! the returned fields remain in global axes, with engineering shear strains.

use fieldcraft_model::{Elastic, MaterialOrientation};

use crate::dense::Matrix;
use crate::element::{Formulation, PointStress, Solid, accumulate_symmetric, b_block};
use crate::material::{elasticity_in_axes, stress};
use crate::{Error, Result};

impl Solid<'_> {
    /// Linear stiffness with material directions defined in space. Rectangular axes use
    /// the existing constant-elasticity formulation. Cylindrical axes vary at every
    /// stiffness point and currently support full and reduced displacement formulations;
    /// their incompatible-mode and hybrid formulations are refused explicitly.
    pub fn stiffness_in_material_axes(&self, elastic: &Elastic, orientation: &MaterialOrientation) -> Result<Matrix> {
        elastic.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        orientation.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        let k = match orientation {
            MaterialOrientation::Rectangular { .. } => {
                let axes = orientation.axes_at([0.0; 3]).map_err(|e| Error::Invalid(e.to_string()))?;
                self.stiffness(&elasticity_in_axes(elastic, &axes)?)?
            }
            MaterialOrientation::Cylindrical { .. } => {
                self.check_varying_material_axes()?;
                let mut k = Matrix::zeros(3 * self.ty.nodes);
                for p in self.rule().points() {
                    let g = self.gradients(p.xi)?;
                    let axes = orientation.axes_at(self.position_of(g.shape.values())).map_err(|e| Error::Invalid(e.to_string()))?;
                    let d = elasticity_in_axes(elastic, &axes)?;
                    accumulate_symmetric(&mut k.data, g.dx.get(..self.ty.nodes).unwrap_or_default(), &d, g.det * p.weight);
                }
                k
            }
        };
        if k.data.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("solid stiffness with material axes is not finite".into()));
        }
        Ok(k)
    }

    /// Global strain and stress at the stiffness points, using the material frame at each
    /// point. Strains are measured from the reference configuration, without thermal strain.
    /// `initial`, when present, is a global stress at each node, interpolated by the element
    /// shape functions: σ = D(x) ε + σ₀. No constitutive state or material damage is evolved.
    pub fn stresses_in_material_axes(
        &self,
        elastic: &Elastic,
        orientation: &MaterialOrientation,
        u: &[[f64; 3]],
        initial: Option<&[[f64; 6]]>,
    ) -> Result<Vec<PointStress>> {
        elastic.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        orientation.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        if u.len() != self.ty.nodes || u.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid(format!("material-axis solid stresses need one finite displacement per node ({} nodes)", self.ty.nodes)));
        }
        if let Some(s) = initial
            && (s.len() != self.ty.nodes || s.iter().flatten().any(|v| !v.is_finite()))
        {
            return Err(Error::Invalid("material-axis solid stresses need one finite initial stress per node".into()));
        }
        let points = match orientation {
            MaterialOrientation::Rectangular { .. } => {
                let axes = orientation.axes_at([0.0; 3]).map_err(|e| Error::Invalid(e.to_string()))?;
                self.stresses_from(&elasticity_in_axes(elastic, &axes)?, u, None, initial)?
            }
            MaterialOrientation::Cylindrical { .. } => {
                self.check_varying_material_axes()?;
                let mut points = Vec::with_capacity(self.rule().points().len());
                for p in self.rule().points() {
                    let g = self.gradients(p.xi)?;
                    let xyz = self.position_of(g.shape.values());
                    let axes = orientation.axes_at(xyz).map_err(|e| Error::Invalid(e.to_string()))?;
                    let d = elasticity_in_axes(elastic, &axes)?;
                    let mut strain = [0.0; 6];
                    for (gradient, displacement) in g.dx.iter().zip(u) {
                        for (value, row) in strain.iter_mut().zip(b_block(gradient)) {
                            *value += row.iter().zip(displacement).map(|(b, u)| b * u).sum::<f64>();
                        }
                    }
                    let mut sigma = stress(&d, &strain);
                    if let Some(s) = initial {
                        for (v, s0) in sigma.iter_mut().zip(Self::initial_stress_at(s, g.shape.values())) {
                            *v += s0;
                        }
                    }
                    points.push(PointStress { xyz, strain, stress: sigma });
                }
                points
            }
        };
        if points.iter().any(|p| p.xyz.iter().chain(&p.strain).chain(&p.stress).any(|v| !v.is_finite())) {
            return Err(Error::Invalid("solid strain or stress with material axes is not finite".into()));
        }
        Ok(points)
    }

    fn check_varying_material_axes(&self) -> Result<()> {
        match self.formulation {
            Formulation::Full | Formulation::Reduced => Ok(()),
            formulation => Err(Error::NotSupported(format!("spatially varying material axes with the {formulation:?} solid formulation"))),
        }
    }
}
