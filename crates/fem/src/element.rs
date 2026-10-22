//! Solid continuum elements: the element table, stiffness, mass, body forces, stresses at the
//! integration points and their extrapolation to the nodes.
//!
//! The isoparametric displacement formulation (Zienkiewicz, Taylor & Zhu 2013, ch. 6; Hughes
//! 2000, ch. 3): K = ∫ Bᵀ D B dV and M = ∫ ρ Nᵀ N dV, integrated by Gauss quadrature on the
//! reference element through the Jacobian of the map from natural to physical coordinates.
//!
//! Formulations:
//! - **full integration**: the rule that integrates the stiffness of an undistorted element
//!   exactly (1 point for tet4, 4 for tet10, 3 triangle points × 2 for wedge6, 2 × 2 × 2 for
//!   hex8, 3 × 3 × 3 for hex20);
//! - **reduced integration** for hex20 (2 × 2 × 2), which avoids locking at the cost of one
//!   spurious mode that a mesh of more than one element restrains. Reduced hex8 needs hourglass
//!   control and is not offered here yet (it arrives with the explicit solver, Flanagan &
//!   Belytschko 1981);
//! - **incompatible modes** for hex8 (E. L. Wilson, R. L. Taylor, W. P. Doherty, J. Ghaboussi,
//!   *Incompatible displacement models*, in Numerical and Computer Methods in Structural
//!   Mechanics, Academic Press 1973, 43–57): nine internal displacement modes (1 − ξ²),
//!   (1 − η²), (1 − ζ²) per direction, condensed out of the element. Their derivatives use the
//!   Jacobian at the element centre scaled by det J₀ / det J, the correction of R. L. Taylor,
//!   P. J. Beresford, E. L. Wilson, *A non-conforming element for stress analysis*, IJNME
//!   10(6):1211–1219, 1976, so the element passes the patch test when distorted. It bends
//!   correctly with one element through the thickness.
//!
//! Lumped mass is the diagonal scaling of Hinton, Rock and Zienkiewicz (E. Hinton, T. Rock,
//! O. C. Zienkiewicz, *A note on mass lumping and related processes in the finite element
//! method*, Earthquake Eng. Struct. Dyn. 4(3):245–249, 1976): the consistent diagonal scaled to
//! the element's mass, positive for quadratic elements too.
//!
//! Nodal stresses extrapolate the integration-point values with the element's linear corner
//! functions fitted by least squares (E. Hinton, J. S. Campbell, *Local and global smoothing of
//! discontinuous finite element functions using a least squares method*, IJNME 8(3):461–480,
//! 1974); mid-side nodes take the mean of their edge's corners.

use fieldcraft_model::ElementKind;

use crate::basis::Basis;
use crate::dense::{Cholesky, Matrix, positive_inverse};
use crate::dof::{DofSet, ElementDofs, Internal};
use crate::material::{Elasticity, stress};
use crate::quadrature::{Point, Rule};
use crate::shape::{self, MAX_NODES, Shape};
use crate::{Error, Result};

/// The reference shape of an element type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    Tetrahedron,
    Wedge,
    Pyramid,
    Hexahedron,
    /// The 2D elements of plane and axisymmetric models ([`crate::plane`]).
    Triangle,
    Quadrilateral,
}

/// One row of the element table: everything the kernels need to know about an element kind.
/// Adding an element (wedge6, shells, beams) is adding a row and its shape functions.
#[derive(Debug)]
pub struct ElementType {
    pub kind: ElementKind,
    pub family: Family,
    pub nodes: usize,
    /// The corner nodes come first in the node list.
    pub corners: usize,
    pub shape: fn([f64; 3]) -> Shape,
    /// Full integration of the stiffness.
    pub full: Rule,
    /// Reduced integration, where the element is stable with it.
    pub reduced: Option<Rule>,
    /// Exact for the consistent mass of an undistorted element.
    pub mass: Rule,
    /// Whether the element has an incompatible-modes formulation.
    pub incompatible_modes: bool,
}

/// The element table.
pub static ELEMENT_TYPES: [ElementType; 6] = [
    ElementType {
        kind: ElementKind::Tet4,
        family: Family::Tetrahedron,
        nodes: 4,
        corners: 4,
        shape: shape::tet4,
        full: Rule::TetCentroid,
        reduced: None,
        mass: Rule::Tet4,
        incompatible_modes: false,
    },
    ElementType {
        kind: ElementKind::Tet10,
        family: Family::Tetrahedron,
        nodes: 10,
        corners: 4,
        shape: shape::tet10,
        full: Rule::Tet4,
        reduced: None,
        mass: Rule::TetConical(4),
        incompatible_modes: false,
    },
    ElementType {
        kind: ElementKind::Wedge6,
        family: Family::Wedge,
        nodes: 6,
        corners: 6,
        shape: shape::wedge6,
        full: Rule::Wedge,
        reduced: None,
        mass: Rule::Wedge,
        incompatible_modes: false,
    },
    ElementType {
        kind: ElementKind::Pyramid5,
        family: Family::Pyramid,
        nodes: 5,
        corners: 5,
        shape: shape::pyramid5,
        full: Rule::Pyramid(2),
        reduced: None,
        mass: Rule::Pyramid(3),
        incompatible_modes: false,
    },
    ElementType {
        kind: ElementKind::Hex8,
        family: Family::Hexahedron,
        nodes: 8,
        corners: 8,
        shape: shape::hex8,
        full: Rule::Hex(2),
        reduced: None,
        mass: Rule::Hex(2),
        incompatible_modes: true,
    },
    ElementType {
        kind: ElementKind::Hex20,
        family: Family::Hexahedron,
        nodes: 20,
        corners: 8,
        shape: shape::hex20,
        full: Rule::Hex(3),
        reduced: Some(Rule::Hex(2)),
        mass: Rule::Hex(3),
        incompatible_modes: false,
    },
];

/// The table row of an element kind.
pub fn element_type(kind: ElementKind) -> Result<&'static ElementType> {
    ELEMENT_TYPES
        .iter()
        .find(|t| t.kind == kind)
        .ok_or_else(|| Error::NotSupported(format!("{} elements (the solver has tet4, tet10, wedge6, pyramid5, hex8 and hex20 so far)", kind.name())))
}

/// How a solid element is formulated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Formulation {
    Full,
    Reduced,
    IncompatibleModes,
    /// Mixed displacement–pressure with `pressures` discontinuous pressure unknowns condensed
    /// out of the element ([`crate::hybrid`]).
    Hybrid {
        pressures: usize,
    },
}

impl Formulation {
    /// The unknowns of a solid element so formulated: three translations per node, and the
    /// incompatible modes it condenses out.
    pub fn dofs(self) -> ElementDofs {
        let incompatible_modes = if self == Formulation::IncompatibleModes { INCOMPATIBLE_DOFS } else { 0 };
        let pressure = if let Formulation::Hybrid { pressures } = self { pressures } else { 0 };
        ElementDofs { node: DofSet::TRANSLATIONS, internal: Internal { incompatible_modes, pressure } }
    }

    /// The formulation a model section asks for, if this element type has it.
    pub fn from_section(kind: ElementKind, f: &fieldcraft_model::Formulation) -> Result<Formulation> {
        let ty = element_type(kind)?;
        if f.hybrid {
            if f.incompatible_modes || f.integration == fieldcraft_model::Integration::Reduced {
                return Err(Error::NotSupported("hybrid elements with incompatible modes or reduced integration".into()));
            }
            return crate::hybrid::pressures(kind).map(|pressures| Formulation::Hybrid { pressures }).ok_or_else(|| {
                Error::NotSupported(format!(
                    "hybrid {} elements (with one pressure per element the linear tetrahedron and wedge still lock; use hex8, hex20 or tet10)",
                    kind.name()
                ))
            });
        }
        if f.incompatible_modes {
            return if ty.incompatible_modes {
                Ok(Formulation::IncompatibleModes)
            } else {
                Err(Error::Invalid(format!("incompatible modes apply to hex8 elements, not {}", kind.name())))
            };
        }
        match f.integration {
            fieldcraft_model::Integration::Full => Ok(Formulation::Full),
            fieldcraft_model::Integration::Reduced if ty.reduced.is_some() => Ok(Formulation::Reduced),
            fieldcraft_model::Integration::Reduced => Err(Error::NotSupported(format!(
                "reduced integration of {} elements (hex8 needs hourglass control; use full integration or incompatible modes)",
                kind.name()
            ))),
        }
    }
}

/// Strain, stress and position at one integration point.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PointStress {
    pub xyz: [f64; 3],
    pub strain: [f64; 6],
    pub stress: [f64; 6],
}

/// Shape function gradients in physical coordinates at one point.
pub(crate) struct Gradients {
    pub(crate) shape: Shape,
    pub(crate) dx: [[f64; 3]; MAX_NODES],
    pub(crate) det: f64,
}

/// A thermal strain field over an element: the expansion coefficients along X, Y, Z and the
/// nodal temperature changes T − T_ref.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalStrain {
    pub alpha: [f64; 3],
    pub delta: Vec<f64>,
}

impl ThermalStrain {
    pub(crate) fn check(&self, nodes: usize) -> Result<()> {
        if self.delta.len() != nodes {
            return Err(Error::Invalid(format!("{} temperatures for a {nodes}-node element", self.delta.len())));
        }
        if self.alpha.iter().chain(&self.delta).any(|v| !v.is_finite()) {
            return Err(Error::Invalid("temperatures and expansion coefficients must be finite".into()));
        }
        Ok(())
    }
}

/// A solid element: its type, its node coordinates and its formulation.
#[derive(Debug, Clone, Copy)]
pub struct Solid<'a> {
    pub ty: &'static ElementType,
    pub(crate) xyz: &'a [[f64; 3]],
    pub(crate) formulation: Formulation,
}

/// The stiffness blocks of the incompatible-modes element.
struct IncompatibleParts {
    kuu: Matrix,
    kua: Vec<f64>,
    kaa: Matrix,
    j0: ([[f64; 3]; 3], f64),
}

/// The incompatible modes' count of internal degrees of freedom.
const INCOMPATIBLE_DOFS: usize = 9;

impl<'a> Solid<'a> {
    pub fn new(kind: ElementKind, xyz: &'a [[f64; 3]], formulation: Formulation) -> Result<Solid<'a>> {
        let ty = element_type(kind)?;
        if xyz.len() != ty.nodes {
            return Err(Error::Invalid(format!("a {} element has {} nodes, got {}", kind.name(), ty.nodes, xyz.len())));
        }
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        let ok = match formulation {
            Formulation::Full => true,
            Formulation::Reduced => ty.reduced.is_some(),
            Formulation::IncompatibleModes => ty.incompatible_modes,
            Formulation::Hybrid { pressures } => crate::hybrid::pressures(kind) == Some(pressures),
        };
        if !ok {
            return Err(Error::NotSupported(format!("{formulation:?} formulation of {} elements", kind.name())));
        }
        Ok(Solid { ty, xyz, formulation })
    }

    /// The element's unknowns, in the order of its matrices.
    pub fn dofs(&self) -> ElementDofs {
        self.formulation.dofs()
    }

    /// The stiffness integration rule.
    pub fn rule(&self) -> Rule {
        match (self.formulation, self.ty.reduced) {
            (Formulation::Reduced, Some(r)) => r,
            _ => self.ty.full,
        }
    }

    /// The stiffness integration points (where stresses are reported).
    pub fn integration_points(&self) -> Vec<Point> {
        self.rule().points()
    }

    /// The Jacobian matrix J[a][b] = ∂x_b/∂ξ_a.
    fn jacobian(&self, s: &Shape) -> [[f64; 3]; 3] {
        let mut j = [[0.0; 3]; 3];
        for (d, x) in s.derivatives().iter().zip(self.xyz) {
            for a in 0..3 {
                for b in 0..3 {
                    j[a][b] += d[a] * x[b];
                }
            }
        }
        j
    }

    pub(crate) fn gradients(&self, xi: [f64; 3]) -> Result<Gradients> {
        let shape = (self.ty.shape)(xi);
        let j = self.jacobian(&shape);
        let (inv, det) = positive_inverse(&j).ok_or_else(|| self.degenerate())?;
        let mut dx = [[0.0; 3]; MAX_NODES];
        for (g, d) in dx.iter_mut().zip(shape.derivatives()) {
            for b in 0..3 {
                g[b] = inv[b][0] * d[0] + inv[b][1] * d[1] + inv[b][2] * d[2];
            }
        }
        Ok(Gradients { shape, dx, det })
    }

    fn degenerate(&self) -> Error {
        Error::Degenerate(format!("a {} element is inverted or degenerate (its Jacobian is not positive)", self.ty.kind.name()))
    }

    /// Position of a natural point.
    pub fn position(&self, xi: [f64; 3]) -> [f64; 3] {
        self.position_of((self.ty.shape)(xi).values())
    }

    /// The position Σ Nᵢ xᵢ for shape values `n`.
    pub(crate) fn position_of(&self, n: &[f64]) -> [f64; 3] {
        let mut p = [0.0; 3];
        for (ni, x) in n.iter().zip(self.xyz) {
            for k in 0..3 {
                p[k] += ni * x[k];
            }
        }
        p
    }

    /// The element's volume.
    pub fn volume(&self) -> Result<f64> {
        let mut v = 0.0;
        for p in self.ty.mass.points() {
            v += self.gradients(p.xi)?.det * p.weight;
        }
        Ok(v)
    }

    /// The incompatible modes' gradients at a point: det J₀/det J · J₀⁻¹ ∂P/∂ξ.
    fn incompatible_gradients(&self, xi: [f64; 3], det: f64, j0: &([[f64; 3]; 3], f64)) -> [[f64; 3]; 3] {
        let (inv0, det0) = j0;
        let mut g = [[0.0; 3]; 3];
        for (k, gk) in g.iter_mut().enumerate() {
            // ∂(1 − ξ_k²)/∂ξ_k = −2 ξ_k, other derivatives zero.
            let dp = -2.0 * xi[k];
            for b in 0..3 {
                gk[b] = det0 / det * inv0[b][k] * dp;
            }
        }
        g
    }

    fn centre_jacobian(&self) -> Result<([[f64; 3]; 3], f64)> {
        let s = (self.ty.shape)([0.0; 3]);
        positive_inverse(&self.jacobian(&s)).ok_or_else(|| self.degenerate())
    }

    /// The parts of the stiffness with the incompatible modes: K_uu (when `with_uu`), K_ua
    /// (3n × 9), K_aa, and the centre Jacobian they were built with.
    fn incompatible_parts(&self, d: &Elasticity, with_uu: bool) -> Result<IncompatibleParts> {
        let n = 3 * self.ty.nodes;
        let j0 = self.centre_jacobian()?;
        let mut kuu = Matrix::zeros(if with_uu { n } else { 0 });
        let mut kua = vec![0.0; n * INCOMPATIBLE_DOFS];
        let mut kaa = Matrix::zeros(INCOMPATIBLE_DOFS);
        for p in self.rule().points() {
            let g = self.gradients(p.xi)?;
            let w = g.det * p.weight;
            let ga = self.incompatible_gradients(p.xi, g.det, &j0);
            let gu = g.dx.get(..self.ty.nodes).unwrap_or_default();
            if with_uu {
                accumulate_symmetric(&mut kuu.data, gu, d, w);
            }
            accumulate(&mut kua, INCOMPATIBLE_DOFS, gu, &ga, d, w);
            accumulate_symmetric(&mut kaa.data, &ga, d, w);
        }
        Ok(IncompatibleParts { kuu, kua, kaa, j0 })
    }

    /// The stiffness matrix, 3n × 3n, degrees of freedom ordered node by node (x, y, z).
    pub fn stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        if let Formulation::Hybrid { .. } = self.formulation {
            return self.hybrid_stiffness(d);
        }
        let n = 3 * self.ty.nodes;
        if self.formulation == Formulation::IncompatibleModes {
            let IncompatibleParts { kuu: mut k, kua, kaa, .. } = self.incompatible_parts(d, true)?;
            let chol = Cholesky::new(&kaa)?;
            // K = K_uu − K_ua K_aa⁻¹ K_au, column by column of K_au.
            for i in 0..n {
                let mut col: Vec<f64> = kua.get(i * INCOMPATIBLE_DOFS..(i + 1) * INCOMPATIBLE_DOFS).unwrap_or_default().to_vec();
                chol.solve(&mut col);
                for j in 0..n {
                    let row = kua.get(j * INCOMPATIBLE_DOFS..(j + 1) * INCOMPATIBLE_DOFS).unwrap_or_default();
                    let v: f64 = row.iter().zip(&col).map(|(a, b)| a * b).sum();
                    k.add(j, i, -v);
                }
            }
            symmetrize(&mut k);
            return Ok(k);
        }
        let mut k = Matrix::zeros(n);
        for p in self.rule().points() {
            let g = self.gradients(p.xi)?;
            let gu = g.dx.get(..self.ty.nodes).unwrap_or_default();
            accumulate_symmetric(&mut k.data, gu, d, g.det * p.weight);
        }
        Ok(k)
    }

    /// The consistent mass matrix per displacement component, n × n: m_ij = ∫ ρ Nᵢ Nⱼ dV (the
    /// 3n × 3n matrix has it on each component's block).
    pub fn mass(&self, density: f64) -> Result<Matrix> {
        let n = self.ty.nodes;
        let mut m = Matrix::zeros(n);
        for p in self.ty.mass.points() {
            let g = self.gradients(p.xi)?;
            let w = density * g.det * p.weight;
            let v = g.shape.values();
            for (i, a) in v.iter().enumerate() {
                for (j, b) in v.iter().enumerate() {
                    m.add(i, j, w * a * b);
                }
            }
        }
        Ok(m)
    }

    /// The lumped (diagonal) mass per node: the consistent diagonal scaled to the element's mass.
    pub fn lumped_mass(&self, density: f64) -> Result<Vec<f64>> {
        lump(&self.mass(density)?).ok_or_else(|| self.degenerate())
    }

    /// Consistent nodal forces of a uniform force per volume: fᵢ = ∫ Nᵢ b dV.
    pub fn body_force(&self, b: [f64; 3]) -> Result<Vec<[f64; 3]>> {
        let mut f = vec![[0.0; 3]; self.ty.nodes];
        for p in self.ty.mass.points() {
            let g = self.gradients(p.xi)?;
            for (fi, n) in f.iter_mut().zip(g.shape.values()) {
                for k in 0..3 {
                    fi[k] += n * b[k] * g.det * p.weight;
                }
            }
        }
        Ok(f)
    }

    /// Consistent nodal forces of a force per volume that varies with position (such as a
    /// centrifugal load): fᵢ = ∫ Nᵢ b(x) dV.
    pub fn body_force_field(&self, b: impl Fn([f64; 3]) -> [f64; 3]) -> Result<Vec<[f64; 3]>> {
        let mut f = vec![[0.0; 3]; self.ty.nodes];
        for p in self.ty.mass.points() {
            let g = self.gradients(p.xi)?;
            let b = b(self.position(p.xi));
            for (fi, n) in f.iter_mut().zip(g.shape.values()) {
                for k in 0..3 {
                    fi[k] += n * b[k] * g.det * p.weight;
                }
            }
        }
        Ok(f)
    }

    /// The thermal strain α ΔT at a point with shape values `n` (xx, yy, zz; no shear).
    pub(crate) fn thermal_strain(thermal: &ThermalStrain, n: &[f64]) -> [f64; 6] {
        let dt: f64 = n.iter().zip(&thermal.delta).map(|(a, b)| a * b).sum();
        [thermal.alpha[0] * dt, thermal.alpha[1] * dt, thermal.alpha[2] * dt, 0.0, 0.0, 0.0]
    }

    /// The equivalent nodal forces of a thermal strain, ∫ Bᵀ D ε_th dV, with the stiffness rule
    /// (the initial-strain load of Zienkiewicz, Taylor & Zhu 2013, §6.2.3); incompatible modes
    /// are condensed out as in the stiffness, f = f_u − K_ua K_aa⁻¹ f_a.
    pub fn thermal_load(&self, d: &Elasticity, thermal: &ThermalStrain) -> Result<Vec<[f64; 3]>> {
        thermal.check(self.ty.nodes)?;
        if let Formulation::Hybrid { .. } = self.formulation {
            return self.hybrid_thermal_load(d, thermal);
        }
        if self.formulation == Formulation::IncompatibleModes {
            let IncompatibleParts { kua, kaa, j0, .. } = self.incompatible_parts(d, false)?;
            let (mut f, mut fa) = self.thermal_forces(d, thermal, Some(&j0))?;
            Cholesky::new(&kaa)?.solve(&mut fa);
            for (i, fi) in f.iter_mut().flatten().enumerate() {
                *fi -= (0..INCOMPATIBLE_DOFS).map(|a| kua.get(i * INCOMPATIBLE_DOFS + a).copied().unwrap_or(0.0) * fa[a]).sum::<f64>();
            }
            return Ok(f);
        }
        Ok(self.thermal_forces(d, thermal, None)?.0)
    }

    /// ∫ Bᵀ D ε_th dV on the nodes and, with the centre Jacobian `j0`, on the incompatible modes.
    fn thermal_forces(
        &self,
        d: &Elasticity,
        thermal: &ThermalStrain,
        j0: Option<&([[f64; 3]; 3], f64)>,
    ) -> Result<(Vec<[f64; 3]>, [f64; INCOMPATIBLE_DOFS])> {
        let mut f = vec![[0.0; 3]; self.ty.nodes];
        let mut fa = [0.0; INCOMPATIBLE_DOFS];
        let add = |out: &mut [f64], g: &[f64; 3], sig: &[f64; 6], w: f64| {
            let b = b_block(g);
            for (c, o) in out.iter_mut().enumerate() {
                *o += w * (0..6).map(|s| b[s][c] * sig[s]).sum::<f64>();
            }
        };
        for p in self.rule().points() {
            let g = self.gradients(p.xi)?;
            let w = g.det * p.weight;
            let sig = stress(d, &Self::thermal_strain(thermal, g.shape.values()));
            for (fi, gi) in f.iter_mut().zip(g.dx.get(..self.ty.nodes).unwrap_or_default()) {
                add(fi, gi, &sig, w);
            }
            if let Some(j0) = j0 {
                for (k, gk) in self.incompatible_gradients(p.xi, g.det, j0).iter().enumerate() {
                    add(fa.get_mut(3 * k..3 * k + 3).unwrap_or_default(), gk, &sig, w);
                }
            }
        }
        Ok((f, fa))
    }

    /// Strain and stress at every integration point for the nodal displacements `u`.
    pub fn stresses(&self, d: &Elasticity, u: &[[f64; 3]]) -> Result<Vec<PointStress>> {
        self.stresses_with(d, u, None)
    }

    /// Strain and stress at every integration point for the nodal displacements `u` and, if
    /// given, a thermal strain: the stress is D (ε − ε_th), the strain reported the total ε.
    pub fn stresses_with(&self, d: &Elasticity, u: &[[f64; 3]], thermal: Option<&ThermalStrain>) -> Result<Vec<PointStress>> {
        if u.len() != self.ty.nodes {
            return Err(Error::Invalid(format!("{} displacements for a {}-node element", u.len(), self.ty.nodes)));
        }
        if let Some(t) = thermal {
            t.check(self.ty.nodes)?;
        }
        // Hybrid: the element's pressure field, and σ = D_dev (ε − ε_th) + p m.
        let hybrid = if let Formulation::Hybrid { .. } = self.formulation { Some(self.pressure_field(d, u, thermal)?) } else { None };
        // Incompatible modes: recover the internal amplitudes from K_aa α + K_au u = f_a (the
        // thermal load on the modes; zero without one).
        let alpha = if self.formulation == Formulation::IncompatibleModes {
            let IncompatibleParts { kua, kaa, j0, .. } = self.incompatible_parts(d, false)?;
            let flat: Vec<f64> = u.iter().flatten().copied().collect();
            let mut rhs = match thermal {
                Some(t) => self.thermal_forces(d, t, Some(&j0))?.1,
                None => [0.0; INCOMPATIBLE_DOFS],
            };
            for (i, ui) in flat.iter().enumerate() {
                for (a, r) in rhs.iter_mut().enumerate() {
                    *r -= kua.get(i * INCOMPATIBLE_DOFS + a).copied().unwrap_or(0.0) * ui;
                }
            }
            Cholesky::new(&kaa)?.solve(&mut rhs);
            Some((rhs, j0))
        } else {
            None
        };
        let mut out = Vec::new();
        for p in self.rule().points() {
            let g = self.gradients(p.xi)?;
            let mut grad = [[0.0; 3]; 3]; // ∂u_c/∂x_b
            for (dn, ui) in g.dx.iter().zip(u) {
                for c in 0..3 {
                    for b in 0..3 {
                        grad[c][b] += ui[c] * dn[b];
                    }
                }
            }
            if let Some((a, j0)) = &alpha {
                let ga = self.incompatible_gradients(p.xi, g.det, j0);
                for (k, gk) in ga.iter().enumerate() {
                    for c in 0..3 {
                        for b in 0..3 {
                            grad[c][b] += a[3 * k + c] * gk[b];
                        }
                    }
                }
            }
            let strain = [grad[0][0], grad[1][1], grad[2][2], grad[0][1] + grad[1][0], grad[0][2] + grad[2][0], grad[1][2] + grad[2][1]];
            let xyz = self.position_of(g.shape.values());
            // The mechanical strain, less the thermal strain when there is one.
            let mechanical = match thermal {
                Some(t) => {
                    let eth = Self::thermal_strain(t, g.shape.values());
                    [0, 1, 2, 3, 4, 5].map(|k| strain[k] - eth[k])
                }
                None => strain,
            };
            let sigma = match &hybrid {
                Some(h) => h.stress(&mechanical, xyz),
                None => stress(d, &mechanical),
            };
            out.push(PointStress { xyz, strain, stress: sigma });
        }
        Ok(out)
    }

    /// Extrapolate values at the integration points to the element's nodes.
    pub fn extrapolate(&self, values: &[[f64; 6]]) -> Result<Vec<[f64; 6]>> {
        let e = extrapolation(self.ty, self.rule())?;
        if values.len() != e.points {
            return Err(Error::Invalid(format!("{} values for {} integration points", values.len(), e.points)));
        }
        Ok(e.apply(values))
    }
}

/// Make a matrix exactly symmetric (the mean of it and its transpose), so assembly can rely on it.
/// HRZ lumping of a consistent mass matrix: its diagonal scaled to the element's mass (the
/// entries sum to it, the shape functions summing to one); `None` when the diagonal has no mass.
pub(crate) fn lump(m: &Matrix) -> Option<Vec<f64>> {
    let total: f64 = m.data.iter().sum();
    let diag: Vec<f64> = (0..m.n).map(|i| m.get(i, i)).collect();
    let sum: f64 = diag.iter().sum();
    (sum > 0.0).then(|| diag.iter().map(|d| d * total / sum).collect())
}

pub(crate) fn symmetrize(k: &mut Matrix) {
    let n = k.n;
    for i in 0..n {
        for j in i + 1..n {
            let v = 0.5 * (k.get(i, j) + k.get(j, i));
            if let Some(x) = k.data.get_mut(i * n + j) {
                *x = v;
            }
            if let Some(x) = k.data.get_mut(j * n + i) {
                *x = v;
            }
        }
    }
}

/// The strain-displacement block of one gradient g: ε = B u for a single node's (x, y, z).
pub(crate) fn b_block(g: &[f64; 3]) -> [[f64; 3]; 6] {
    [[g[0], 0.0, 0.0], [0.0, g[1], 0.0], [0.0, 0.0, g[2]], [g[1], g[0], 0.0], [g[2], 0.0, g[0]], [0.0, g[2], g[1]]]
}

/// out[3i + r][3j + c] += w (Bᵢᵀ D Bⱼ)[r][c] for gradients `ga` (rows) and `gb` (columns), with
/// `cols` columns in the row-major `out`.
fn accumulate(out: &mut [f64], cols: usize, ga: &[[f64; 3]], gb: &[[f64; 3]], d: &Elasticity, w: f64) {
    for (j, gj) in gb.iter().enumerate() {
        let bj = b_block(gj);
        // D Bⱼ: 6 × 3.
        let mut db = [[0.0; 3]; 6];
        for s in 0..6 {
            for c in 0..3 {
                db[s][c] = (0..6).map(|t| d[s][t] * bj[t][c]).sum();
            }
        }
        for (i, gi) in ga.iter().enumerate() {
            let bi = b_block(gi);
            for r in 0..3 {
                for c in 0..3 {
                    let v: f64 = (0..6).map(|s| bi[s][r] * db[s][c]).sum();
                    if let Some(x) = out.get_mut((3 * i + r) * cols + 3 * j + c) {
                        *x += w * v;
                    }
                }
            }
        }
    }
}

/// The symmetric case of [`accumulate`] (rows and columns from the same gradients, `out` square):
/// computes the blocks with j ≥ i and mirrors them, so the result is exactly symmetric.
pub(crate) fn accumulate_symmetric(out: &mut [f64], g: &[[f64; 3]], d: &Elasticity, w: f64) {
    let cols = 3 * g.len();
    for (j, gj) in g.iter().enumerate() {
        let bj = b_block(gj);
        let mut db = [[0.0; 3]; 6];
        for s in 0..6 {
            for c in 0..3 {
                db[s][c] = (0..6).map(|t| d[s][t] * bj[t][c]).sum();
            }
        }
        for (i, gi) in g.iter().enumerate().take(j + 1) {
            let bi = b_block(gi);
            for r in 0..3 {
                // On a diagonal block only the upper triangle is computed, then mirrored.
                for c in if i == j { r } else { 0 }..3 {
                    let v = w * (0..6).map(|s| bi[s][r] * db[s][c]).sum::<f64>();
                    let (row, col) = (3 * i + r, 3 * j + c);
                    if let Some(x) = out.get_mut(row * cols + col) {
                        *x += v;
                    }
                    if row != col
                        && let Some(x) = out.get_mut(col * cols + row)
                    {
                        *x += v;
                    }
                }
            }
        }
    }
}

/// The map from integration-point values to nodal values of an element type and rule.
#[derive(Debug, Clone)]
pub struct Extrapolation {
    pub points: usize,
    /// nodes × points, row-major.
    pub weights: Vec<f64>,
}

impl Extrapolation {
    pub fn apply(&self, values: &[[f64; 6]]) -> Vec<[f64; 6]> {
        self.weights
            .chunks_exact(self.points.max(1))
            .map(|row| {
                let mut v = [0.0; 6];
                for (w, x) in row.iter().zip(values) {
                    for k in 0..6 {
                        v[k] += w * x[k];
                    }
                }
                v
            })
            .collect()
    }
}

/// The least-squares extrapolation of a rule's points to the nodes with the linear corner
/// functions (Hinton & Campbell 1974); the mean when the rule has fewer points than corners.
pub fn extrapolation(ty: &ElementType, rule: Rule) -> Result<Extrapolation> {
    let pts = rule.points();
    let g = pts.len();
    let c = ty.corners;
    let corner_basis = |xi: [f64; 3]| match ty.family {
        Family::Tetrahedron => shape::tet4(xi),
        Family::Wedge => shape::wedge6(xi),
        Family::Pyramid => shape::pyramid5(xi),
        Family::Hexahedron => shape::hex8(xi),
        Family::Triangle => Basis::Tri3.eval(xi),
        Family::Quadrilateral => Basis::Quad4.eval(xi),
    };
    let mut corner_rows = vec![0.0; c * g];
    if g < c {
        corner_rows.iter_mut().for_each(|w| *w = 1.0 / g.max(1) as f64);
    } else {
        // (AᵀA)⁻¹ Aᵀ with A[p][j] = corner function j at point p.
        let a: Vec<Vec<f64>> = pts.iter().map(|p| corner_basis(p.xi).values().to_vec()).collect();
        let mut ata = Matrix::zeros(c);
        for row in &a {
            for i in 0..c {
                for j in 0..c {
                    ata.add(i, j, row.get(i).copied().unwrap_or(0.0) * row.get(j).copied().unwrap_or(0.0));
                }
            }
        }
        let chol = Cholesky::new(&ata)?;
        for (p, row) in a.iter().enumerate() {
            let mut col = row.clone();
            chol.solve(&mut col);
            for (i, v) in col.iter().enumerate() {
                if let Some(x) = corner_rows.get_mut(i * g + p) {
                    *x = *v;
                }
            }
        }
    }
    let mut weights = vec![0.0; ty.nodes * g];
    for (w, v) in weights.iter_mut().zip(&corner_rows) {
        *w = *v;
    }
    for edge in ty.kind.edges() {
        if let [a, b, m] = **edge {
            for p in 0..g {
                let v = 0.5 * (corner_rows.get(a * g + p).copied().unwrap_or(0.0) + corner_rows.get(b * g + p).copied().unwrap_or(0.0));
                if let Some(x) = weights.get_mut(m * g + p) {
                    *x = v;
                }
            }
        }
    }
    Ok(Extrapolation { points: g, weights })
}
