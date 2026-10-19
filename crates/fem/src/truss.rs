//! Truss (axial) elements: the two-node rod2 and the three-node rod3, straight or curved, in any
//! orientation (a plane truss is a space truss whose nodes lie in a plane).
//!
//! **Formulation.** The isoparametric bar of R. D. Cook, D. S. Malkus, M. E. Plesha, R. J. Witt,
//! *Concepts and Applications of Finite Element Analysis*, 4th ed., Wiley 2002, §3.6 and §6.2
//! (and T. J. R. Hughes, *The Finite Element Method*, Dover 2000, §1.10), written along the
//! member's arc in three dimensions. On ξ ∈ [−1, 1] the shape functions are
//! N = ((1 − ξ)/2, (1 + ξ)/2) for rod2 and N = (ξ(ξ − 1)/2, ξ(ξ + 1)/2, 1 − ξ²) for rod3 (the
//! ends, then the middle node). The position is x(ξ) = Σ Nᵢ xᵢ, the arc length element
//! ds = J dξ with J = |x′(ξ)|, and the unit tangent t = x′/J. The axial strain is the stretch
//! along the tangent, ε = t · ∂u/∂s = Σᵢ (Nᵢ′/J) t · uᵢ, so with Bᵢ = (Nᵢ′/J) tᵀ:
//!
//! - stiffness K = ∫ E A Bᵀ B ds,
//! - consistent mass ρ A ∫ Nᵢ Nⱼ ds per direction, lumped mass by HRZ diagonal scaling (E. Hinton,
//!   T. Rock, O. C. Zienkiewicz, *A note on mass lumping and related processes in the finite
//!   element method*, Earthq. Eng. Struct. Dyn. 4(3):245–249, 1976), which gives ρAL/2 at each
//!   end of rod2 and ρAL·(1/6, 1/6, 2/3) on a straight rod3,
//! - body force A ∫ Nᵢ b ds for a force b per volume,
//! - stress σ = E ε along the tangent: in global axes the uniaxial tensor σ t tᵀ (strain ε t tᵀ,
//!   with engineering shear strains), so stresses of trusses and solids read the same way and
//!   the axial force is N = A tᵀ σ t.
//!
//! The rules are Gauss–Legendre with two points for rod2 and three for rod3: exact for the
//! stiffness and mass of a straight member, and the usual full integration of a curved one.
//! Material: the Young's modulus of an isotropic material; the truss knows no other constant.

use fieldcraft_model::{Elastic, ElementKind};

use crate::basis::Basis;
use crate::dense::Matrix;
use crate::dof::{DofSet, ElementDofs};
use crate::element::{PointStress, ThermalStrain};
use crate::quadrature::gauss_legendre;
use crate::{Error, Result};

/// What a truss's nodes carry: the three translations (a plane truss is a space truss whose
/// out-of-plane motion is held; trusses of planar models will declare [`DofSet::PLANE`]).
pub const DOFS: ElementDofs = ElementDofs::new(DofSet::TRANSLATIONS);

/// A member whose tangent is shorter than this fraction of its chord somewhere has folded back
/// or has no length.
const FOLDED: f64 = 1e-9; // unit-factors: relative to the member's chord

/// The Young's modulus a truss uses, or why there is none.
pub fn young_modulus(elastic: &Elastic) -> Result<f64> {
    match *elastic {
        Elastic::Isotropic { young, .. } if young > 0.0 && young.is_finite() => Ok(young),
        Elastic::Isotropic { young, .. } => Err(Error::Invalid(format!("Young's modulus must be positive, got {young}"))),
        Elastic::EngineeringConstants { .. } => {
            Err(Error::NotSupported("orthotropic materials on trusses (a truss takes one Young's modulus: use an isotropic material)".into()))
        }
    }
}

/// Shape functions and their ξ-derivatives ([`Basis::Line2`], [`Basis::Line3`]).
fn shape(n: usize, xi: f64) -> ([f64; 3], [f64; 3]) {
    let s = if n == 2 { Basis::Line2 } else { Basis::Line3 }.eval([xi, 0.0, 0.0]);
    let (mut v, mut d) = ([0.0; 3], [0.0; 3]);
    for i in 0..n.min(3) {
        v[i] = s.n[i];
        d[i] = s.dn[i][0];
    }
    (v, d)
}

/// The geometry at one point: shape functions, their arc-length derivatives Nᵢ′/J, the unit
/// tangent and the arc length per unit ξ.
#[derive(Debug, Clone)]
struct Point {
    n: [f64; 3],
    ds: [f64; 3],
    tangent: [f64; 3],
    jacobian: f64,
    weight: f64,
    xyz: [f64; 3],
}

/// A truss element on its node coordinates (model node order), with the geometry at its
/// integration points worked out once.
#[derive(Debug, Clone)]
pub struct Truss<'a> {
    pub kind: ElementKind,
    xyz: &'a [[f64; 3]],
    rule: Vec<Point>,
}

impl<'a> Truss<'a> {
    pub fn new(kind: ElementKind, xyz: &'a [[f64; 3]]) -> Result<Truss<'a>> {
        if !matches!(kind, ElementKind::Rod2 | ElementKind::Rod3) {
            return Err(Error::Invalid(format!("{} is not a truss element", kind.name())));
        }
        if xyz.len() != kind.node_count() {
            return Err(Error::Invalid(format!("{} needs {} nodes, got {}", kind.name(), kind.node_count(), xyz.len())));
        }
        // The tangent must not vanish or turn back along the member (a middle node off its span).
        let scale = chord(xyz).max(f64::MIN_POSITIVE);
        let reference = at(xyz, 0.0, 0.0)?.tangent;
        for xi in [-1.0, -0.5, 0.5, 1.0] {
            let p = at(xyz, xi, 0.0)?;
            if !(p.jacobian > FOLDED * scale) || p.tangent.iter().zip(&reference).map(|(a, b)| a * b).sum::<f64>() <= 0.0 {
                return Err(Error::Degenerate("the member folds back on itself or has no length".into()));
            }
        }
        let rule = gauss_legendre(kind.node_count()).iter().map(|(xi, w)| at(xyz, *xi, *w)).collect::<Result<_>>()?;
        Ok(Truss { kind, xyz, rule })
    }

    fn nodes(&self) -> usize {
        self.kind.node_count()
    }

    /// The integration points: two for rod2, three for rod3.
    fn points(&self) -> Result<&[Point]> {
        Ok(&self.rule)
    }

    /// The natural coordinates of the integration points, in order.
    pub fn abscissae(&self) -> Vec<f64> {
        gauss_legendre(self.nodes()).iter().map(|(xi, _)| *xi).collect()
    }

    /// The arc length.
    pub fn length(&self) -> Result<f64> {
        Ok(self.points()?.iter().map(|p| p.weight * p.jacobian).sum())
    }

    /// The stiffness (3n × 3n, translations node by node) for Young's modulus `young` and
    /// cross-sectional area `area`.
    pub fn stiffness(&self, young: f64, area: f64) -> Result<Matrix> {
        let n = self.nodes();
        let mut k = Matrix::zeros(3 * n);
        for p in self.points()? {
            let b: Vec<f64> = (0..3 * n).map(|a| p.ds[a / 3] * p.tangent[a % 3]).collect();
            let factor = young * area * p.weight * p.jacobian;
            for (i, bi) in b.iter().enumerate() {
                for (j, bj) in b.iter().enumerate() {
                    k.add(i, j, factor * bi * bj);
                }
            }
        }
        Ok(k)
    }

    /// The consistent mass per direction (n × n).
    pub fn mass(&self, density: f64, area: f64) -> Result<Matrix> {
        let n = self.nodes();
        let mut m = Matrix::zeros(n);
        for p in self.points()? {
            let factor = density * area * p.weight * p.jacobian;
            for i in 0..n {
                for j in 0..n {
                    m.add(i, j, factor * p.n[i] * p.n[j]);
                }
            }
        }
        Ok(m)
    }

    /// The HRZ lumped mass per node: the consistent diagonal scaled to the member's mass.
    pub fn lumped_mass(&self, density: f64, area: f64) -> Result<Vec<f64>> {
        let m = self.mass(density, area)?;
        let n = self.nodes();
        let diagonal: Vec<f64> = (0..n).map(|i| m.get(i, i)).collect();
        let total = density * area * self.length()?;
        let sum: f64 = diagonal.iter().sum();
        if !(sum > 0.0) {
            return Err(Error::Degenerate("the member has no mass to lump".into()));
        }
        Ok(diagonal.iter().map(|d| d * total / sum).collect())
    }

    /// Consistent nodal forces of a force `b` per volume.
    pub fn body_force(&self, b: [f64; 3], area: f64) -> Result<Vec<[f64; 3]>> {
        let mut f = vec![[0.0; 3]; self.nodes()];
        for p in self.points()? {
            for (i, fi) in f.iter_mut().enumerate() {
                let w = area * p.weight * p.jacobian * p.n[i];
                for k in 0..3 {
                    fi[k] += w * b[k];
                }
            }
        }
        Ok(f)
    }

    /// At each integration point, the row that gives the axial strain from the element's
    /// displacements (3n, node by node), and the tangent there: ε = row · u.
    pub fn strain_operator(&self) -> Result<Vec<(Vec<f64>, [f64; 3])>> {
        let n = self.nodes();
        Ok(self.points()?.iter().map(|p| ((0..3 * n).map(|a| p.ds[a / 3] * p.tangent[a % 3]).collect(), p.tangent)).collect())
    }

    /// The axial strain at each integration point for nodal displacements `u`.
    pub fn axial_strains(&self, u: &[[f64; 3]]) -> Result<Vec<f64>> {
        if u.len() != self.nodes() {
            return Err(Error::Invalid(format!("{} displacements for a {}", u.len(), self.kind.name())));
        }
        Ok(self.points()?.iter().map(|p| (0..self.nodes()).map(|i| p.ds[i] * (0..3).map(|k| p.tangent[k] * u[i][k]).sum::<f64>()).sum()).collect())
    }

    /// Strain and stress at the integration points as global tensors (see the module docs).
    pub fn stresses(&self, young: f64, u: &[[f64; 3]]) -> Result<Vec<PointStress>> {
        self.stresses_with(young, u, None)
    }

    /// The axial thermal strain at each integration point: the thermal strain tensor
    /// diag(α) ΔT along the tangent, Σ αₖ tₖ² ΔT, with ΔT interpolated from the nodes.
    fn thermal_strains(&self, thermal: &ThermalStrain) -> Result<Vec<f64>> {
        thermal.check(self.nodes())?;
        Ok(self
            .points()?
            .iter()
            .map(|p| {
                let dt: f64 = p.n.iter().zip(&thermal.delta).map(|(a, b)| a * b).sum();
                (0..3).map(|k| thermal.alpha[k] * p.tangent[k] * p.tangent[k]).sum::<f64>() * dt
            })
            .collect())
    }

    /// The equivalent nodal forces of a thermal strain, ∫ Bᵀ E A ε_th ds.
    pub fn thermal_load(&self, young: f64, area: f64, thermal: &ThermalStrain) -> Result<Vec<[f64; 3]>> {
        let strains = self.thermal_strains(thermal)?;
        let mut f = vec![[0.0; 3]; self.nodes()];
        for (p, e) in self.points()?.iter().zip(strains) {
            let force = young * area * e * p.weight * p.jacobian;
            for (fi, ds) in f.iter_mut().zip(p.ds) {
                for k in 0..3 {
                    fi[k] += force * ds * p.tangent[k];
                }
            }
        }
        Ok(f)
    }

    /// Strain and stress at the integration points with an optional thermal strain: the
    /// stress is E (ε − ε_th), the strain reported the total.
    pub fn stresses_with(&self, young: f64, u: &[[f64; 3]], thermal: Option<&ThermalStrain>) -> Result<Vec<PointStress>> {
        let strains = self.axial_strains(u)?;
        let thermal = match thermal {
            Some(t) => self.thermal_strains(t)?,
            None => vec![0.0; strains.len()],
        };
        Ok(self
            .points()?
            .iter()
            .zip(strains.into_iter().zip(thermal))
            .map(|(p, (e, th))| PointStress { xyz: p.xyz, strain: uniaxial(e, p.tangent, 2.0), stress: uniaxial(young * (e - th), p.tangent, 1.0) })
            .collect())
    }

    /// Integration-point values extrapolated to the nodes (the polynomial through the points,
    /// evaluated at the nodes' ξ: −1, 1, then 0).
    pub fn extrapolate(&self, values: &[[f64; 6]]) -> Result<Vec<[f64; 6]>> {
        let xs = self.abscissae();
        if values.len() != xs.len() {
            return Err(Error::Invalid(format!("{} values for {} integration points", values.len(), xs.len())));
        }
        let at = |xi: f64| {
            let mut out = [0.0; 6];
            for (a, va) in xs.iter().zip(values) {
                let l: f64 = xs.iter().filter(|b| *b != a).map(|b| (xi - b) / (a - b)).product();
                for k in 0..6 {
                    out[k] += l * va[k];
                }
            }
            out
        };
        Ok([-1.0, 1.0, 0.0].into_iter().take(self.nodes()).map(at).collect())
    }

    /// The tangent at each integration point.
    pub fn tangents(&self) -> Result<Vec<[f64; 3]>> {
        Ok(self.points()?.iter().map(|p| p.tangent).collect())
    }

    /// An upper bound of the highest natural frequency squared of the element with its lumped
    /// mass, for explicit time integration: exact (4c²/L²) for rod2, and Gershgorin's bound on
    /// M^{-1/2} K M^{-1/2} for rod3.
    pub fn frequency_bound(&self, young: f64, density: f64) -> Result<f64> {
        if self.kind == ElementKind::Rod2 {
            let l = chord(self.xyz);
            return Ok(4.0 * young / (density * l * l));
        }
        // The area cancels: K ∝ A and M ∝ A.
        let k = self.stiffness(young, 1.0)?;
        let m = self.lumped_mass(density, 1.0)?;
        let mut worst: f64 = 0.0;
        for i in 0..k.n {
            let mi = m.get(i / 3).copied().unwrap_or(0.0);
            let row: f64 = (0..k.n).map(|j| k.get(i, j).abs() / (mi * m.get(j / 3).copied().unwrap_or(0.0)).sqrt()).sum();
            worst = worst.max(row);
        }
        Ok(worst)
    }
}

/// The distance between the ends.
fn chord(xyz: &[[f64; 3]]) -> f64 {
    match (xyz.first(), xyz.get(1)) {
        (Some(a), Some(b)) => ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt(),
        _ => 0.0,
    }
}

/// The geometry at `xi` of the element on `xyz` (two or three nodes).
fn at(xyz: &[[f64; 3]], xi: f64, weight: f64) -> Result<Point> {
    let (n, dn) = shape(xyz.len(), xi);
    let mut d = [0.0; 3];
    let mut x = [0.0; 3];
    for (i, p) in xyz.iter().enumerate() {
        for k in 0..3 {
            d[k] += dn[i] * p[k];
            x[k] += n[i] * p[k];
        }
    }
    let jacobian = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if !(jacobian > 0.0 && jacobian.is_finite()) {
        return Err(Error::Degenerate("the member has no length".into()));
    }
    Ok(Point { n, ds: dn.map(|v| v / jacobian), tangent: d.map(|c| c / jacobian), jacobian, weight, xyz: x })
}

/// The tensor s t tᵀ as a 6-vector (xx, yy, zz, xy, xz, yz), shear terms times `shear` (2 for
/// engineering shear strains).
pub fn uniaxial(s: f64, t: [f64; 3], shear: f64) -> [f64; 6] {
    [s * t[0] * t[0], s * t[1] * t[1], s * t[2] * t[2], shear * s * t[0] * t[1], shear * s * t[0] * t[2], shear * s * t[1] * t[2]]
}

/// The axial value (stress or strain) of a uniaxial tensor along a unit tangent: tᵀ σ t, with
/// engineering shear strains halved when `engineering`.
pub fn axial(v: &[f64; 6], t: [f64; 3], engineering: bool) -> f64 {
    let h = if engineering { 0.5 } else { 1.0 };
    v[0] * t[0] * t[0] + v[1] * t[1] * t[1] + v[2] * t[2] * t[2] + 2.0 * h * (v[3] * t[0] * t[1] + v[4] * t[0] * t[2] + v[5] * t[1] * t[2])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn norm(a: [f64; 3]) -> f64 {
        (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
    }

    /// A rod along a skew direction, and a rod3 with its middle node at the midpoint.
    fn straight(kind: ElementKind) -> Vec<[f64; 3]> {
        let (a, b) = ([0.5, -1.0, 2.0], [2.5, 1.0, 3.0]);
        let mut x = vec![a, b];
        if kind == ElementKind::Rod3 {
            x.push([0, 1, 2].map(|k| (a[k] + b[k]) / 2.0));
        }
        x
    }

    #[test]
    fn straight_rods_match_the_textbook_matrices() {
        let (young, area, density) = (200.0, 3.0, 2.0);
        for kind in [ElementKind::Rod2, ElementKind::Rod3] {
            let x = straight(kind);
            let t = Truss::new(kind, &x).unwrap();
            let l = norm([2.0, 2.0, 1.0]);
            assert!((t.length().unwrap() - l).abs() < 1e-14);
            // Rigid translations and the rotation about the axis' normals: zero force.
            let k = t.stiffness(young, area).unwrap();
            for d in 0..3 {
                let u: Vec<f64> = (0..k.n).map(|a| if a % 3 == d { 1.0 } else { 0.0 }).collect();
                assert!(k.mul_vec(&u).iter().all(|f| f.abs() < 1e-12), "{kind:?}");
            }
            // A uniform stretch ε along the axis: strain ε at every point, end forces ±EAε.
            let eps = 1e-3;
            let dir = [2.0 / l, 2.0 / l, 1.0 / l];
            let u: Vec<[f64; 3]> =
                x.iter().map(|p| [0, 1, 2].map(|c| eps * dir[c] * ((p[0] - 0.5) * dir[0] + (p[1] + 1.0) * dir[1] + (p[2] - 2.0) * dir[2]))).collect();
            for e in t.axial_strains(&u).unwrap() {
                assert!((e - eps).abs() < 1e-15, "{kind:?} {e}");
            }
            let flat: Vec<f64> = u.iter().flatten().copied().collect();
            let f = k.mul_vec(&flat);
            for c in 0..3 {
                assert!((f[3 + c] - young * area * eps * dir[c]).abs() < 1e-12, "{kind:?}");
            }
            let s = t.stresses(young, &u).unwrap();
            for p in &s {
                assert!((axial(&p.stress, dir, false) - young * eps).abs() < 1e-12);
                assert!((axial(&p.strain, dir, true) - eps).abs() < 1e-15);
            }
            let nodal = t.extrapolate(&s.iter().map(|p| p.stress).collect::<Vec<_>>()).unwrap();
            assert!(nodal.iter().all(|v| (axial(v, dir, false) - young * eps).abs() < 1e-12));
            // Mass: total ρAL per direction; lumped HRZ.
            let m = t.mass(density, area).unwrap();
            let total = density * area * l;
            assert!((m.data.iter().sum::<f64>() - total).abs() < 1e-12);
            let lumped = t.lumped_mass(density, area).unwrap();
            let want: Vec<f64> = if kind == ElementKind::Rod2 { vec![0.5, 0.5] } else { vec![1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0] };
            for (a, b) in lumped.iter().zip(want) {
                assert!((a - b * total).abs() < 1e-12, "{kind:?} {lumped:?}");
            }
            let fb = t.body_force([0.0, 0.0, -1.0], area).unwrap();
            assert!((fb.iter().map(|f| f[2]).sum::<f64>() + area * l).abs() < 1e-12);
            // The highest frequency of the free element with its lumped mass: exact for rod2,
            // an upper bound of the exact one for rod3.
            let bound = t.frequency_bound(young, density).unwrap();
            let c2 = young / density;
            if kind == ElementKind::Rod2 {
                assert!((bound - 4.0 * c2 / (l * l)).abs() < 1e-9 * bound);
            } else {
                // Exact: axial eigenproblem of K = EA/(3L)[7 1 −8; 1 7 −8; −8 −8 16] with
                // M = ρAL diag(1/6, 1/6, 2/3): ω²max = 24 c²/L² (antisymmetric 12, symmetric 24).
                assert!(bound >= 24.0 * c2 / (l * l) * (1.0 - 1e-12), "{bound} vs {}", 24.0 * c2 / (l * l));
            }
        }
    }

    #[test]
    fn bad_geometry_is_an_error() {
        assert!(Truss::new(ElementKind::Rod2, &[[1.0; 3], [1.0; 3]]).is_err());
        // The middle node beyond an end folds the member back.
        assert!(Truss::new(ElementKind::Rod3, &[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]).is_err());
        assert!(Truss::new(ElementKind::Rod3, &[[0.0; 3], [1.0, 0.0, 0.0]]).is_err());
        assert!(Truss::new(ElementKind::Hex8, &[[0.0; 3]; 8]).is_err());
        assert!(Truss::new(ElementKind::Rod2, &[[0.0; 3], [f64::NAN, 0.0, 0.0]]).is_err());
        assert!(young_modulus(&Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0; 3] }).is_err());
    }

    /// A curved member: a circular arc of radius R (middle node on the arc). Any linear
    /// displacement field is reproduced exactly: a uniform expansion u = s x (a ring growing
    /// radially by δ = s R) strains it by s everywhere, an infinitesimal rotation not at all.
    #[test]
    fn curved_rod3_reproduces_linear_fields() {
        let (radius, young) = (3.0, 70.0);
        for span in [0.3f64, 1.0, 2.0] {
            let at = |a: f64| [radius * a.cos(), radius * a.sin(), 0.5];
            let x = [at(0.0), at(span), at(span / 2.0)];
            let t = Truss::new(ElementKind::Rod3, &x).unwrap();
            let (s, w) = (2e-3, [1e-3, -2e-3, 3e-3]);
            let u: Vec<[f64; 3]> =
                x.iter().map(|p| [0, 1, 2].map(|k| s * p[k] + w[(k + 1) % 3] * p[(k + 2) % 3] - w[(k + 2) % 3] * p[(k + 1) % 3])).collect();
            let points = t.stresses(young, &u).unwrap();
            for (p, tangent) in points.iter().zip(t.tangents().unwrap()) {
                assert!((axial(&p.strain, tangent, true) - s).abs() < 1e-15, "{span}");
                assert!((axial(&p.stress, tangent, false) - young * s).abs() < 1e-13, "{span}");
            }
            // The arc length: exact for a quadratic arc only as the span shrinks.
            let err = (t.length().unwrap() - radius * span).abs() / (radius * span);
            assert!(err < 0.01 * span.powi(4), "{span}: {err}");
        }
    }

    #[test]
    fn thermal_strain_follows_the_axis() {
        // A skew rod3 with orthotropic α: the free expansion along the member, u = ε_t (x − x₀)
        // with ε_t = Σ αₖ tₖ² ΔT, is stress-free and balances the thermal load (K u = f_th).
        let t = [2.0f64, 1.0, 2.0].map(|c| c / 3.0);
        let xyz = [[0.0; 3], t.map(|c| 3.0 * c), t.map(|c| 1.5 * c)];
        let truss = Truss::new(ElementKind::Rod3, &xyz).unwrap();
        let alpha = [1.0e-5, 2.0e-5, 3.0e-5];
        let th = ThermalStrain { alpha, delta: vec![40.0; 3] };
        let eps: f64 = (0..3).map(|k| alpha[k] * t[k] * t[k]).sum::<f64>() * 40.0;
        let u: Vec<[f64; 3]> = xyz.iter().map(|x| x.map(|c| eps * c)).collect();
        for p in truss.stresses_with(2e5, &u, Some(&th)).unwrap() {
            assert!(p.stress.iter().all(|s| s.abs() < 1e-9), "{:?}", p.stress);
        }
        let f = truss.thermal_load(2e5, 2.0, &th).unwrap();
        let ku = truss.stiffness(2e5, 2.0).unwrap().mul_vec(&u.iter().flatten().copied().collect::<Vec<_>>());
        for (a, b) in f.iter().flatten().zip(&ku) {
            assert!((a - b).abs() < 1e-9 * 2e5 * 2.0 * eps, "{a} vs {b}");
        }
        assert!(truss.thermal_load(2e5, 2.0, &ThermalStrain { alpha, delta: vec![1.0] }).is_err());
    }
}
