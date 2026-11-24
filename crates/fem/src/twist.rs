//! Axisymmetric elements with twist: the meridian element of [`crate::plane`] (u_r along X, u_z
//! along Y) with a third unknown per node, the hoop displacement u_θ (the model's Z), interpolated
//! by the same shape functions. It gives the torsion of bodies of revolution.
//!
//! Sources: O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu, *The Finite Element Method: Its Basis and
//! Fundamentals*, 7th ed., 2013, §6.6; S. P. Timoshenko, J. N. Goodier, *Theory of Elasticity*,
//! 3rd ed., 1970, ch. 10 (torsion of circular shafts of variable diameter).

use crate::dense::Matrix;
use crate::dof::{DofSet, ElementDofs};
use crate::element::{PointStress, ThermalStrain};
use crate::material::Elasticity;
use crate::plane::{At, Plane, Theory};
use crate::shape::MAX_NODES;
use crate::{Error, Result};

/// The components of the 6-vector (xx, yy, zz, xy, xz, yz) the twist strains: xz is rθ, yz zθ.
const TWIST: [usize; 2] = [4, 5];

/// The components the meridian element strains: rr, axial, hoop and r–axial.
const MERIDIAN: [usize; 4] = [0, 1, 2, 3];

/// An axisymmetric element with twist: the meridian element and the hoop displacement.
#[derive(Debug, Clone, Copy)]
pub struct Twisted<'a> {
    plane: Plane<'a>,
}

impl<'a> Twisted<'a> {
    /// The element with twist on an axisymmetric meridian element.
    pub fn new(plane: Plane<'a>) -> Result<Twisted<'a>> {
        if plane.theory != Theory::Axisymmetric {
            return Err(Error::Invalid("twist is for axisymmetric elements, not plane stress or plane strain ones".into()));
        }
        Ok(Twisted { plane })
    }

    /// The meridian element: its mass, loads, integration rule and type are the element's.
    pub fn plane(&self) -> &Plane<'a> {
        &self.plane
    }

    /// The element's unknowns: u_r, u_z and u_θ at every node ([`DofSet::TRANSLATIONS`]), and
    /// what the meridian element condenses out.
    pub fn dofs(&self) -> ElementDofs {
        ElementDofs { node: DofSet::TRANSLATIONS, internal: self.plane.dofs().internal }
    }

    /// The elasticity of the twist shears (rθ, zθ), refusing one that couples them to the
    /// meridian strains.
    fn shear(d: &Elasticity) -> Result<[[f64; 2]; 2]> {
        if MERIDIAN.iter().any(|&i| TWIST.iter().any(|&j| d[i][j] != 0.0 || d[j][i] != 0.0)) {
            return Err(Error::NotSupported(
                "a material whose elasticity couples the twist shears to the meridian strains, in axisymmetric elements with twist".into(),
            ));
        }
        Ok(TWIST.map(|i| TWIST.map(|j| d[i][j])))
    }

    /// The twist columns (∂N/∂r − N/r, ∂N/∂z) of the nodes at one integration point.
    fn columns(a: &At) -> [[f64; 2]; MAX_NODES] {
        a.g.map(|g| [g[0] - g[2], g[1]])
    }

    /// The stiffness of the hoop displacement, n × n: ∫ Bᵀ C B dV.
    pub fn twist_stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        let c = Self::shear(d)?;
        let n = self.plane.ty.nodes;
        let mut k = Matrix::zeros(n);
        for (_, a) in self.plane.points()? {
            let b = Self::columns(&a);
            for (j, bj) in b.iter().enumerate().take(n) {
                let cb = [c[0][0] * bj[0] + c[0][1] * bj[1], c[1][0] * bj[0] + c[1][1] * bj[1]];
                for (i, bi) in b.iter().enumerate().take(j + 1) {
                    // Computed for i ≤ j and mirrored, so the matrix is exactly symmetric.
                    let v = a.dv * (bi[0] * cb[0] + bi[1] * cb[1]);
                    k.add(i, j, v);
                    if i != j {
                        k.add(j, i, v);
                    }
                }
            }
        }
        Ok(k)
    }

    /// The stiffness matrix, 3n × 3n, degrees of freedom node by node (u_r, u_z, u_θ).
    pub fn stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        let n = self.plane.ty.nodes;
        let twist = self.twist_stiffness(d)?;
        let meridian = self.plane.stiffness(d)?;
        let mut k = Matrix::zeros(3 * n);
        for i in 0..n {
            for j in 0..n {
                for r in 0..2 {
                    for c in 0..2 {
                        k.add(3 * i + r, 3 * j + c, meridian.get(2 * i + r, 2 * j + c));
                    }
                }
                k.add(3 * i + 2, 3 * j + 2, twist.get(i, j));
            }
        }
        Ok(k)
    }

    /// Strain and stress at every integration point for the nodal displacements `u`
    /// (u_r, u_z, u_θ) and, if given, a thermal strain (which has no shear): the meridian
    /// element's, with γ_rθ, γ_zθ and their stresses in components xz and yz.
    pub fn stresses_with(&self, d: &Elasticity, u: &[[f64; 3]], thermal: Option<&ThermalStrain>) -> Result<Vec<PointStress>> {
        let c = Self::shear(d)?;
        let meridian: Vec<[f64; 2]> = u.iter().map(|v| [v[0], v[1]]).collect();
        let mut out = self.plane.stresses_with(d, &meridian, thermal)?;
        for (p, (_, a)) in out.iter_mut().zip(self.plane.points()?) {
            let (mut rt, mut zt) = (0.0, 0.0);
            for (b, v) in Self::columns(&a).iter().zip(u) {
                rt += b[0] * v[2];
                zt += b[1] * v[2];
            }
            p.strain[TWIST[0]] = rt;
            p.strain[TWIST[1]] = zt;
            p.stress[TWIST[0]] = c[0][0] * rt + c[0][1] * zt;
            p.stress[TWIST[1]] = c[1][0] * rt + c[1][1] * zt;
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Formulation;
    use crate::material::{isotropic, orthotropic};
    use fieldcraft_model::ElementKind;

    /// A distorted element of `kind` with corners in x ∈ [1, 3]; mid-side nodes a little off the
    /// chord midpoints, so quadratic elements are curved (the elements of `plane`'s tests).
    fn element(kind: ElementKind) -> Vec<[f64; 3]> {
        let corners: Vec<[f64; 3]> = match kind {
            ElementKind::Tri3 | ElementKind::Tri6 => vec![[1.0, 0.2, 0.0], [3.0, 0.0, 0.0], [1.6, 1.9, 0.0]],
            _ => vec![[1.0, 0.1, 0.0], [3.0, 0.0, 0.0], [2.7, 2.2, 0.0], [1.2, 1.6, 0.0]],
        };
        let mut xyz = corners.clone();
        for (k, e) in kind.edges().iter().enumerate() {
            if let [a, b, _] = **e {
                let bump = 0.05 * (k as f64 + 1.0);
                xyz.push([(corners[a][0] + corners[b][0]) / 2.0 + bump, (corners[a][1] + corners[b][1]) / 2.0 - bump, 0.0]);
            }
        }
        xyz
    }

    fn variants() -> Vec<(ElementKind, Formulation)> {
        vec![
            (ElementKind::Tri3, Formulation::Full),
            (ElementKind::Tri6, Formulation::Full),
            (ElementKind::Quad4, Formulation::Full),
            (ElementKind::Quad4, Formulation::IncompatibleModes),
            (ElementKind::Quad8, Formulation::Full),
            (ElementKind::Quad8, Formulation::Reduced),
        ]
    }

    #[test]
    fn a_rigid_rotation_strains_nothing_and_the_twist_is_uncoupled() {
        let d = isotropic(1000.0, 0.3).unwrap();
        for (kind, f) in variants() {
            let xyz = element(kind);
            let plane = Plane::new(kind, &xyz, f, Theory::Axisymmetric).unwrap();
            let e = Twisted::new(plane).unwrap();
            assert_eq!(e.dofs().node, DofSet::TRANSLATIONS);
            assert_eq!(e.dofs().internal, plane.dofs().internal);
            let n = xyz.len();
            let k = e.stiffness(&d).unwrap();
            assert_eq!(k.n, 3 * n);
            let meridian = plane.stiffness(&d).unwrap();
            let twist = e.twist_stiffness(&d).unwrap();
            for i in 0..3 * n {
                for j in 0..3 * n {
                    assert_eq!(k.get(i, j), k.get(j, i), "{kind:?} {f:?}");
                    match (i % 3 == 2, j % 3 == 2) {
                        (true, true) => assert_eq!(k.get(i, j), twist.get(i / 3, j / 3)),
                        (false, false) => assert_eq!(k.get(i, j), meridian.get(2 * (i / 3) + i % 3, 2 * (j / 3) + j % 3)),
                        _ => assert_eq!(k.get(i, j), 0.0),
                    }
                }
            }
            let scale = k.data.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            // The rotation about the axis, u_θ = r, and the translation along it.
            let rotation: Vec<f64> = xyz.iter().flat_map(|p| [0.0, 0.0, p[0]]).collect();
            let axial: Vec<f64> = xyz.iter().flat_map(|_| [0.0, 1.0, 0.0]).collect();
            for u in [rotation, axial] {
                let f = k.mul_vec(&u);
                assert!(f.iter().all(|v| v.abs() < 1e-12 * scale), "{kind:?} {f:?}");
            }
            // A uniform hoop displacement is not rigid: γ_rθ = −u_θ / r.
            let uniform: Vec<f64> = xyz.iter().flat_map(|_| [0.0, 0.0, 1.0]).collect();
            let energy: f64 = k.mul_vec(&uniform).iter().zip(&uniform).map(|(a, b)| a * b).sum();
            assert!(energy > 1e-3 * scale, "{kind:?}: {energy}");
        }
    }

    #[test]
    fn a_linear_hoop_displacement_gives_the_exact_shears_on_any_element() {
        // u_θ = a x + b y: γ_rθ = ∂u_θ/∂r − u_θ/r = −b y / x, γ_zθ = b.
        let (young, nu) = (1000.0, 0.3);
        let shear = young / (2.0 * (1.0 + nu));
        let d = isotropic(young, nu).unwrap();
        let (a, b) = (0.003, -0.002);
        // The hybrid meridian element too: the twist is pure shear, outside its volume constraint.
        for (kind, f) in variants().into_iter().chain([(ElementKind::Quad4, Formulation::Hybrid { pressures: 1 })]) {
            let xyz = element(kind);
            let e = Twisted::new(Plane::new(kind, &xyz, f, Theory::Axisymmetric).unwrap()).unwrap();
            let u: Vec<[f64; 3]> = xyz.iter().map(|p| [0.0, 0.0, a * p[0] + b * p[1]]).collect();
            let points = e.stresses_with(&d, &u, None).unwrap();
            assert_eq!(points.len(), e.plane().rule().points().len());
            for p in points {
                let want = [-b * p.xyz[1] / p.xyz[0], b];
                assert!((p.strain[4] - want[0]).abs() < 1e-13 && (p.strain[5] - want[1]).abs() < 1e-13, "{kind:?} {f:?}: {:?}", p.strain);
                assert!((p.stress[4] - shear * want[0]).abs() < 1e-10 && (p.stress[5] - shear * want[1]).abs() < 1e-10, "{kind:?}: {:?}", p.stress);
                assert!(p.stress.iter().take(4).all(|s| s.abs() < 1e-10) && p.strain.iter().take(4).all(|s| s.abs() < 1e-13), "{kind:?}: {p:?}");
            }
        }
    }

    #[test]
    fn uniform_torsion_is_exact_where_the_element_holds_r_times_z() {
        // u_θ = κ r z: γ_zθ = κ r, γ_rθ = 0. Orthotropic: G13 acts on rθ, G23 on zθ.
        let d = orthotropic([1000.0, 1200.0, 900.0], [0.2, 0.25, 0.3], [300.0, 400.0, 500.0]).unwrap();
        let kappa = 0.01;
        let mid = |a: [f64; 3], b: [f64; 3]| [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, 0.0];
        let sq = [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 2.0, 0.0], [1.0, 2.0, 0.0]];
        let quad8 = [sq[0], sq[1], sq[2], sq[3], mid(sq[0], sq[1]), mid(sq[1], sq[2]), mid(sq[2], sq[3]), mid(sq[3], sq[0])];
        let tri = [[1.0, 0.2, 0.0], [3.0, 0.0, 0.0], [1.6, 1.9, 0.0]];
        let tri6 = [tri[0], tri[1], tri[2], mid(tri[0], tri[1]), mid(tri[1], tri[2]), mid(tri[2], tri[0])];
        let cases: [(ElementKind, &[[f64; 3]], Formulation); 5] = [
            (ElementKind::Quad4, &sq, Formulation::Full),
            (ElementKind::Quad4, &sq, Formulation::IncompatibleModes),
            (ElementKind::Quad8, &quad8, Formulation::Full),
            (ElementKind::Quad8, &quad8, Formulation::Reduced),
            (ElementKind::Tri6, &tri6, Formulation::Full),
        ];
        for (kind, xyz, f) in cases {
            let e = Twisted::new(Plane::new(kind, xyz, f, Theory::Axisymmetric).unwrap()).unwrap();
            let torsion: Vec<[f64; 3]> = xyz.iter().map(|p| [0.0, 0.0, kappa * p[0] * p[1]]).collect();
            for p in e.stresses_with(&d, &torsion, None).unwrap() {
                assert!(p.strain[4].abs() < 1e-13 && (p.strain[5] - kappa * p.xyz[0]).abs() < 1e-13, "{kind:?} {f:?}: {:?}", p.strain);
                assert!(p.stress[4].abs() < 1e-10 && (p.stress[5] - 500.0 * kappa * p.xyz[0]).abs() < 1e-10, "{kind:?} {f:?}: {:?}", p.stress);
            }
            // u_θ = b y strains rθ too: τ_rθ = G13 γ_rθ.
            let axial: Vec<[f64; 3]> = xyz.iter().map(|p| [0.0, 0.0, 0.004 * p[1]]).collect();
            for p in e.stresses_with(&d, &axial, None).unwrap() {
                assert!((p.strain[4] + 0.004 * p.xyz[1] / p.xyz[0]).abs() < 1e-13, "{kind:?}: {:?}", p.strain);
                assert!((p.stress[4] - 400.0 * p.strain[4]).abs() < 1e-10 && (p.stress[5] - 500.0 * 0.004).abs() < 1e-10, "{kind:?}: {:?}", p.stress);
            }
        }
    }

    #[test]
    fn what_the_twist_refuses() {
        let d = isotropic(1000.0, 0.3).unwrap();
        let sq = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        let slab = Plane::new(ElementKind::Quad4, &sq, Formulation::Full, Theory::PlaneStrain { thickness: 1.0 }).unwrap();
        assert!(matches!(Twisted::new(slab), Err(Error::Invalid(_))));
        let e = Twisted::new(Plane::new(ElementKind::Quad4, &sq, Formulation::Full, Theory::Axisymmetric).unwrap()).unwrap();
        // A material whose twist shears act on the meridian strains.
        let mut coupled = d;
        coupled[0][4] = 1.0;
        coupled[4][0] = 1.0;
        assert!(matches!(e.stiffness(&coupled), Err(Error::NotSupported(_))));
        assert!(matches!(e.stresses_with(&coupled, &[[0.0; 3]; 4], None), Err(Error::NotSupported(_))));
        assert!(e.stresses_with(&d, &[[0.0; 3]; 3], None).is_err());
        // An element touching the axis is fine (its integration points are off it).
        let on_axis = sq.map(|p| [p[0] - 1.0, p[1], 0.0]);
        let e = Twisted::new(Plane::new(ElementKind::Quad4, &on_axis, Formulation::IncompatibleModes, Theory::Axisymmetric).unwrap()).unwrap();
        assert!(e.stiffness(&d).is_ok());
    }
}
