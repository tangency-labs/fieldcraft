//! Heat conduction kernels: conductivity, heat capacity, internal heat sources, surface fluxes and
//! convection, and the heat flux at the integration points.
//!
//! The Galerkin form of transient conduction ρc Ṫ = ∇·(k∇T) + Q with flux and convection on the
//! boundary (Zienkiewicz, Taylor & Zhu 2013, §3.2 and ch. 7 of vol. 1; Hughes 2000, ch. 1, §8.1)
//! gives, element by element, with one temperature per node ([`ElementDofs`] `TEMPERATURE`):
//! - conductivity K = ∫ ∇Nᵀ k ∇N dV, k the diagonal conductivity in the global axes (the
//!   material axes until orientations land), integrated with the element's full rule;
//! - capacity C = ∫ ρc N Nᵀ dV with the mass rule (consistent);
//! - source f = ∫ N Q dV;
//! - on faces: a flux q into the body gives f = ∫ N q dA; convection with film coefficient h to a
//!   sink at T∞ adds h ∫ N Nᵀ dA to K and h T∞ ∫ N dA to f.
//!
//! Faces are integrated on their own (possibly curved) geometry with the model's face rule
//! (the one pressures use), exact for N Nᵀ on flat faces.

use fieldcraft_model::ElementKind;

use crate::basis::Basis;
use crate::dense::Matrix;
use crate::dof::{DofSet, ElementDofs};
use crate::element::{ElementType, element_type};
use crate::quadrature::Rule;
use crate::{Error, Result};

/// What a heat conduction element carries: one temperature per node.
pub const HEAT_DOFS: ElementDofs = ElementDofs::new(DofSet::TEMPERATURE);

/// The heat flux −k∇T at one integration point.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PointFlux {
    pub xyz: [f64; 3],
    pub flux: [f64; 3],
}

/// A solid element conducting heat.
#[derive(Debug, Clone, Copy)]
pub struct Conduction<'a> {
    pub ty: &'static ElementType,
    basis: Basis,
    xyz: &'a [[f64; 3]],
}

impl<'a> Conduction<'a> {
    pub fn new(kind: ElementKind, xyz: &'a [[f64; 3]]) -> Result<Conduction<'a>> {
        let ty = element_type(kind)?;
        if xyz.len() != ty.nodes {
            return Err(Error::Invalid(format!("a {} element has {} nodes, got {}", kind.name(), ty.nodes, xyz.len())));
        }
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        Ok(Conduction { ty, basis: Basis::of(kind), xyz })
    }

    /// Visit the points of `rule`: shape values, physical gradients and the weight w det J.
    fn each(&self, rule: Rule, mut visit: impl FnMut(&[f64], &[[f64; 3]], f64)) -> Result<()> {
        for p in rule.points() {
            let s = self.basis.eval(p.xi);
            let g = self.basis.gradients(&s, self.xyz)?;
            visit(s.values(), g.gradients(), g.det * p.weight);
        }
        Ok(())
    }

    /// K = ∫ ∇Nᵀ k ∇N dV, n × n.
    pub fn conductivity(&self, k: [f64; 3]) -> Result<Matrix> {
        let mut m = Matrix::zeros(self.ty.nodes);
        self.each(self.ty.full, |_, g, w| {
            for (i, gi) in g.iter().enumerate() {
                for (j, gj) in g.iter().enumerate() {
                    m.add(i, j, w * (k[0] * gi[0] * gj[0] + k[1] * gi[1] * gj[1] + k[2] * gi[2] * gj[2]));
                }
            }
        })?;
        Ok(m)
    }

    /// C = ∫ ρc N Nᵀ dV, n × n; `capacity` is ρc (energy per volume per temperature).
    pub fn capacity(&self, capacity: f64) -> Result<Matrix> {
        let mut m = Matrix::zeros(self.ty.nodes);
        self.each(self.ty.mass, |n, _, w| {
            for (i, a) in n.iter().enumerate() {
                for (j, b) in n.iter().enumerate() {
                    m.add(i, j, w * capacity * a * b);
                }
            }
        })?;
        Ok(m)
    }

    /// f = ∫ N Q dV for a uniform heat source Q per volume.
    pub fn source(&self, q: f64) -> Result<Vec<f64>> {
        let mut f = vec![0.0; self.ty.nodes];
        self.each(self.ty.mass, |n, _, w| {
            for (fi, ni) in f.iter_mut().zip(n) {
                *fi += w * q * ni;
            }
        })?;
        Ok(f)
    }

    /// The heat flux −k∇T at the full rule's points for nodal temperatures `t`.
    pub fn flux(&self, k: [f64; 3], t: &[f64]) -> Result<Vec<PointFlux>> {
        if t.len() != self.ty.nodes {
            return Err(Error::Invalid(format!("{} temperatures for a {}-node element", t.len(), self.ty.nodes)));
        }
        let mut out = Vec::new();
        self.each(self.ty.full, |n, g, _| {
            let (mut grad, mut xyz) = ([0.0; 3], [0.0; 3]);
            for ((gi, ti), (ni, xi)) in g.iter().zip(t).zip(n.iter().zip(self.xyz)) {
                for c in 0..3 {
                    grad[c] += gi[c] * ti;
                    xyz[c] += ni * xi[c];
                }
            }
            out.push(PointFlux { xyz, flux: [0, 1, 2].map(|c| -k[c] * grad[c]) });
        })?;
        Ok(out)
    }
}

/// Visit the points of a face (node coordinates in the face's node order): shape values and
/// the area weight, with the model's face rule ([`fieldcraft_model::nodal::integrate_face`], the
/// one pressures and tractions use).
fn each_face_point(xyz: &[[f64; 3]], mut visit: impl FnMut(&[f64], f64)) -> Result<()> {
    let mut collapsed = false;
    fieldcraft_model::nodal::integrate_face(xyz, |n, nda| {
        let area = (nda[0] * nda[0] + nda[1] * nda[1] + nda[2] * nda[2]).sqrt();
        collapsed |= !(area > 0.0 && area.is_finite());
        visit(n, area);
    })
    .map_err(|e| Error::Invalid(e.to_string()))?;
    if collapsed {
        return Err(Error::Degenerate("a face is collapsed (zero area)".into()));
    }
    Ok(())
}

/// ∫ N dA for every node of a face: the nodal shares of a uniform flux.
pub fn face_weights(xyz: &[[f64; 3]]) -> Result<Vec<f64>> {
    let mut f = vec![0.0; xyz.len()];
    each_face_point(xyz, |n, w| {
        for (fi, ni) in f.iter_mut().zip(n) {
            *fi += w * ni;
        }
    })?;
    Ok(f)
}

/// h ∫ N Nᵀ dA for a face: the conductivity convection adds.
pub fn face_film(xyz: &[[f64; 3]], film: f64) -> Result<Matrix> {
    let mut m = Matrix::zeros(xyz.len());
    each_face_point(xyz, |n, w| {
        for (i, a) in n.iter().enumerate() {
            for (j, b) in n.iter().enumerate() {
                m.add(i, j, w * film * a * b);
            }
        }
    })?;
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shape::{HEX_CORNERS, TET_CORNERS};

    fn box_hex(a: [f64; 3]) -> Vec<[f64; 3]> {
        HEX_CORNERS.iter().map(|c| [0, 1, 2].map(|k| (c[k] + 1.0) / 2.0 * a[k])).collect()
    }

    #[test]
    fn conductivity_rows_sum_to_zero_and_reproduce_a_linear_field() {
        // A uniform gradient G in a 2 × 1 × 3 box: K T = boundary flux, zero at interior
        // balance; here check K·1 = 0 and Tᵀ K T = k |G|² V.
        let xyz = box_hex([2.0, 1.0, 3.0]);
        for kind in [ElementKind::Hex8] {
            let e = Conduction::new(kind, &xyz).unwrap();
            let k = [2.0, 3.0, 5.0];
            let m = e.conductivity(k).unwrap();
            let ones = vec![1.0; 8];
            assert!(m.mul_vec(&ones).iter().all(|v| v.abs() < 1e-12));
            let g = [0.3, -0.7, 0.2];
            let t: Vec<f64> = xyz.iter().map(|x| g[0] * x[0] + g[1] * x[1] + g[2] * x[2]).collect();
            let kt = m.mul_vec(&t);
            let energy: f64 = kt.iter().zip(&t).map(|(a, b)| a * b).sum();
            let exact = (k[0] * g[0] * g[0] + k[1] * g[1] * g[1] + k[2] * g[2] * g[2]) * 6.0;
            assert!((energy - exact).abs() < 1e-12 * exact, "{energy} vs {exact}");
            let fl = e.flux(k, &t).unwrap();
            for p in fl {
                for c in 0..3 {
                    assert!((p.flux[c] + k[c] * g[c]).abs() < 1e-13);
                }
            }
        }
    }

    #[test]
    fn capacity_and_source_integrate_the_volume() {
        let tet: Vec<[f64; 3]> = TET_CORNERS.iter().map(|c| [2.0 * c[0], 3.0 * c[1], c[2]]).collect();
        let e = Conduction::new(ElementKind::Tet4, &tet).unwrap();
        let v = 1.0;
        let c = e.capacity(4.0).unwrap();
        assert!((c.data.iter().sum::<f64>() - 4.0 * v).abs() < 1e-13);
        assert!((c.get(0, 0) - 4.0 * v / 10.0).abs() < 1e-13 && (c.get(0, 1) - 4.0 * v / 20.0).abs() < 1e-13);
        let f = e.source(3.0).unwrap();
        assert!(f.iter().all(|x| (x - 0.75).abs() < 1e-13), "{f:?}");
        assert!(Conduction::new(ElementKind::Tet4, &tet[..3]).is_err());
        assert!(e.flux([1.0; 3], &[0.0; 3]).is_err());
    }

    #[test]
    fn face_weights_and_film() {
        let quad = [[0., 0., 0.], [2., 0., 0.], [2., 3., 0.], [0., 3., 0.]];
        let w = face_weights(&quad).unwrap();
        assert!(w.iter().all(|x| (x - 1.5).abs() < 1e-14));
        let m = face_film(&quad, 2.0).unwrap();
        // h A / 9 on the diagonal of a bilinear rectangle, sum h A.
        assert!((m.data.iter().sum::<f64>() - 12.0).abs() < 1e-12);
        assert!((m.get(0, 0) - 2.0 * 6.0 / 9.0).abs() < 1e-13);
        let mut t6 = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            t6.push([0, 1, 2].map(|k| (t6[a][k] + t6[b][k]) / 2.0));
        }
        let w = face_weights(&t6).unwrap();
        for (i, x) in w.iter().enumerate() {
            let want = if i < 3 { 0.0 } else { 0.5 / 3.0 };
            assert!((x - want).abs() < 1e-14, "{i}: {x}");
        }
        assert!(face_weights(&quad[..2]).is_err());
        assert!(face_weights(&[[0.0; 3]; 4]).is_err());
        assert!(face_film(&[[f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], 1.0).is_err());
    }
}
