//! Hybrid (mixed displacement–pressure) solids for nearly incompressible materials.
//!
//! The isotropic elasticity splits into a deviatoric part and the bulk modulus,
//! D = D_dev + K m mᵀ with m = (1, 1, 1, 0, 0, 0) and K = mᵀ D m / 9. The pressure (here the
//! mean stress, positive in tension) is an independent field p = N_p · p̂, discontinuous between
//! elements, and the volumetric strain is enforced weakly (the perturbed Lagrangian of
//! T. J. R. Hughes, *The Finite Element Method*, Dover 2000, §4.3–4.4):
//!
//! K_dev = ∫ Bᵀ D_dev B dV, G = ∫ Bᵀ m N_pᵀ dV, M = ∫ N_p N_pᵀ dV,
//! [K_dev G; Gᵀ −M/K] [u; p̂] = [f; 0],
//!
//! so p̂ = K M⁻¹ Gᵀ u is condensed out element by element and the element stiffness is
//! K_dev + K G M⁻¹ Gᵀ: symmetric positive definite for any ν < ½, and bounded volumetric
//! stiffness per element however close ν is to ½ — the constraint count no longer grows with
//! the nodes, which is what locks the displacement elements. With a constant pressure this is
//! the B-bar / mean dilatation element (Hughes 2000, §4.5; J. C. Nagtegaal, D. M. Parks,
//! J. R. Rice, *On numerically accurate finite element solutions in the fully plastic range*,
//! CMAME 4(2):153–177, 1974).
//!
//! Pressure fields and stability (the inf-sup or LBB condition, F. Brezzi, M. Fortin, *Mixed
//! and Hybrid Finite Element Methods*, Springer 1991, ch. VI; K.-J. Bathe, *Finite Element
//! Procedures*, 2nd ed., 2014, §4.5):
//! - **hex20**: linear pressure 1, x, y, z (Q2–P1disc, inf-sup stable);
//! - **tet10**: constant pressure (P2–P0, inf-sup stable; P2–P1disc is not on general meshes);
//! - **tet10**: a single element has three more zero-energy modes (dilatations with a linear
//!   volumetric strain, which its deviatoric stiffness does not see and one pressure controls
//!   only in the mean); neighbouring elements restrain them, as they do the reduced hex20's;
//! - **hex8**: constant pressure (Q1–P0). It passes the patch test and does not lock, but is
//!   not inf-sup stable: on some meshes and constraints the element pressures show a
//!   checkerboard mode while displacements stay good; nodal stresses are averaged anyway;
//! - **tet4, wedge6**: refused — their cross-section is the linear triangle, and one pressure
//!   per element still leaves too few displacement unknowns per constraint: they lock (verified
//!   on the thick cylinder; a nodal-pressure formulation would be needed).
//!
//! A continuous nodal pressure (Taylor–Hood, P on the corner nodes) makes the global system a
//! symmetric indefinite saddle point that the Cholesky factorization cannot solve; it waits for
//! an indefinite (LDLᵀ) factorization.

use fieldcraft_model::ElementKind;

use crate::dense::{Cholesky, Matrix};
use crate::element::{Formulation, Solid, ThermalStrain, accumulate_symmetric, b_block};
use crate::material::Elasticity;
use crate::{Error, Result};

/// The pressure unknowns of the hybrid formulation of an element kind, if it has one.
pub fn pressures(kind: ElementKind) -> Option<usize> {
    match kind {
        ElementKind::Hex20 => Some(4),
        ElementKind::Hex8 | ElementKind::Tet10 => Some(1),
        _ => None,
    }
}

/// Two moduli are equal when they differ by less than this fraction of the largest.
const SAME_MODULUS: f64 = 1e-12; // unit-factors: relative to the largest modulus

/// The deviatoric elasticity and the bulk modulus of an isotropic elasticity; orthotropic
/// materials have no such split and are refused.
pub fn split(d: &Elasticity) -> Result<(Elasticity, f64)> {
    let scale = d.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    let same = |a: f64, b: f64| (a - b).abs() <= SAME_MODULUS * scale;
    let (c11, c12, c44) = (d[0][0], d[0][1], d[3][3]);
    let isotropic = (0..3).all(|i| same(d[i][i], c11) && same(d[i + 3][i + 3], c44))
        && (0..3).all(|i| (0..3).all(|j| i == j || same(d[i][j], c12)))
        && same(c44, 0.5 * (c11 - c12))
        && (0..6).all(|i| (0..6).all(|j| (i < 3 && j < 3) || i == j || d[i][j] == 0.0));
    if !isotropic {
        return Err(Error::NotSupported("hybrid elements of orthotropic materials".into()));
    }
    let bulk = (c11 + 2.0 * c12) / 3.0;
    let mut dev = *d;
    for row in dev.iter_mut().take(3) {
        for v in row.iter_mut().take(3) {
            *v -= bulk;
        }
    }
    Ok((dev, bulk))
}

/// The parts of a hybrid element: G, 3n × np row-major; the factor of M; and K_dev when
/// asked for.
struct Mixed {
    kdev: Option<Matrix>,
    g: Vec<f64>,
    m: Cholesky,
    /// With a thermal strain: ∫ N_p tr ε_th dV (np) and ∫ Bᵀ D_dev ε_th dV (3n).
    h: Vec<f64>,
    f_dev: Vec<f64>,
}

/// The pressure field of a hybrid element under nodal displacements, for its stresses.
pub(crate) struct PressureField {
    dev: Elasticity,
    bulk: f64,
    /// M⁻¹ Gᵀ u, the coefficients of the volumetric strain's projection.
    coefficients: [f64; 4],
    np: usize,
    frame: ([f64; 3], f64),
}

impl PressureField {
    /// D_dev ε + p m at `x`.
    pub(crate) fn stress(&self, strain: &[f64; 6], x: [f64; 3]) -> [f64; 6] {
        let basis = pressure_basis(self.np, x, &self.frame);
        let p = self.bulk * self.coefficients.iter().zip(basis).map(|(a, b)| a * b).sum::<f64>();
        let mut s = crate::material::stress(&self.dev, strain);
        for v in s.iter_mut().take(3) {
            *v += p;
        }
        s
    }
}

/// The pressure functions at a point: 1, then (x − c)/h, (y − c)/h, (z − c)/h for the linear
/// field (`np` = 4).
pub(crate) fn pressure_basis(np: usize, x: [f64; 3], frame: &([f64; 3], f64)) -> [f64; 4] {
    let (c, h) = frame;
    let mut out = [1.0, 0.0, 0.0, 0.0];
    if np == 4 && *h > 0.0 {
        for k in 0..3 {
            out[k + 1] = (x[k] - c[k]) / h;
        }
    }
    out
}

impl Solid<'_> {
    /// The pressure unknowns of this (hybrid) element.
    pub(crate) fn np(&self) -> usize {
        match self.formulation {
            Formulation::Hybrid { pressures } => pressures.min(4),
            _ => 0,
        }
    }

    /// The element's centre (mean of its corners) and size (largest corner distance from it),
    /// which scale the linear pressure functions.
    pub(crate) fn frame(&self) -> ([f64; 3], f64) {
        let corners = self.xyz.get(..self.ty.corners).unwrap_or_default();
        let n = corners.len().max(1) as f64;
        let c = corners.iter().fold([0.0; 3], |a, x| [a[0] + x[0] / n, a[1] + x[1] / n, a[2] + x[2] / n]);
        let size = corners.iter().map(|x| ((x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2) + (x[2] - c[2]).powi(2)).sqrt()).fold(0.0, f64::max);
        (c, size)
    }

    /// G and M (and K_dev with `kdev`, and the thermal terms with `thermal`) in one pass over
    /// the stiffness rule.
    fn mixed(&self, dev: &Elasticity, kdev: bool, thermal: Option<&ThermalStrain>) -> Result<Mixed> {
        let (np, nodes, frame) = (self.np(), self.ty.nodes, self.frame());
        let mut kdev = kdev.then(|| Matrix::zeros(3 * nodes));
        let mut g = vec![0.0; 3 * nodes * np];
        let mut m = Matrix::zeros(np);
        let mut h = vec![0.0; np];
        let mut f_dev = vec![0.0; 3 * nodes];
        for p in self.rule().points() {
            let gr = self.gradients(p.xi)?;
            let w = gr.det * p.weight;
            let gu = gr.dx.get(..nodes).unwrap_or_default();
            if let Some(k) = kdev.as_mut() {
                accumulate_symmetric(&mut k.data, gu, dev, w);
            }
            let x = if np == 4 { self.position_of(gr.shape.values()) } else { [0.0; 3] };
            let np_at = pressure_basis(np, x, &frame);
            if let Some(t) = thermal {
                let eth = Self::thermal_strain(t, gr.shape.values());
                let volumetric = eth[0] + eth[1] + eth[2];
                for (hk, nk) in h.iter_mut().zip(np_at) {
                    *hk += w * nk * volumetric;
                }
                let sig = crate::material::stress(dev, &eth);
                for (i, gi) in gu.iter().enumerate() {
                    let b = b_block(gi);
                    for c in 0..3 {
                        if let Some(f) = f_dev.get_mut(3 * i + c) {
                            *f += w * (0..6).map(|s| b[s][c] * sig[s]).sum::<f64>();
                        }
                    }
                }
            }
            // Bᵢᵀ m is node i's gradient.
            for (i, gi) in gu.iter().enumerate() {
                for c in 0..3 {
                    for k in 0..np {
                        if let Some(v) = g.get_mut((3 * i + c) * np + k) {
                            *v += w * gi[c] * np_at[k];
                        }
                    }
                }
            }
            for a in 0..np {
                for b in 0..np {
                    m.add(a, b, w * np_at[a] * np_at[b]);
                }
            }
        }
        Ok(Mixed { kdev, g, m: Cholesky::new(&m)?, h, f_dev })
    }

    /// K_dev + K G M⁻¹ Gᵀ, symmetric by construction (H = M⁻¹ Gᵀ once, upper triangle mirrored).
    pub(crate) fn hybrid_stiffness(&self, d: &Elasticity) -> Result<Matrix> {
        let (dev, bulk) = split(d)?;
        let Mixed { kdev, g, m, .. } = self.mixed(&dev, true, None)?;
        let mut k = kdev.ok_or_else(|| Error::Invalid("hybrid stiffness".into()))?;
        let (n, np) = (k.n, self.np());
        let row = |i: usize| g.get(i * np..(i + 1) * np).unwrap_or_default();
        let h: Vec<Vec<f64>> = (0..n)
            .map(|j| {
                let mut y = row(j).to_vec();
                m.solve(&mut y);
                y
            })
            .collect();
        for i in 0..n {
            for (j, hj) in h.iter().enumerate().skip(i) {
                let v = bulk * row(i).iter().zip(hj).map(|(a, b)| a * b).sum::<f64>();
                k.add(i, j, v);
                if i != j {
                    k.add(j, i, v);
                }
            }
        }
        Ok(k)
    }

    /// The equivalent nodal forces of a thermal strain: ∫ Bᵀ D_dev ε_th dV for the deviatoric
    /// part, and K G M⁻¹ ∫ N_p tr ε_th dV for the volumetric one through the pressure.
    pub(crate) fn hybrid_thermal_load(&self, d: &Elasticity, thermal: &ThermalStrain) -> Result<Vec<[f64; 3]>> {
        let (dev, bulk) = split(d)?;
        let np = self.np();
        let Mixed { g, m, mut h, f_dev, .. } = self.mixed(&dev, false, Some(thermal))?;
        m.solve(&mut h);
        let f: Vec<f64> = f_dev
            .iter()
            .enumerate()
            .map(|(i, fd)| fd + bulk * (0..np).map(|k| g.get(i * np + k).copied().unwrap_or(0.0) * h.get(k).copied().unwrap_or(0.0)).sum::<f64>())
            .collect();
        Ok(f.as_chunks::<3>().0.to_vec())
    }

    /// The element's pressure field under the nodal displacements `u` and a thermal strain:
    /// p̂ = K M⁻¹ (Gᵀ u − ∫ N_p tr ε_th dV).
    pub(crate) fn pressure_field(&self, d: &Elasticity, u: &[[f64; 3]], thermal: Option<&ThermalStrain>) -> Result<PressureField> {
        let (dev, bulk) = split(d)?;
        let np = self.np();
        let Mixed { g, m, h, .. } = self.mixed(&dev, false, thermal)?;
        let mut pc: Vec<f64> = h.iter().map(|x| -x).collect();
        for (i, ui) in u.iter().flatten().enumerate() {
            for (k, x) in pc.iter_mut().enumerate() {
                *x += g.get(i * np + k).copied().unwrap_or(0.0) * ui;
            }
        }
        m.solve(&mut pc);
        let mut coefficients = [0.0; 4];
        for (c, v) in coefficients.iter_mut().zip(pc) {
            *c = v;
        }
        Ok(PressureField { dev, bulk, coefficients, np, frame: self.frame() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{isotropic, orthotropic};

    #[test]
    fn the_split_is_deviatoric_plus_bulk() {
        let d = isotropic(1000.0, 0.3).unwrap();
        let (dev, k) = split(&d).unwrap();
        assert!((k - 1000.0 / (3.0 * (1.0 - 0.6))).abs() < 1e-9);
        // D_dev m = 0.
        for r in 0..6 {
            assert!((0..3).map(|c| dev[r][c]).sum::<f64>().abs() < 1e-9);
        }
        assert!(split(&orthotropic([1000.0, 800.0, 600.0], [0.3, 0.25, 0.2], [400.0, 300.0, 350.0]).unwrap()).is_err());
        assert_eq!(pressures(ElementKind::Tet4), None);
        assert_eq!(pressures(ElementKind::Hex20), Some(4));
    }
}
