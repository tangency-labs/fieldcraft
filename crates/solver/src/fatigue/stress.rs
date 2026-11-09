//! Scalar nominal stress from complete global symmetric tensors, in Pa.
//! Order xx,yy,zz,xy,xz,yz; shear entries are tensor stresses, not engineering strains.
//! Signed VM follows Altair's public Uniaxial Fatigue Analysis, Equivalent Nominal
//! Stress: sign of the absolute-largest principal stress. Opposite-sign ties are
//! refused instead of inventing a reversal. Principal eigenvalues use Smith (1961)
//! through fieldcraft-fem. These uniaxial S-N reductions are restricted to verified
//! proportional tensors, not a general nonproportional/critical-plane fatigue law.
use crate::{Error, Result};

const PROPORTIONAL_TOL: f64 = 1e-9; // unit-factors: dimensionless normalized proportionality residual
const SIGN_TIE_TOL: f64 = 1e-10;

#[derive(Debug, Clone, Copy, serde::Deserialize, serde::Serialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Criterion {
    Projection {
        weights: [f64; 6],
    },
    SignedVonMises,
    /// Algebraically largest eigenvalue, NOT the largest absolute principal value.
    MaxPrincipal,
}
impl Criterion {
    pub fn validate(self) -> Result<()> {
        if let Self::Projection { weights } = self
            && (weights.iter().any(|v| !v.is_finite()) || weights.iter().all(|v| *v == 0.))
        {
            return Err(Error::Invalid("need finite nonzero stress projection weights".into()));
        }
        Ok(())
    }
    pub fn scalarizer(self) -> Result<Scalarizer> {
        self.validate()?;
        Ok(Scalarizer { criterion: self, reference: None })
    }
}
pub struct Scalarizer {
    criterion: Criterion,
    reference: Option<[f64; 6]>,
}
impl Scalarizer {
    pub fn push(&mut self, stress: &[f64; 6]) -> Result<f64> {
        if stress.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("need six finite stress components".into()));
        }
        let value = match self.criterion {
            Criterion::Projection { weights } => stress.iter().zip(weights).map(|(v, w)| v * w).sum(),
            _ => {
                let scale = stress.iter().map(|v| v.abs()).fold(0.0, f64::max);
                if scale == 0. {
                    return Ok(0.);
                }
                let normalized = stress.map(|v| v / scale);
                if let Some(reference) = self.reference {
                    let pivot =
                        reference.iter().position(|v| v.abs() == 1.).ok_or_else(|| Error::Invalid("invalid proportional reference".into()))?;
                    let r = reference.get(pivot).copied().ok_or_else(|| Error::Invalid("missing reference component".into()))?;
                    let factor = normalized.get(pivot).copied().ok_or_else(|| Error::Invalid("missing tensor component".into()))? / r;
                    if normalized.iter().zip(reference).any(|(v, r)| (*v - factor * r).abs() > PROPORTIONAL_TOL) {
                        return Err(Error::Invalid("nonproportional tensor history: signed VM/max principal nominal S-N reduction is not supported; use a justified fixed projection or a multiaxial method".into()));
                    }
                } else {
                    self.reference = Some(normalized);
                }
                let [largest, _, smallest] = fieldcraft_fem::principal(&normalized);
                match self.criterion {
                    Criterion::MaxPrincipal => largest * scale,
                    Criterion::SignedVonMises => {
                        let vm = fieldcraft_fem::von_mises(&normalized);
                        if vm == 0. {
                            return Ok(0.);
                        }
                        let size = largest.abs().max(smallest.abs());
                        if largest > 0. && smallest < 0. && (largest.abs() - smallest.abs()).abs() <= SIGN_TIE_TOL * size {
                            return Err(Error::Invalid("signed von Mises has ambiguous opposite-sign absolute-principal tie; supply a justified fixed projection, not a PSD-derived sign".into()));
                        }
                        let signed = if smallest.abs() > largest.abs() { smallest } else { largest };
                        vm.copysign(signed) * scale
                    }
                    Criterion::Projection { .. } => return Err(Error::Invalid("invalid tensor criterion".into())),
                }
            }
        };
        if !value.is_finite() {
            return Err(Error::Invalid("scalar fatigue stress overflow".into()));
        }
        Ok(value)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fatigue::{Basquin, life::StressMean, miner, rainflow};
    fn signal(c: Criterion, tensor: [f64; 6], factors: &[f64]) -> Vec<f64> {
        let mut scalar = c.scalarizer().unwrap();
        factors.iter().map(|f| scalar.push(&tensor.map(|v| v * f)).unwrap()).collect()
    }
    #[test]
    fn proportional_tension_compression_and_biaxial_shear_to_miner() {
        let curve = Basquin { amplitude_ref: 1., cycles_ref: 100., exponent: 2., endurance_limit: 0. };
        // Plane-stress characteristic polynomial: p(p-4)(p+1), VM²=3²+3*2²=21.
        for tensor in [[3., 0., 0., 2., 0., 0.], [-0.5, 3.5, 0., 1.5, 0., 0.]] {
            let s = signal(Criterion::SignedVonMises, tensor, &[-1., 1., -1.]);
            let cycles = rainflow(&s).unwrap();
            assert!((miner(&cycles, &curve, StressMean::None).unwrap() - 0.21).abs() < 1e-13);
            let s = signal(Criterion::MaxPrincipal, tensor, &[-1., 1., -1.]);
            for (a, b) in s.iter().zip([1., 4., 1.]) {
                assert!((a - b).abs() < 1e-13);
            }
            let cycles = rainflow(&s).unwrap();
            assert!(cycles.iter().all(|c| (c.mean - 2.5).abs() < 1e-13));
            assert!((miner(&cycles, &curve, StressMean::Goodman { ultimate: 10. }).unwrap() - 0.04).abs() < 1e-13);
        }
        assert_eq!(signal(Criterion::SignedVonMises, [2., 0., 0., 0., 0., 0.], &[-1., 1., -1.]), vec![-2., 2., -2.]);
        assert_eq!(signal(Criterion::MaxPrincipal, [2., 0., 0., 0., 0., 0.], &[-1., 1., -1.]), vec![0., 2., 0.]);
        assert_eq!(signal(Criterion::MaxPrincipal, [3., 2., 1., 0., 0., 0.], &[-1., 1., -1.]), vec![-1., 3., -1.]);
    }
    #[test]
    fn ambiguous_shear_nonproportional_and_nonfinite_are_refused() {
        let mut s = Criterion::SignedVonMises.scalarizer().unwrap();
        assert!(s.push(&[0., 0., 0., 1., 0., 0.]).is_err());
        for c in [Criterion::SignedVonMises, Criterion::MaxPrincipal] {
            let mut s = c.scalarizer().unwrap();
            s.push(&[1., 0., 0., 0., 0., 0.]).unwrap();
            assert!(s.push(&[0., 1., 0., 0., 0., 0.]).is_err());
            assert!(s.push(&[f64::NAN; 6]).is_err());
        }
    }
}
