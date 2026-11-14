//! Backend-neutral linear products in an operator's declared degree-of-freedom space.

use crate::sparse::Csr;
use crate::{Error, Result};

/// A square linear map. `apply` checks both vector lengths against `dimension` and overwrites
/// the complete output on success; output after an error is unspecified. Applying an operator
/// must not commit material history or change the physical state represented by the operator.
pub trait LinearOperator {
    fn dimension(&self) -> usize;
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()>;
}

impl LinearOperator for Csr {
    fn dimension(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
        if x.len() != self.n || y.len() != self.n {
            return Err(Error::Invalid("linear operator vector dimensions do not match".into()));
        }
        self.mul(x, y);
        Ok(())
    }
}
