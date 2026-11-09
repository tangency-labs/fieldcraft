//! One immutable, exact-input memo of an already validated modal basis.
//! Step Rayleigh/loss damping and sweep points are evaluated later by `Modal::solve`;
//! they do not change the original-pencil eigenpairs or load-dependent Ritz vectors.
use super::{Matrices, Modal, ModalData};
use crate::{Limits, Result, complex::Vector, sparse::Csr};
use std::sync::Arc;

#[derive(Debug, Clone, Default)]
pub struct ModalCache {
    entry: Option<Arc<Entry>>,
}

#[derive(Debug)]
struct Entry {
    matrices: [Option<Csr>; 4],
    force: Vector,
    count: usize,
    cycles: usize,
    limits: Limits,
    data: Arc<ModalData>,
}

fn same_values(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}
fn same_matrix(a: Option<&Csr>, b: Option<&Csr>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => a.n == b.n && a.ptr == b.ptr && a.idx == b.idx && same_values(&a.val, &b.val),
        _ => false,
    }
}
impl ModalCache {
    /// Drop retained inputs and basis. No cache is serialized with solver results.
    pub fn clear(&mut self) {
        self.entry = None;
    }

    pub(super) fn get<'a>(&mut self, matrices: &Matrices<'a>, count: usize, force: &Vector, cycles: usize, limits: &Limits) -> Result<Modal<'a>> {
        let inputs = [Some(matrices.stiffness), Some(matrices.mass), matrices.viscous, matrices.loss];
        let lookup = crate::timing::Span::new("harmonic_cache_lookup");
        if let Some(e) = &self.entry
            && e.count == count
            && e.cycles == cycles
            && e.limits == *limits
            && same_values(&e.force.real, &force.real)
            && same_values(&e.force.imaginary, &force.imaginary)
            && inputs.iter().zip(&e.matrices).all(|(a, b)| same_matrix(*a, b.as_ref()))
        {
            drop(lookup);
            let _hit = crate::timing::Span::new("harmonic_cache_hit");
            return Ok(Modal::from_data(matrices, Arc::clone(&e.data)));
        }
        drop(lookup);
        self.clear();
        let modal = Modal::with_residual_vectors(matrices, count, force, cycles, limits)?;
        // Account for retained values and indices in conservative eight-byte slots.
        let storage = inputs
            .iter()
            .flatten()
            .try_fold(0usize, |n, m| n.checked_add(m.ptr.len())?.checked_add(m.idx.len())?.checked_add(m.val.len()))
            .and_then(|n| n.checked_add(force.real.len())?.checked_add(force.imaginary.len()))
            .and_then(|n| modal.data.shapes.iter().try_fold(n, |n, v| n.checked_add(v.len())))
            .and_then(|n| n.checked_add(modal.data.stiffness.len().checked_mul(4)?));
        if storage.is_some_and(|n| n <= limits.max_factor_entries) {
            let _store = crate::timing::Span::new("harmonic_cache_store");
            self.entry = Some(Arc::new(Entry {
                matrices: inputs.map(|m| m.cloned()),
                force: force.clone(),
                count,
                cycles,
                limits: *limits,
                data: Arc::clone(&modal.data),
            }));
        }
        Ok(modal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn matrices() -> (Csr, Csr) {
        let matrix = |values| Csr { n: 8, ptr: (0..=8).collect(), idx: (0..8).collect(), val: values };
        (matrix((1..=8).map(|x| (x * x) as f64).collect()), matrix(vec![1.0; 8]))
    }
    #[test]
    fn exact_basis_reuse_keeps_response_and_rejects_changed_inputs() -> Result<()> {
        let (k, m) = matrices();
        let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
        let force = Vector { real: vec![1.0; 8], imaginary: vec![0.0; 8] };
        let limits = Limits::default();
        let mut cache = ModalCache::default();
        let cold = cache.get(&matrices, 3, &force, 2, &limits)?;
        let warm = cache.get(&matrices, 3, &force, 2, &limits)?;
        assert!(Arc::ptr_eq(&cold.data, &warm.data));
        for omega in [0.0, 0.5, 3.5] {
            let d = fieldcraft_model::harmonic::Damping { alpha: 0.1, beta: 0.02, loss_factor: 0.01 };
            let a = cold.solve(omega, d, &[], &force, &limits)?;
            let b = warm.solve(omega, d, &[], &force, &limits)?;
            assert!(same_values(&a.displacement.real, &b.displacement.real));
            assert!(same_values(&a.displacement.imaginary, &b.displacement.imaginary));
        }
        let original = cache.clone();
        let mut different_force = force.clone();
        different_force.real[0] = 2.0;
        assert!(!Arc::ptr_eq(&cold.data, &cache.get(&matrices, 3, &different_force, 2, &limits)?.data));
        cache = original.clone();
        let mut changed = k.clone();
        changed.val[0] = 1.5;
        assert!(!Arc::ptr_eq(&cold.data, &cache.get(&Matrices { stiffness: &changed, ..matrices }, 3, &force, 2, &limits)?.data));
        cache = original.clone();
        let mut changed = m.clone();
        changed.val[0] = 1.5;
        assert!(!Arc::ptr_eq(&cold.data, &cache.get(&Matrices { mass: &changed, ..matrices }, 3, &force, 2, &limits)?.data));
        for loss in [false, true] {
            cache = original.clone();
            let matrices = if loss { Matrices { loss: Some(&m), ..matrices } } else { Matrices { viscous: Some(&m), ..matrices } };
            assert!(!Arc::ptr_eq(&cold.data, &cache.get(&matrices, 3, &force, 2, &limits)?.data));
        }
        cache = original.clone();
        assert!(!Arc::ptr_eq(&cold.data, &cache.get(&matrices, 4, &force, 2, &limits)?.data));
        cache = original.clone();
        assert!(!Arc::ptr_eq(&cold.data, &cache.get(&matrices, 3, &force, 1, &limits)?.data));
        cache = original;
        assert!(cache.get(&matrices, 3, &force, 2, &Limits { max_factor_entries: 0, ..limits }).is_err());
        assert!(cache.entry.is_none());
        Ok(())
    }
}
