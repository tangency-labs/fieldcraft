//! The local edge-winding certificate requested by issue #685. Positive tetrahedra round an
//! interior edge fill one turn; round a boundary edge they leave a gap. A positive, overlapping
//! star can cover two turns while its summed volume still equals its winding-weighted boundary
//! volume, so the existing volume check cannot replace this check.
//!
//! These are floating-point measures, as in the volume gate, not replacement orientation or
//! intersection predicates. Dihedral measures follow Shewchuk 2002 / Klingner--Shewchuk 2007
//! ([`super::measure`]); `atan2` is the existing platform-independent implementation. The budget
//! is 1e-8 radians per incident tetrahedron, allowing measurement roundoff plus accumulation.
//! This is an angular numerical tolerance, unrelated to solver or volume tolerances. Interior
//! sums must be within that budget of 2pi; boundary sums must be below 2pi by more than the budget.
//! A boundary gap too small to distinguish from a full turn is conservatively uncertified.
//! This is a declared numerical acceptance policy, not a rigorous error bound for every
//! ill-conditioned tetrahedron; exact/adaptive orientation checks remain in the existing gate.
//!
//! Work/storage are bounded by the input: six edge incidences per tetrahedron, no pairwise
//! intersection search. Canonical tetrahedron order and compensated sums give a fixed reduction
//! order even if the caller permutes tetrahedra or their corner slots. This local certificate
//! does not prove the absence of intersections between nonadjacent parts of the mesh.

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

use fieldcraft_model::progress;

use super::{measure, sorted2};
use crate::{Error, Result};

// The same order as measure::{dihedral_sines, dihedral_cosines}.
const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
const ANGULAR_BUDGET: f64 = 1e-8; // radians per incident tetrahedron

#[derive(Default)]
struct Coverage {
    sum: f64,
    correction: f64,
    count: usize,
}

impl Coverage {
    fn add(&mut self, angle: f64) {
        let y = angle - self.correction;
        let next = self.sum + y;
        self.correction = (next - self.sum) - y;
        self.sum = next;
        self.count += 1;
    }

    fn budget(&self) -> f64 {
        ANGULAR_BUDGET * self.count as f64
    }

    fn passes(&self, boundary: bool) -> bool {
        if boundary { self.sum < TAU - self.budget() } else { (self.sum - TAU).abs() <= self.budget() }
    }
}

/// Only called after the existing element-admissibility and volume checks pass. `boundary`
/// contains exactly the canonical faces held once, already checked by `linear_gate`.
pub(super) fn check(points: &[[f64; 3]], tets: &[[u32; 4]], boundary: &[[u32; 3]]) -> Result<()> {
    let boundary_edges: BTreeSet<_> = boundary.iter().flat_map(|[a, b, c]| [sorted2(*a, *b), sorted2(*b, *c), sorted2(*c, *a)]).collect();
    let mut ordered = tets.to_vec();
    for t in &mut ordered {
        t.sort_unstable();
    }
    ordered.sort_unstable();
    let mut coverage: BTreeMap<(u32, u32), Coverage> = BTreeMap::new();
    for (i, t) in ordered.iter().enumerate() {
        if i % 1024 == 0 {
            progress::check()?;
        }
        let point = |v: u32| points.get(v as usize).copied().ok_or_else(|| Error::Failed("gate: vertex id".into()));
        let x = [point(t[0])?, point(t[1])?, point(t[2])?, point(t[3])?];
        let sines = measure::dihedral_sines(&x);
        let cosines = measure::dihedral_cosines(&x);
        for (k, (a, b)) in EDGES.into_iter().enumerate() {
            let (s, c) = (sines[k].abs(), cosines[k]);
            // The measure helpers define degenerate/nonfinite inputs as zero sine. Never let
            // atan2's defined zero fallback turn an unmeasurable incidence into a certificate.
            if !(s > 0.0 && s.is_finite() && c.is_finite()) {
                return Err(Error::failed(crate::message![
                    format!("mesh validity gate: cannot measure the dihedral at edge ({}, {}) near ", t[a], t[b]),
                    crate::message::point(x[a])
                ]));
            }
            coverage.entry((t[a], t[b])).or_default().add(crate::vec3::atan2(s, c));
        }
    }
    for ((a, b), edge) in coverage {
        let boundary = boundary_edges.contains(&(a, b));
        if !edge.passes(boundary) {
            let at = points.get(a as usize).copied().ok_or_else(|| Error::Failed("gate: vertex id".into()))?;
            let kind = if boundary { "boundary" } else { "interior" };
            let expected = if boundary { "less than one turn, outside the budget" } else { "one turn, within the budget" };
            return Err(Error::failed(crate::message![
                format!("mesh validity gate: {kind} edge ({a}, {b}) near "),
                crate::message::point(at),
                format!(
                    " has angular coverage {:.12} radians from {} tetrahedra (one turn is {TAU:.12}; budget {:.3e} radians): expected {expected}; the edge winding is not certified",
                    edge.sum,
                    edge.count,
                    edge.budget()
                )
            ]));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
