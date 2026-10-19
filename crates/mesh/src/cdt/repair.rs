//! Repairing the tetrahedra that rounding to doubles leaves flat, inverted or inadmissible.
//!
//! The exact tetrahedralization is valid, but its Steiner points are rounded to doubles, which
//! can flip a nearly flat tetrahedron (Diazzi et al. 2023, §5.1), and nearly cospherical input
//! (a tessellated sphere, a cylinder) gives tetrahedra that are positive yet flat to rounding:
//! their normalized Jacobian is below FieldCraft's validity criterion
//! ([`fieldcraft_model::admissible_jacobian`]), so a solver would refuse them. Following §5.1,
//! such tetrahedra are removed by flips that never touch the boundary, on the rounded
//! coordinates and decided by exact `orient3d`: edge removal (which includes the 3-2 and 4-4
//! flips) around each of their edges and face removal (the 2-3 flip and multi-face removal)
//! across each of their faces, each applied only when it improves the worst tetrahedron it
//! changes ([`crate::tet::ops`]). Every boundary edge is held fixed, so every input triangle
//! stays a union of tetrahedron faces. The passes are capped; whatever is left is reported by the
//! validity gate, never hidden.

use fieldcraft_model::{normalized_jacobian, tet4_jacobian};

use super::BoundaryFace;
use crate::tet::complex::{EDGES, NONE, TetComplex, VertexKind};
use crate::tet::ops::{remove_edge_at, remove_face};
use crate::tet::refine::feature_curves;
use crate::{Error, Result};

/// A tetrahedron is repaired when its normalized Jacobian is at most this (a margin of a
/// thousand above the validity criterion, so rounding differences cannot bring it back).
const REPAIR_BELOW: f64 = 1e3 * fieldcraft_model::MIN_NORMALIZED_JACOBIAN; // unit-factors: a multiple of a normalized Jacobian
/// Passes over the tetrahedra still to repair.
const MAX_PASSES: usize = 8;

fn needs_repair(points: &[[f64; 3]], t: &[u32; 4]) -> bool {
    let p = |v: u32| points.get(v as usize).copied().unwrap_or([f64::NAN; 3]);
    !(normalized_jacobian(&tet4_jacobian(&[p(t[0]), p(t[1]), p(t[2]), p(t[3])])) > REPAIR_BELOW)
}

/// Repair the rounded mesh in place; the number of tetrahedra that needed it.
pub(super) fn repair(points: &[[f64; 3]], tets: &mut Vec<[u32; 4]>, boundary: &mut Vec<BoundaryFace>, max_tets: usize) -> Result<usize> {
    let bad = tets.iter().filter(|t| needs_repair(points, t)).count();
    if bad == 0 {
        return Ok(0);
    }
    let tags: Vec<(u32, usize, u32)> = boundary.iter().map(|b| (b.tet, usize::from(b.opposite), b.triangle)).collect();
    let kinds = vec![VertexKind::Corner; points.len()];
    let mut cx = TetComplex::with_invalid(points, kinds, tets, &tags, &[], points.len(), max_tets.max(tets.len() * 2))?;
    // The segments: boundary edges between two different input triangles (never by angle).
    // Edges inside one input triangle (between its pieces) may be flipped; the triangle stays
    // covered.
    let (_, segments) = feature_curves(&cx, -2.0)?;
    for ((a, b), c) in segments {
        cx.set_curve(a, b, Some(c));
    }
    for _ in 0..MAX_PASSES {
        let targets: Vec<(u32, [u32; 4])> = cx.live_tets().filter(|(_, vs)| needs_repair(points, vs)).collect();
        if targets.is_empty() {
            break;
        }
        let mut changed = false;
        for (t, vs) in targets {
            if !cx.is_live(t) || cx.tet(t)? != vs {
                continue;
            }
            for &(i, j) in &EDGES {
                if cx.is_live(t) && cx.tet(t)? == vs && remove_edge_at(&mut cx, t, vs[i], vs[j])? {
                    changed = true;
                    break;
                }
            }
            // A flat tetrahedron lying on one input triangle (two of its faces are pieces of
            // it): delete it; its two other faces cover the same piece of the triangle.
            if cx.is_live(t) && cx.tet(t)? == vs && cx.delete_boundary_tet(t)? {
                changed = true;
                continue;
            }
            if cx.is_live(t) && cx.tet(t)? == vs {
                for k in 0..4 {
                    if remove_face(&mut cx, t, k)? {
                        changed = true;
                        break;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    let c = cx.compact()?;
    if c.source.len() != points.len() || c.source.iter().enumerate().any(|(i, &s)| s as usize != i) {
        return Err(Error::Failed("repair removed a vertex".into()));
    }
    *tets = c.tets;
    *boundary = c.boundary.into_iter().map(|(tet, k, tag)| BoundaryFace { tet, opposite: k as u8, triangle: tag }).collect();
    if boundary.iter().any(|b| b.triangle == NONE) {
        return Err(Error::Failed("repair lost a boundary face's triangle".into()));
    }
    Ok(bad)
}
