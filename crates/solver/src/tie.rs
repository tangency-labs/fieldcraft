//! Ties: non-matching meshes bonded by mortar constraints with dual shape functions.
//!
//! Each node A of a tie's secondary surface follows the main surface under it through the
//! mortar condition ∫ ψ_A (u_s − u_m) dA = 0 per direction, integrated over the overlaps of the
//! secondary faces with the main faces in the reference configuration ([`fieldcraft_fem::mortar`]).
//! The multiplier space uses the dual shape functions of each secondary face, ψ_a = Σ_b A_ab N_b
//! with A = D M⁻¹ (M_ab = ∫ N_a N_b over the face, D = diag ∫ N_a), biorthogonal to the face's
//! own shape functions, ∫ ψ_a N_b = δ_ab ∫ N_b (B. I. Wohlmuth, *A mortar finite element method
//! using dual spaces for the Lagrange multiplier*, SIAM J. Numer. Anal. 38(3):989–1012, 2000;
//! A. Popp, B. I. Wohlmuth, M. W. Gee, W. A. Wall, *Dual quadratic mortar finite element methods
//! for 3D finite deformation contact*, SIAM J. Sci. Comput. 34(4):B421–B446, 2012). The mortar
//! matrix of the secondary side is then diagonal, and the condition makes each secondary node a
//! weighted average of the main nodes under it,
//!
//!   u_A = Σ_J (M_AJ / Σ_K M_AK) u_J,   M_AJ = Σ_overlaps ∫ ψ_A N_J dA,
//!
//! a linear constraint row per direction that the constraint elimination takes like a coupling's
//! ([`crate::constraints`]): the secondary node is the intended dependent. Where a secondary
//! face is wholly covered, Σ_K M_AK = ∫ N_A, the plain mortar condition; dividing by the sum keeps
//! rigid motions exact at the surface's edges, where it is not. A uniform stress crosses a tie on
//! flat faces exactly (the tie patch test), because the dual mortar conditions are satisfied by
//! the linear field on both sides.
//!
//! **Nodes in several ties.** A node on the secondary surfaces of several ties (two secondary
//! surfaces meeting at an edge: a partition's cut faces and a body's outer faces, a chip's
//! underside and its sides) gets one constraint, over the main faces of all of them together, as
//! if they were one tie. Two constraints, one per tie, would nearly repeat each other, and
//! eliminating the second divides by round-off: the solution is then silently wrong (#389).
//! Two ties holding one interface both ways (a node following a node that follows it) are an
//! error naming them, and so are constraints that still nearly repeat each other (a tie and a
//! coupling on one node, longer loops; [`crate::constraints`]).
//!
//! **Pairing.** Each secondary face takes the main faces whose boxes, grown by the position
//! tolerance, overlap its own, and that it overlaps along its normal with a mean gap within the
//! tolerance (by default half the thickness of the main element). A secondary node with no main
//! face under it stays free, with a warning. The constraints are small-displacement (linear)
//! constraints in nonlinear steps too.
//!
//! **Temperature.** The same weights tie each secondary node's temperature to the main nodes'
//! (T_A = Σ_J w_AJ T_J) in heat transfer and in the heat part of coupled temperature–displacement
//! steps: the temperature analogue of the displacement tie, so a linear temperature field crosses
//! a flat tie exactly, and since the constraint is eliminated (no interface resistance), the heat
//! leaving one side enters the other.
//!
//! **Scope.** Linear solid faces (tri3, quad4) in 3D models; explicit dynamics refuses ties by
//! name.

use std::collections::BTreeMap;

use fieldcraft_fem::DofKind;
use fieldcraft_fem::contact::{Aabb, Grid, bounds};
use fieldcraft_fem::mortar::integrate;
use fieldcraft_model::{Model, Resolved};

use crate::constraints::{Row, Source};
use crate::contact::faces;
use crate::system::Mesh;
use crate::{Error, Result};

type Vec3 = fieldcraft_fem::frame::Vec3;

/// Faces per main grid cell budget (the grid's cap is this times the faces).
const GRID_CELLS_PER_FACE: usize = 4;

/// A secondary coefficient below this fraction of its node's total is round-off.
const NEGLIGIBLE: f64 = 1e-12; // unit-factors: relative to the node's total

/// x = M⁻¹ b for a small symmetric positive definite M (Gauss–Jordan with partial pivoting);
/// `None` when it is singular.
fn solve(m: &[[f64; 4]; 4], n: usize, b: [f64; 4]) -> Option<[f64; 4]> {
    let mut a = *m;
    let mut x = b;
    for c in 0..n {
        let p = (c..n).max_by(|i, j| a[*i][c].abs().total_cmp(&a[*j][c].abs()))?;
        if !(a[p][c].abs() > 0.0) {
            return None;
        }
        a.swap(c, p);
        x.swap(c, p);
        for r in 0..n {
            if r != c {
                let f = a[r][c] / a[c][c];
                for k in c..n {
                    a[r][k] -= f * a[c][k];
                }
                x[r] -= f * x[c];
            }
        }
    }
    for c in 0..n {
        x[c] /= a[c][c];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// The dual coefficients A = D M⁻¹ of a face (rows a: ψ_a = Σ_b A_ab N_b).
fn dual(x: &[Vec3]) -> Option<[[f64; 4]; 4]> {
    let n = x.len();
    let mirror: Vec<Vec3> = x.iter().rev().copied().collect();
    let mut mass = [[0.0; 4]; 4];
    integrate(x, &mirror, |p| {
        for a in 0..4 {
            for b in 0..4 {
                mass[a][b] += p.weight * p.secondary[a] * p.secondary[b];
            }
        }
    })?;
    // A M = D, M symmetric: row a of A solves M (row a)ᵀ = D_aa e_a.
    let mut out = [[0.0; 4]; 4];
    for a in 0..n {
        let d: f64 = mass[a].iter().take(n).sum();
        let mut e = [0.0; 4];
        e[a] = d;
        out[a] = solve(&mass, n, e)?;
    }
    Some(out)
}

/// One face of a tie surface: its nodes (slots of the caller's numbering, in the face's order)
/// and the depth of the element under it (the default position tolerance is half of it).
pub(crate) struct TieFace {
    pub nodes: Vec<usize>,
    pub thickness: f64,
}

/// The nodes a tie works on: their coordinates and labels by slot, and the faces of a set
/// (`faces(owner, set)`, refusing what the analysis cannot tie).
pub(crate) struct TieNodes<'a, F> {
    pub xyz: &'a [Vec3],
    pub labels: &'a [u32],
    pub faces: F,
}

/// The constraint rows of the ties, one per secondary node (of one tie or several) and kind in
/// `kinds` (secondary nodes by slot, kinds in order); `dof(slot, kind)` numbers the degrees of
/// freedom. The same mortar weights tie displacements (the translations) and temperatures (one
/// scalar per node). Warnings name the secondary nodes left free and those several ties share.
pub(crate) fn rows_of<'m, F>(
    model: &'m Model,
    nodes: &TieNodes<'_, F>,
    kinds: &[DofKind],
    dof: impl Fn(usize, DofKind) -> Result<usize>,
    warnings: &mut Vec<String>,
) -> Result<Vec<Row<'m>>>
where
    F: Fn(&str, &str) -> Result<Vec<TieFace>>,
{
    // Per secondary node (slot): the ties it is a secondary node of, and its main nodes'
    // weights M_AJ summed over them.
    let mut held: BTreeMap<usize, (Vec<&'m str>, BTreeMap<usize, f64>)> = BTreeMap::new();
    let corners = |f: &TieFace| -> Vec<Vec3> { f.nodes.iter().map(|n| nodes.xyz.get(*n).copied().unwrap_or_default()).collect() };
    let label = |n: usize| nodes.labels.get(n).copied().unwrap_or(0);
    for (name, tie) in model.ties() {
        let owner = format!("tie `{name}`");
        let secondary = (nodes.faces)(&owner, &tie.secondary)?;
        let main = (nodes.faces)(&owner, &tie.main)?;
        let main_x: Vec<Vec<Vec3>> = main.iter().map(corners).collect();
        let reach = |f: &TieFace| tie.position_tolerance.unwrap_or(0.5 * f.thickness);
        let boxes: Vec<Aabb> = main.iter().zip(&main_x).map(|(f, x)| bounds(x, reach(f))).collect();
        let grid = Grid::new(&boxes, GRID_CELLS_PER_FACE * boxes.len() + 64);
        let mut candidates = Vec::new();
        for f in &secondary {
            let xs = corners(f);
            let first = f.nodes.first().map_or(0, |n| label(*n));
            let coefficients = dual(&xs).ok_or_else(|| Error::Invalid(format!("{owner}: a secondary face at node {first} has no area")))?;
            for n in &f.nodes {
                let (ties, _) = held.entry(*n).or_default();
                if ties.last() != Some(&name.as_str()) {
                    ties.push(name.as_str());
                }
            }
            grid.candidates_in(&bounds(&xs, 0.0), &mut candidates);
            for &c in &candidates {
                let (Some(m), Some(xm)) = (main.get(c as usize), main_x.get(c as usize)) else { continue };
                if m.nodes.iter().any(|n| f.nodes.contains(n)) {
                    continue;
                }
                let mut local = [[0.0; 4]; 4];
                let (mut area, mut gap) = (0.0, 0.0);
                integrate(&xs, xm, |p| {
                    area += p.weight;
                    gap += p.weight * p.gap().abs();
                    for a in 0..4 {
                        let psi: f64 = (0..4).map(|b| coefficients[a][b] * p.secondary[b]).sum();
                        for j in 0..4 {
                            local[a][j] += p.weight * psi * p.main[j];
                        }
                    }
                });
                if !(area > 0.0) || gap / area > reach(m) {
                    continue;
                }
                for (a, s) in f.nodes.iter().enumerate() {
                    let (_, row) = held.entry(*s).or_default();
                    for (j, mn) in m.nodes.iter().enumerate() {
                        *row.entry(*mn).or_insert(0.0) += local[a][j];
                    }
                }
            }
        }
    }
    // One constraint per secondary node, over the main faces of all its ties together: two ties
    // holding one node (secondary surfaces meeting at an edge) would otherwise give it two
    // constraints that nearly repeat each other (#389).
    let mut free: BTreeMap<&str, Vec<u32>> = BTreeMap::new();
    let mut shared: BTreeMap<Vec<&str>, Vec<u32>> = BTreeMap::new();
    // Per held node: its ties and the main nodes it follows, with their weights.
    let mut follows: BTreeMap<usize, (&Vec<&'m str>, Vec<(usize, f64)>)> = BTreeMap::new();
    for (node, (ties, row)) in &held {
        let Some(&name) = ties.first() else { continue };
        let total: f64 = row.values().sum();
        if !(total.abs() > 0.0 && total.is_finite()) {
            free.entry(name).or_default().push(label(*node));
            continue;
        }
        if ties.len() > 1 {
            shared.entry(ties.clone()).or_default().push(label(*node));
        }
        let terms: Vec<(usize, f64)> = row.iter().map(|(j, w)| (*j, w / total)).filter(|(_, w)| w.abs() > NEGLIGIBLE).collect();
        follows.insert(*node, (ties, terms));
    }
    // A node may not follow a node that follows it (a tie each way on one interface): the two
    // constraints would hold both twice.
    for (node, (ties, terms)) in &follows {
        for (j, _) in terms {
            let Some((back, theirs)) = follows.get(j) else { continue };
            if theirs.iter().any(|(k, _)| k == node) {
                let names = |t: &Vec<&str>| t.iter().map(|n| format!("tie `{n}`")).collect::<Vec<_>>().join(" and ");
                let at = nodes.xyz.get(*node).copied().unwrap_or_default();
                return Err(Error::Invalid(format!(
                    "{} and {} hold one interface both ways: node {} at ({:.6}, {:.6}, {:.6}) follows node {}, which follows it back, over-constraining both; tie each interface one way only",
                    names(ties),
                    names(back),
                    label(*node),
                    at[0],
                    at[1],
                    at[2],
                    label(*j)
                )));
            }
        }
    }
    let mut out = Vec::new();
    for (node, (ties, terms)) in &follows {
        let Some(&name) = ties.first() else { continue };
        for &kind in kinds {
            let pivot = dof(*node, kind)?;
            let mut row = vec![(pivot, 1.0)];
            for (j, w) in terms {
                row.push((dof(*j, kind)?, -w));
            }
            out.push(Row { terms: row, pivot, coupling: name, source: Source::Tie });
        }
    }
    for (name, list) in &free {
        warnings.push(format!(
            "tie `{name}`: {} secondary node(s) have no main face within the position tolerance and stay free (node {} first)",
            list.len(),
            list.first().copied().unwrap_or(0)
        ));
    }
    for (ties, list) in &shared {
        let names: Vec<String> = ties.iter().map(|t| format!("`{t}`")).collect();
        warnings.push(format!(
            "ties {}: {} node(s) are secondary nodes of each, and follow their main surfaces together (node {} first)",
            names.join(", "),
            list.len(),
            list.first().copied().unwrap_or(0)
        ));
    }
    Ok(out)
}

/// The displacement constraint rows of every tie in a structural analysis (three translations
/// per secondary node).
pub fn rows<'m>(model: &'m Model, resolved: &Resolved, mesh: &Mesh, warnings: &mut Vec<String>) -> Result<Vec<Row<'m>>> {
    if model.ties().is_empty() {
        return Ok(Vec::new());
    }
    if mesh.theory.is_some() {
        let name = model.ties().keys().next().map_or("", String::as_str);
        return Err(Error::NotSupported(format!("ties (tie `{name}`) in 2D models")));
    }
    let nodes = TieNodes {
        xyz: &mesh.xyz,
        labels: &mesh.labels,
        faces: |owner: &str, set: &str| -> Result<Vec<TieFace>> {
            Ok(faces(model, resolved, mesh, owner, set)?.into_iter().map(|f| TieFace { nodes: f.nodes, thickness: f.thickness }).collect())
        },
    };
    let translations = [DofKind::Ux, DofKind::Uy, DofKind::Uz];
    rows_of(model, &nodes, &translations, |slot, kind| mesh.dof(slot, kind), warnings)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dual_functions_are_biorthogonal_to_the_face_s_own() {
        let quad = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.2, 1.0, 0.0], [0.1, 1.3, 0.0]];
        for face in [&quad[..], &quad[..3]] {
            let a = dual(face).unwrap();
            let mirror: Vec<Vec3> = face.iter().rev().copied().collect();
            let mut check = [[0.0; 4]; 4];
            let mut own = [0.0; 4];
            integrate(face, &mirror, |p| {
                for i in 0..4 {
                    own[i] += p.weight * p.secondary[i];
                    let psi: f64 = (0..4).map(|b| a[i][b] * p.secondary[b]).sum();
                    for j in 0..4 {
                        check[i][j] += p.weight * psi * p.secondary[j];
                    }
                }
            });
            for i in 0..face.len() {
                for j in 0..face.len() {
                    let want = if i == j { own[i] } else { 0.0 };
                    assert!((check[i][j] - want).abs() < 1e-12, "{i} {j}: {} vs {want}", check[i][j]);
                }
            }
        }
        assert!(dual(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]).is_none());
    }
}
