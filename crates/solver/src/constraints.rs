//! Couplings of reference points to the mesh, and the elimination of every linear constraint
//! (couplings and prescribed values) from the system.
//!
//! A coupled reference point is a node with six degrees of freedom (three translations, three
//! rotations) after the mesh's translational ones. Each coupling becomes linear constraint rows
//! Σ cⱼ uⱼ = 0:
//!
//! - **Kinematic** (a rigid link): every node s of the set follows the reference point r rigidly
//!   in the coupled translations, u_s = u_r + θ_r × (x_s − x_r) (small rotations; R. D. Cook,
//!   D. S. Malkus, M. E. Plesha, R. J. Witt, *Concepts and Applications of Finite Element
//!   Analysis*, 4th ed., Wiley 2002, §9.2). The reference point's rotations take part through the
//!   lever arms, so a moment on it reaches the face as a couple even when only translations are
//!   coupled. Coupled rotations tie the rotations of shell and beam nodes, θ_s = θ_r; solid
//!   nodes have none.
//! - **Distributing** (an interpolation constraint): the reference point moves with the rigid
//!   motion that fits the set's nodal translations best in the least-squares sense, with equal
//!   weights (as the model and the deck writers define it). About the nodes' centroid c, with
//!   rᵢ = xᵢ − c, the fit is u_c = Σ uᵢ / N and θ = J⁻¹ Σ rᵢ × uᵢ with
//!   J = Σ (|rᵢ|² I − rᵢ rᵢᵀ); the reference point's motion is θ_r = θ and
//!   u_r = u_c + θ × (x_r − c). Only the coupled degrees of freedom are tied. By virtual work a
//!   force or moment on the point is spread over the nodes as the forces whose resultant and
//!   moment about the point equal it, weighted like the fit, and the region stays free to deform
//!   (G. H. K. Heirman, W. Desmet, *Interface reduction of flexible bodies for efficient modeling
//!   of body flexibility in multibody dynamics*, Multibody Syst. Dyn. 24(2):219–234, 2010, on
//!   interpolation constraints). Nodes on one line cannot carry a rotation about that line: that
//!   is an error when rotations are coupled, and their mean motion is used otherwise.
//!
//! **Elimination** (master–slave by transformation, Cook et al. 2002 §9.2; C. A. Felippa,
//! *Introduction to Finite Element Methods*, ch. 8, University of Colorado lecture notes): every
//! degree of freedom is free, prescribed (a boundary condition), dependent on others, or inactive
//! (a reference point direction no coupling ties). Rows are taken one at a time in a fixed order
//! (boundary conditions first, then couplings by name, nodes by label, directions in order). A
//! row is first rewritten in terms of the current free and prescribed degrees of freedom; one
//! free degree of freedom is then chosen to depend on the rest: the row's intended one (the
//! coupled node, or the distributing coupling's reference point) while it is free, otherwise the
//! largest coefficient, lowest index first (the automatic choice of J. F. Abel, M. S. Shephard,
//! *An algorithm for multipoint constraints in finite element analysis*, IJNME 14(3):464–467,
//! 1979). A row with no free degree of freedom left is redundant when the prescribed values
//! satisfy it and contradictory otherwise. Every dependent degree of freedom ends up as a
//! combination of free and prescribed ones, u = T_f q + T_p u_p, so the reduced system is
//! T_fᵀ K T_f q = T_fᵀ (f − K T_p u_p) and the reactions at prescribed degrees of freedom are
//! T_pᵀ (K u − f): no penalty numbers and no Lagrange multipliers, so the matrix stays symmetric
//! positive definite for the Cholesky factorization. With no coupling, T is the identity on the
//! free degrees of freedom and the result is the plain elimination of prescribed values, bit for
//! bit.
//!
//! **Memory.** A held distributing coupling point makes one coupled node depend on all the others,
//! a dense (3N)² block in the reduced matrix; its size is estimated from the rows and refused
//! against the factorization's memory budget before anything is built, and every other reduction
//! counts its entries (stopping at the budget) before it allocates them.

use std::collections::BTreeMap;

use fieldcraft_fem::dense::inverse3;
use fieldcraft_fem::{DofKind, DofLayout, DofSet};
use fieldcraft_model::{BoundaryCondition, CouplingKind, Model, Resolved, Step};

use crate::sparse::Csr;
use crate::system::Mesh;
use crate::{Error, Result};

/// One linear constraint Σ cⱼ uⱼ = 0 and the degree of freedom meant to be eliminated with it.
#[derive(Debug, Clone)]
pub struct Row<'m> {
    pub terms: Vec<(usize, f64)>,
    pub pivot: usize,
    /// The name of the coupling or tie it comes from (for messages).
    pub coupling: &'m str,
    /// What made it.
    pub source: Source,
}

/// What a constraint row comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Coupling,
    /// A symmetry boundary condition on an inclined plane: n · u = 0, and θ ∥ n at nodes with
    /// rotations.
    Symmetry,
    /// A tie: its secondary node follows a weighted average of main nodes.
    Tie,
}

impl Source {
    fn what(self) -> &'static str {
        match self {
            Source::Coupling => "coupling",
            Source::Symmetry => "boundary condition",
            Source::Tie => "tie",
        }
    }
}

/// A coefficient below this fraction of a row's largest one is round-off, not a constraint.
const NEGLIGIBLE: f64 = 1e-10;

/// A tie row whose free coefficients, once the constraints before it are substituted, are all
/// below this fraction of the coefficients that cancelled to give them nearly repeats those
/// constraints: eliminating it would divide by round-off (#389), so it is a conflict.
const NEARLY_REPEATED: f64 = 1e-4;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// [a]×, the matrix of a × ·.
fn cross_matrix(a: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]]
}

fn mat_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

/// The constraint rows of every coupling: couplings by name, nodes by label, directions in order.
pub fn coupling_rows<'m>(model: &'m Model, resolved: &Resolved, mesh: &Mesh) -> Result<Vec<Row<'m>>> {
    let mut rows = Vec::new();
    for (name, c) in model.couplings() {
        let point = mesh.point(&c.reference_point).ok_or_else(|| Error::Internal(format!("reference point `{}`", c.reference_point)))?;
        let xr = mesh.points.get(point).map(|p| p.xyz).ok_or_else(|| Error::Internal("reference point index".into()))?;
        let slot = mesh.layout.point_slot(point);
        // Reference point direction c (0–5), and translation k of mesh node n.
        let rp = |c: usize| DofKind::from_direction(c).ok_or_else(|| Error::Internal("direction".into())).and_then(|kind| mesh.dof(slot, kind));
        let u = |n: usize, k: usize| DofKind::translation(k).ok_or_else(|| Error::Internal("axis".into())).and_then(|kind| mesh.dof(n, kind));
        let set = resolved.sets.get(&c.set).ok_or_else(|| Error::Invalid(format!("coupling `{name}`: set `{}` is not resolved", c.set)))?;
        if set.nodes.is_empty() {
            return Err(Error::Invalid(format!("coupling `{name}`: set `{}` has no nodes", c.set)));
        }
        let nodes: Vec<usize> = set
            .nodes
            .iter()
            .map(|l| mesh.node(*l).ok_or_else(|| Error::Invalid(format!("coupling `{name}`: node {l} belongs to no solid element"))))
            .collect::<Result<_>>()?;
        let xyz = |n: usize| mesh.xyz.get(n).copied().unwrap_or_default();
        // The coupled directions of the reference point, 0–5.
        let coupled = |range: std::ops::Range<usize>| range.filter(|k| c.dofs.contains(*k as u8 + 1));
        match c.kind {
            CouplingKind::Kinematic => {
                let rotations = |n: usize| mesh.layout.kinds(n).includes(DofSet::ROTATIONS);
                if coupled(0..3).next().is_none() && !(nodes.iter().any(|n| rotations(*n)) && coupled(3..6).next().is_some()) {
                    return Err(Error::Invalid(format!(
                        "coupling `{name}` couples only rotations, which the nodes of solid elements do not have; couple translations too"
                    )));
                }
                for &n in &nodes {
                    let r = sub(xyz(n), xr);
                    for k in coupled(0..3) {
                        // u_s,k − u_r,k − (θ × r)_k = 0, (θ × r)_k = θ_a r_b − θ_b r_a.
                        let (a, b) = ((k + 1) % 3, (k + 2) % 3);
                        let pivot = u(n, k)?;
                        let mut terms = vec![(pivot, 1.0), (rp(k)?, -1.0)];
                        if r[b] != 0.0 {
                            terms.push((rp(3 + a)?, -r[b]));
                        }
                        if r[a] != 0.0 {
                            terms.push((rp(3 + b)?, r[a]));
                        }
                        rows.push(Row { terms, pivot, coupling: name, source: Source::Coupling });
                    }
                    // Shell and beam nodes also turn with the point: θ_s = θ_r.
                    if rotations(n) {
                        for k in coupled(3..6) {
                            let pivot =
                                DofKind::from_direction(k).ok_or_else(|| Error::Internal("direction".into())).and_then(|kind| mesh.dof(n, kind))?;
                            rows.push(Row { terms: vec![(pivot, 1.0), (rp(k)?, -1.0)], pivot, coupling: name, source: Source::Coupling });
                        }
                    }
                }
            }
            CouplingKind::Distributing => {
                let count = nodes.len() as f64;
                let mut centroid = [0.0; 3];
                for &n in &nodes {
                    let x = xyz(n);
                    for k in 0..3 {
                        centroid[k] += x[k] / count;
                    }
                }
                let arms: Vec<[f64; 3]> = nodes.iter().map(|n| sub(xyz(*n), centroid)).collect();
                let mut j = [[0.0; 3]; 3];
                for r in &arms {
                    let rr = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    for a in 0..3 {
                        for b in 0..3 {
                            j[a][b] += if a == b { rr } else { 0.0 } - r[a] * r[b];
                        }
                    }
                }
                let trace = j[0][0] + j[1][1] + j[2][2];
                let inverse = inverse3(&j).filter(|(_, det)| *det > 1e-10 * trace.powi(3) / 27.0).map(|(inv, _)| inv);
                if inverse.is_none() && coupled(3..6).next().is_some() {
                    return Err(Error::Invalid(format!(
                        "coupling `{name}`: the nodes of set `{}` lie on one line (or at one point), so a distributing coupling cannot carry rotations through them; couple translations only or use a kinematic coupling",
                        c.set
                    )));
                }
                // θ = Σ Gᵢ uᵢ with Gᵢ = J⁻¹ [rᵢ]×; u_r = Σ (I/N − [d]× Gᵢ) uᵢ with d = x_r − c.
                let d = cross_matrix(sub(xr, centroid));
                let g: Vec<[[f64; 3]; 3]> = arms.iter().map(|r| inverse.map_or([[0.0; 3]; 3], |inv| mat_mul(&inv, &cross_matrix(*r)))).collect();
                for k in coupled(0..6) {
                    let mut terms = Vec::with_capacity(1 + 3 * nodes.len());
                    let pivot = rp(k)?;
                    terms.push((pivot, 1.0));
                    for (&n, gi) in nodes.iter().zip(&g) {
                        for comp in 0..3 {
                            let a = if k < 3 {
                                let dg: f64 = (0..3).map(|m| d[k][m] * gi[m][comp]).sum();
                                if k == comp { 1.0 / count - dg } else { -dg }
                            } else {
                                gi[k - 3][comp]
                            };
                            if a != 0.0 {
                                terms.push((u(n, comp)?, -a));
                            }
                        }
                    }
                    rows.push(Row { terms, pivot, coupling: name, source: Source::Coupling });
                }
            }
        }
    }
    Ok(rows)
}

/// The constraint rows of a step's symmetry conditions on inclined planes (those on a plane
/// normal to an axis hold global degrees of freedom and are prescribed values instead), in the
/// step's order, nodes by label. At every target node the translation along the unit normal n
/// is zero, n · u = 0, eliminating the translation of the largest |nₖ|; at a node with
/// rotations (shells, beams, a reference point whose rotations exist) the rotation about each of
/// two in-plane axes t₁, t₂ is zero too, tᵢ · θ = 0: a mirror image turns only about the normal.
/// `exists` says whether a reference point direction exists (mesh node ones always do).
pub fn symmetry_rows<'m>(model: &Model, resolved: &Resolved, mesh: &Mesh, step: &'m Step, exists: impl Fn(usize) -> bool) -> Result<Vec<Row<'m>>> {
    let mut rows = Vec::new();
    for name in &step.boundary_conditions {
        let Some(BoundaryCondition::Symmetry { target, plane }) = model.boundary_condition(name) else { continue };
        let (Some(n), Some([t1, t2])) = (plane.inclined(), plane.in_plane()) else { continue };
        let rp = target.reference_point();
        if rp.is_some_and(|r| mesh.point(r).is_none()) {
            // Reported by `prescribed` as holding nothing, like any condition on such a point.
            continue;
        }
        for &label in resolved.target_nodes(target) {
            let slot = match rp {
                Some(r) => mesh.point(r).map(|k| mesh.layout.point_slot(k)),
                None => mesh.layout.slot(label),
            };
            let Some(slot) = slot else { continue };
            let kinds = mesh.layout.kinds(slot);
            let dof = |kind: Option<DofKind>| kind.and_then(|k| kinds.contains(k).then(|| mesh.layout.dof(slot, k)).flatten()).filter(|d| exists(*d));
            let row = |axis: [f64; 3], rotation: bool| -> Result<Option<Row<'m>>> {
                let mut terms = Vec::with_capacity(3);
                let mut missing = false;
                for (k, c) in axis.into_iter().enumerate().filter(|(_, c)| *c != 0.0) {
                    match dof(if rotation { DofKind::rotation(k) } else { DofKind::translation(k) }) {
                        Some(d) => terms.push((d, c)),
                        None => missing = true,
                    }
                }
                if missing {
                    // No rotation exists at a solid node: only its translation is held.
                    if rotation && terms.is_empty() {
                        return Ok(None);
                    }
                    return Err(Error::Invalid(format!(
                        "boundary condition `{name}`: node {label} cannot move along the symmetry plane's normal {n:?} (it has no {} along some axis the normal leans on); in a 2D model the plane must contain the Z axis",
                        if rotation { "rotation" } else { "translation" }
                    )));
                }
                let pivot = terms.iter().fold(None, |best: Option<(usize, f64)>, t| match best {
                    Some(b) if b.1.abs() >= t.1.abs() => Some(b),
                    _ => Some(*t),
                });
                Ok(pivot.map(|(pivot, _)| Row { terms, pivot, coupling: name, source: Source::Symmetry }))
            };
            rows.extend(row(n, false)?);
            rows.extend(row(t1, true)?);
            rows.extend(row(t2, true)?);
        }
    }
    Ok(rows)
}

/// The degrees of freedom that exist: every mesh node's, and the reference point directions
/// some constraint row or `extra` (the discrete elements') involves.
fn active_dofs(layout: &DofLayout, extra: impl Iterator<Item = usize>, rows: &[Row]) -> Vec<bool> {
    let mut active = vec![false; layout.dofs()];
    active.iter_mut().take(layout.node_dofs()).for_each(|a| *a = true);
    for dof in extra {
        if let Some(a) = active.get_mut(dof) {
            *a = true;
        }
    }
    for row in rows {
        for (dof, c) in &row.terms {
            if *c != 0.0
                && let Some(a) = active.get_mut(*dof)
            {
                *a = true;
            }
        }
    }
    active
}

/// A sparse list of (index, coefficient) per row.
#[derive(Debug, Clone, Default)]
struct Terms {
    ptr: Vec<usize>,
    idx: Vec<usize>,
    val: Vec<f64>,
}

impl Terms {
    fn row(&self, r: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (a, b) = (self.ptr.get(r).copied().unwrap_or(0), self.ptr.get(r + 1).copied().unwrap_or(0));
        self.idx.get(a..b).unwrap_or_default().iter().copied().zip(self.val.get(a..b).unwrap_or_default().iter().copied())
    }

    fn is_empty_row(&self, r: usize) -> bool {
        self.ptr.get(r) == self.ptr.get(r + 1)
    }

    fn push_row(&mut self, terms: impl Iterator<Item = (usize, f64)>) {
        if self.ptr.is_empty() {
            self.ptr.push(0);
        }
        for (i, v) in terms {
            self.idx.push(i);
            self.val.push(v);
        }
        self.ptr.push(self.idx.len());
    }

    /// The transpose, with `n` rows.
    fn transpose(&self, n: usize) -> Terms {
        let rows = self.ptr.len().saturating_sub(1);
        let mut count = vec![0usize; n + 1];
        for i in &self.idx {
            if let Some(c) = count.get_mut(i + 1) {
                *c += 1;
            }
        }
        for i in 0..n {
            count[i + 1] += count[i];
        }
        let mut next = count.clone();
        let mut idx = vec![0; self.idx.len()];
        let mut val = vec![0.0; self.idx.len()];
        for r in 0..rows {
            for (i, v) in self.row(r) {
                if let Some(slot) = next.get_mut(i)
                    && let (Some(x), Some(y)) = (idx.get_mut(*slot), val.get_mut(*slot))
                {
                    *x = r;
                    *y = v;
                    *slot += 1;
                }
            }
        }
        Terms { ptr: count, idx, val }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum State {
    Free,
    Fixed,
    Dependent(usize),
    Inactive,
}

/// Sort terms by degree of freedom (stably, so equal ones add in the order given) and add the
/// equal ones, dropping exact zeros.
fn combine(mut terms: Vec<(usize, f64)>) -> Vec<(usize, f64)> {
    terms.sort_by_key(|t| t.0);
    let mut out: Vec<(usize, f64)> = Vec::with_capacity(terms.len());
    for (dof, c) in terms {
        match out.last_mut() {
            Some(last) if last.0 == dof => last.1 += c,
            _ => out.push((dof, c)),
        }
    }
    out.retain(|t| t.1 != 0.0);
    out
}

/// The map from the free equations and prescribed values to every degree of freedom.
#[derive(Debug, Clone)]
pub struct Reduction {
    /// The degree of freedom of each free equation.
    pub eq_dof: Vec<usize>,
    /// Whether a degree of freedom exists (the mesh nodes', and the reference point directions
    /// some coupling ties).
    pub active: Vec<bool>,
    /// T_f: per degree of freedom, (equation, coefficient).
    free: Terms,
    /// T_fᵀ: per equation, (degree of freedom, coefficient).
    columns: Terms,
    /// T_p: per degree of freedom, (prescribed degree of freedom, coefficient).
    fixed: Terms,
    /// T_p u_p, the motion of the prescribed values alone.
    prescribed_motion: Vec<f64>,
    /// Whether any degree of freedom depends on others (T is not a plain selection).
    has_dependents: bool,
    /// The degrees of freedom an inclined symmetry plane supports (those of its rows that are not
    /// prescribed): their share of K u − f is the plane's reaction on them, n λ, and is reported
    /// where it acts rather than carried to prescribed degrees of freedom.
    supports: Vec<bool>,
    budget: Budget,
}

/// The memory the reduced matrices may take, checked before they are built: the same budget as
/// the factorization's (`Limits::max_factor_entries` doubles of 8 bytes).
#[derive(Debug, Clone, Copy)]
pub struct Budget {
    pub max_bytes: usize,
    /// Entries of the assembled mesh pattern.
    pub pattern_entries: usize,
    /// Matrices reduced together (stiffness, and mass for frequency steps).
    pub matrices: usize,
}

impl Budget {
    /// Bytes of `entries` reduced entries: one column index and one value per matrix.
    fn bytes(&self, entries: usize) -> usize {
        entries.saturating_mul(8 * (1 + self.matrices))
    }

    fn exceeded(&self, entries: usize) -> bool {
        entries > self.max_bytes / (8 * (1 + self.matrices))
    }
}

/// An empty vector with room for `n` elements, or a limit error when the memory is not there.
fn reserved<T>(n: usize) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| Error::Limit(format!("cannot allocate {} for the reduced matrices", gigabytes(n.saturating_mul(8)))))?;
    Ok(v)
}

fn gigabytes(bytes: usize) -> String {
    format!("{:.2} GB", bytes as f64 / 1e9)
}

impl Reduction {
    /// Eliminate the prescribed values (`prescribed(active)` gives them) and the constraint rows.
    /// `order` is the mesh nodes' order for numbering the free equations; reference point
    /// equations follow them. `layout` numbers the degrees of freedom (mesh nodes, then reference
    /// points) and `extra` lists the reference point ones that exist without a row (discrete
    /// elements); a heat transfer model passes its temperatures and none.
    ///
    /// A distributing coupling whose reference point is held makes one of its nodes depend on
    /// all the others, which couples every coupled degree of freedom with every other: a dense
    /// (3N)² block in the reduced matrix. Its size is checked against `budget` before anything is
    /// built.
    pub fn new(
        layout: &DofLayout,
        extra: impl Iterator<Item = usize>,
        order: &[usize],
        rows: &[Row],
        prescribed: impl FnOnce(&[bool]) -> Result<Vec<Option<f64>>>,
        budget: Budget,
    ) -> Result<Reduction> {
        let n = layout.dofs();
        let mesh_dofs = layout.node_dofs();
        let active = active_dofs(layout, extra, rows);
        let prescribed = prescribed(&active)?;
        let mut state: Vec<State> = (0..n)
            .map(|d| match (active.get(d).copied().unwrap_or(false), prescribed.get(d).copied().flatten()) {
                (false, _) => State::Inactive,
                (true, Some(_)) => State::Fixed,
                (true, None) => State::Free,
            })
            .collect();
        let values: Vec<f64> = (0..n).map(|d| prescribed.get(d).copied().flatten().unwrap_or(0.0)).collect();
        // Dense blocks: the mesh degrees of freedom of the rows meant to eliminate a held
        // reference point, one block per coupling.
        let mut coupled: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
        for row in rows {
            if row.pivot >= mesh_dofs && state.get(row.pivot) == Some(&State::Fixed) {
                coupled.entry(row.coupling).or_default().extend(row.terms.iter().map(|t| t.0).filter(|d| *d < mesh_dofs));
            }
        }
        let mut blocks: Vec<(&str, usize, usize)> = Vec::new();
        for (name, mut dofs) in coupled {
            dofs.sort_unstable();
            dofs.dedup();
            let mut nodes: Vec<usize> = dofs.iter().filter_map(|d| layout.locate(*d).map(|(slot, _)| slot)).collect();
            nodes.dedup();
            blocks.push((name, nodes.len(), dofs.len().saturating_mul(dofs.len())));
        }
        let dense: usize = blocks.iter().fold(0, |a, b| a.saturating_add(b.2));
        if dense > 0 && budget.exceeded(budget.pattern_entries.saturating_add(dense)) {
            let (name, nodes, _) = blocks.iter().max_by_key(|b| b.2).copied().unwrap_or(("?", 0, 0));
            return Err(Error::Limit(format!(
                "distributing coupling `{name}` couples {nodes} nodes and its reference point is held, which ties every coupled node to every other: the reduced stiffness needs about {} (the limit is {}); split the surface into smaller couplings, hold the surface itself, or use a kinematic coupling",
                gigabytes(budget.bytes(budget.pattern_entries.saturating_add(dense))),
                gigabytes(budget.max_bytes)
            )));
        }
        // Dependent expressions in free and prescribed degrees of freedom, and for each free
        // degree of freedom the expressions that use it (possibly listed twice).
        let mut exprs: Vec<Vec<(usize, f64)>> = Vec::new();
        // The row each expression came from (for messages).
        let mut made_by: Vec<(Source, &str)> = Vec::new();
        let mut users: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for row in rows {
            let scale = row.terms.iter().fold(0.0f64, |m, t| m.max(t.1.abs()));
            if !(scale > 0.0 && scale.is_finite()) {
                return Err(Error::Invalid(format!("{} `{}`: a constraint with no finite coefficient", row.source.what(), row.coupling)));
            }
            // The row in current free and prescribed degrees of freedom.
            let mut expanded = Vec::with_capacity(row.terms.len());
            for &(dof, c) in &row.terms {
                match state.get(dof).copied() {
                    Some(State::Free | State::Fixed) => expanded.push((dof, c)),
                    Some(State::Dependent(e)) => expanded.extend(exprs.get(e).into_iter().flatten().map(|t| (t.0, c * t.1))),
                    Some(State::Inactive) | None => return Err(Error::Internal("a constraint on a missing degree of freedom".into())),
                }
            }
            let gross: f64 = expanded.iter().map(|t| t.1.abs()).sum();
            let terms = combine(expanded);
            let candidate = |t: &&(usize, f64)| state.get(t.0) == Some(&State::Free) && t.1.abs() > NEGLIGIBLE * scale;
            let intended = terms.iter().filter(candidate).find(|t| t.0 == row.pivot);
            let largest = terms.iter().filter(candidate).fold(None, |best: Option<&(usize, f64)>, t| match best {
                Some(b) if b.1.abs() >= t.1.abs() => Some(b),
                _ => Some(t),
            });
            // A tie whose secondary node is held already: what is left of the row must be a
            // constraint of its own, not the difference of two near copies of one.
            let nearly_repeated = row.source == Source::Tie && intended.is_none() && largest.is_some_and(|t| t.1.abs() < NEARLY_REPEATED * gross);
            // Held by a prescribed value (the secondary surface meeting a boundary condition or an
            // axis-aligned symmetry plane: a crosspoint), the row can only constrain main degrees
            // of freedom, and nearly repeated, what it adds is the mortar condition of a node
            // whose motion is already given: its weights on main nodes off the boundary do not
            // vanish (the dual basis spreads over the node's faces), so after the substitutions
            // a small but real constraint is left, not round-off (1.3e-6 of the coefficients that
            // cancelled to give it where a tet face meets a symmetry plane, #540; round-off would
            // be ~1e-14 there). Eliminating it divides by that smallness; dropping it leaves the
            // tie off by that fraction at the crosspoint, below the discretization error. So it
            // is dropped, with the same threshold that calls it a near repeat. A secondary degree
            // of freedom another constraint already eliminated (a node secondary in two ties, a
            // tie and a coupling) is a real conflict and still refused (#389).
            if nearly_repeated && state.get(row.pivot) == Some(&State::Fixed) {
                // Still consistent with the prescribed values, or a contradiction.
                let value = |t: &(usize, f64)| t.1 * values.get(t.0).copied().unwrap_or(0.0);
                let residual: f64 = terms.iter().filter(|t| state.get(t.0) == Some(&State::Fixed)).map(value).sum();
                let size: f64 = terms.iter().map(|t| value(t).abs()).sum();
                if residual.abs() > 1e-9 * size.max(f64::MIN_POSITIVE) && residual.abs() > 0.0 {
                    return Err(Error::Invalid(format!(
                        "{} `{}` contradicts the prescribed values (boundary conditions or couplings hold its nodes in a way it cannot follow)",
                        row.source.what(),
                        row.coupling
                    )));
                }
                continue;
            }
            if nearly_repeated {
                let mut others: Vec<String> = row
                    .terms
                    .iter()
                    .filter_map(|t| match state.get(t.0) {
                        Some(State::Dependent(e)) => made_by.get(*e).map(|(k, n)| format!("{} `{n}`", k.what())),
                        _ => None,
                    })
                    .filter(|o| *o != format!("tie `{}`", row.coupling))
                    .collect();
                others.sort_unstable();
                others.dedup();
                return Err(Error::Invalid(format!(
                    "tie `{}` conflicts with {}: they hold the same nodes almost alike (a node secondary in two ties, or two ties on one interface), which over-constrains the model; keep one of them there",
                    row.coupling,
                    if others.is_empty() { "the constraints before it".to_string() } else { others.join(", ") }
                )));
            }
            let Some(&(pivot, cp)) = intended.or(largest) else {
                // Only prescribed degrees of freedom (or round-off) are left.
                let value = |t: &(usize, f64)| t.1 * values.get(t.0).copied().unwrap_or(0.0);
                let residual: f64 = terms.iter().filter(|t| state.get(t.0) == Some(&State::Fixed)).map(value).sum();
                let size: f64 = terms.iter().map(|t| value(t).abs()).sum();
                if residual.abs() > 1e-9 * size.max(f64::MIN_POSITIVE) && residual.abs() > 0.0 {
                    return Err(Error::Invalid(format!(
                        "{} `{}` contradicts the prescribed values (boundary conditions or couplings hold its nodes in a way it cannot follow)",
                        row.source.what(),
                        row.coupling
                    )));
                }
                continue;
            };
            let expr: Vec<(usize, f64)> = terms.iter().filter(|t| t.0 != pivot).map(|t| (t.0, -t.1 / cp)).collect();
            let free_terms: Vec<usize> = expr.iter().filter(|t| state.get(t.0) == Some(&State::Free)).map(|t| t.0).collect();
            // Substitute the pivot in the expressions that use it.
            if let Some(mut using) = users.remove(&pivot) {
                using.sort_unstable();
                using.dedup();
                for u in using {
                    let Some(e) = exprs.get_mut(u) else { continue };
                    let Ok(at) = e.binary_search_by_key(&pivot, |t| t.0) else { continue };
                    let a = e.remove(at).1;
                    let merged = combine(e.drain(..).chain(expr.iter().map(|t| (t.0, a * t.1))).collect());
                    *e = merged;
                    for &dof in &free_terms {
                        users.entry(dof).or_default().push(u);
                    }
                }
            }
            let index = exprs.len();
            for &dof in &free_terms {
                users.entry(dof).or_default().push(index);
            }
            exprs.push(expr);
            made_by.push((row.source, row.coupling));
            if let Some(s) = state.get_mut(pivot) {
                *s = State::Dependent(index);
            }
        }
        // Equations: free mesh degrees of freedom in node order, then reference point ones.
        let mut eq = vec![usize::MAX; n];
        let mut eq_dof = Vec::new();
        let mesh_order = order.iter().flat_map(|node| layout.range(*node));
        for dof in mesh_order.chain(mesh_dofs..n) {
            if state.get(dof) == Some(&State::Free)
                && let Some(e) = eq.get_mut(dof)
                && *e == usize::MAX
            {
                *e = eq_dof.len();
                eq_dof.push(dof);
            }
        }
        if state.iter().enumerate().any(|(d, s)| *s == State::Free && eq.get(d) == Some(&usize::MAX)) {
            return Err(Error::Internal("the node order misses a degree of freedom".into()));
        }
        let (mut free, mut fixed) = (Terms::default(), Terms::default());
        for dof in 0..n {
            match state.get(dof).copied().unwrap_or(State::Inactive) {
                State::Free => {
                    free.push_row(std::iter::once((eq.get(dof).copied().unwrap_or(0), 1.0)));
                    fixed.push_row(std::iter::empty());
                }
                State::Fixed => {
                    free.push_row(std::iter::empty());
                    fixed.push_row(std::iter::once((dof, 1.0)));
                }
                State::Dependent(e) => {
                    let expr = exprs.get(e).map(Vec::as_slice).unwrap_or_default();
                    let state = &state;
                    let on = |want: State| expr.iter().filter(move |t| state.get(t.0) == Some(&want));
                    free.push_row(on(State::Free).map(|t| (eq.get(t.0).copied().unwrap_or(0), t.1)));
                    fixed.push_row(on(State::Fixed).copied());
                }
                State::Inactive => {
                    free.push_row(std::iter::empty());
                    fixed.push_row(std::iter::empty());
                }
            }
        }
        let prescribed_motion = (0..n).map(|d| fixed.row(d).map(|(p, c)| c * values.get(p).copied().unwrap_or(0.0)).sum()).collect();
        let columns = free.transpose(eq_dof.len());
        let mut supports = vec![false; n];
        for row in rows.iter().filter(|r| r.source == Source::Symmetry) {
            for &(dof, _) in &row.terms {
                if state.get(dof) != Some(&State::Fixed)
                    && let Some(x) = supports.get_mut(dof)
                {
                    *x = true;
                }
            }
        }
        Ok(Reduction { eq_dof, active, free, columns, fixed, prescribed_motion, has_dependents: !exprs.is_empty(), supports, budget })
    }

    /// Free equations.
    pub fn neq(&self) -> usize {
        self.eq_dof.len()
    }

    /// T_fᵀ A T_f for each matrix, all over the mesh's degrees of freedom and on one pattern
    /// (reference points carry no stiffness or mass of their own). Row by row with a dense
    /// accumulator (F. G. Gustavson, *Two fast algorithms for sparse matrices*, ACM TOMS
    /// 4(3):250–269, 1978), in a fixed order, sharing the symbolic work between the matrices;
    /// structural zeros are kept, so the results share one pattern too.
    pub fn reduce<const M: usize>(&self, mats: [&Csr; M]) -> Result<[Csr; M]> {
        if let Some(first) = mats.first()
            && mats.iter().any(|m| m.ptr != first.ptr || m.idx != first.idx)
        {
            return Err(Error::Internal("matrices with different patterns".into()));
        }
        let pattern = mats.first().ok_or_else(|| Error::Internal("no matrix to reduce".into()))?;
        let n = self.neq();
        let mut mark = vec![usize::MAX; n];
        // A plain selection has at most the mesh pattern's entries; otherwise count them first
        // (stopping at the budget) so nothing larger than the budget is ever allocated.
        let entries = if self.has_dependents {
            let budget = Budget { matrices: M, ..self.budget };
            let mut count = 0usize;
            for e in 0..n {
                for (dof, _) in self.columns.row(e) {
                    for b in pattern.row(dof).0 {
                        for (f, _) in self.free.row(*b) {
                            if let Some(m) = mark.get_mut(f)
                                && *m != e
                            {
                                *m = e;
                                count += 1;
                            }
                        }
                    }
                }
                if budget.exceeded(count) {
                    return Err(Error::Limit(format!(
                        "the reduced matrices need more than {} (the limit is {}); couplings tie too many nodes together",
                        gigabytes(budget.bytes(count)),
                        gigabytes(budget.max_bytes)
                    )));
                }
            }
            mark.iter_mut().for_each(|m| *m = usize::MAX);
            count
        } else {
            pattern.idx.len()
        };
        let mut acc = vec![[0.0; M]; n];
        let mut touched: Vec<usize> = Vec::new();
        let mut ptr = Vec::with_capacity(n + 1);
        ptr.push(0);
        let mut idx: Vec<usize> = reserved(entries)?;
        let mut vals: [Vec<f64>; M] = std::array::from_fn(|_| Vec::new());
        for v in &mut vals {
            *v = reserved(entries)?;
        }
        for e in 0..n {
            touched.clear();
            for (dof, ca) in self.columns.row(e) {
                let (start, end) = (pattern.ptr.get(dof).copied().unwrap_or(0), pattern.ptr.get(dof + 1).copied().unwrap_or(0));
                for at in start..end {
                    let Some(&b) = pattern.idx.get(at) else { continue };
                    let v: [f64; M] = std::array::from_fn(|i| mats.get(i).and_then(|m| m.val.get(at)).copied().unwrap_or(0.0));
                    for (f, cb) in self.free.row(b) {
                        let (Some(m), Some(x)) = (mark.get_mut(f), acc.get_mut(f)) else { continue };
                        if *m != e {
                            *m = e;
                            *x = [0.0; M];
                            touched.push(f);
                        }
                        for i in 0..M {
                            x[i] += ca * v[i] * cb;
                        }
                    }
                }
            }
            touched.sort_unstable();
            for &f in &touched {
                idx.push(f);
                let x = acc.get(f).copied().unwrap_or([0.0; M]);
                for i in 0..M {
                    vals[i].push(x[i]);
                }
            }
            ptr.push(idx.len());
        }
        Ok(vals.map(|val| {
            let mut out = Csr { n, ptr: ptr.clone(), idx: idx.clone(), val };
            if self.has_dependents {
                out.mirror_lower();
            }
            out
        }))
    }

    /// The reduced right-hand side T_fᵀ (f − K T_p u_p); `k` is over the mesh's degrees of freedom.
    pub fn rhs(&self, f: &[f64], k: &Csr) -> Vec<f64> {
        let mut ku = Vec::new();
        if self.prescribed_motion.iter().any(|v| *v != 0.0) {
            ku = vec![0.0; k.n];
            k.mul(self.prescribed_motion.get(..k.n).unwrap_or_default(), &mut ku);
        }
        let mut rhs = vec![0.0; self.neq()];
        for (dof, fv) in f.iter().enumerate() {
            let g = fv - ku.get(dof).copied().unwrap_or(0.0);
            for (e, c) in self.free.row(dof) {
                if let Some(r) = rhs.get_mut(e) {
                    *r += c * g;
                }
            }
        }
        rhs
    }

    /// Every degree of freedom from the free equations' values `q` (and the prescribed values
    /// when `with_values`).
    pub fn expand(&self, q: &[f64], with_values: bool) -> Vec<f64> {
        (0..self.prescribed_motion.len())
            .map(|d| {
                let free: f64 = self.free.row(d).map(|(e, c)| c * q.get(e).copied().unwrap_or(0.0)).sum();
                match self.prescribed_motion.get(d) {
                    Some(base) if with_values && !self.fixed.is_empty_row(d) => base + free,
                    _ => free,
                }
            })
            .collect()
    }

    /// T_fᵀ r: a vector over every degree of freedom (a residual) restricted to the free
    /// equations.
    pub fn restrict(&self, r: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.neq()];
        for (dof, v) in r.iter().enumerate() {
            for (e, c) in self.free.row(dof) {
                if let Some(x) = out.get_mut(e) {
                    *x += c * v;
                }
            }
        }
        out
    }

    /// Every degree of freedom from the free equations' values `q` with the prescribed values
    /// scaled by `scale` (a load factor): T_f q + scale T_p u_p.
    pub fn expand_scaled(&self, q: &[f64], scale: f64) -> Vec<f64> {
        (0..self.prescribed_motion.len())
            .map(|d| {
                let free: f64 = self.free.row(d).map(|(e, c)| c * q.get(e).copied().unwrap_or(0.0)).sum();
                free + scale * self.prescribed_motion.get(d).copied().unwrap_or(0.0)
            })
            .collect()
    }

    /// The largest prescribed motion, |T_p u_p|∞.
    pub fn largest_prescribed(&self) -> f64 {
        self.prescribed_motion.iter().fold(0.0, |m, v| m.max(v.abs()))
    }

    /// The largest prescribed motion over the degrees of freedom `select` picks (the
    /// translations, or the rotations, which a nonlinear step measures apart).
    pub fn largest_prescribed_where(&self, select: impl Fn(usize) -> bool) -> f64 {
        self.prescribed_motion.iter().enumerate().filter(|(d, _)| select(*d)).fold(0.0, |m, (_, v)| m.max(v.abs()))
    }

    /// T_pᵀ g for a vector g over every degree of freedom (the reactions of the out-of-balance
    /// forces g = f_int − f_ext), zero where nothing is prescribed; at the degrees of freedom an
    /// inclined symmetry plane supports, g itself.
    pub fn reactions_of(&self, g: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.prescribed_motion.len()];
        for (dof, v) in g.iter().enumerate() {
            if self.supports.get(dof).copied().unwrap_or(false) {
                if let Some(x) = out.get_mut(dof) {
                    *x += v;
                }
                continue;
            }
            for (p, c) in self.fixed.row(dof) {
                if let Some(x) = out.get_mut(p) {
                    *x += c * v;
                }
            }
        }
        out
    }

    /// The reactions at the prescribed degrees of freedom, T_pᵀ (K u − f), and at those an
    /// inclined symmetry plane supports ([`Reduction::reactions_of`]), over every degree of
    /// freedom (zero elsewhere).
    pub fn reactions(&self, k: &Csr, u: &[f64], f: &[f64]) -> Vec<f64> {
        let residual: Vec<f64> = f
            .iter()
            .enumerate()
            .map(|(dof, fv)| {
                if self.fixed.is_empty_row(dof) && !self.supports.get(dof).copied().unwrap_or(false) {
                    return 0.0;
                }
                let (cols, vals) = k.row(dof);
                let ku: f64 = cols.iter().zip(vals).map(|(j, v)| v * u.get(*j).copied().unwrap_or(0.0)).sum();
                ku - fv
            })
            .collect();
        self.reactions_of(&residual)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tie(terms: Vec<(usize, f64)>, pivot: usize) -> Row<'static> {
        Row { terms, pivot, coupling: "glue", source: Source::Tie }
    }

    /// Dofs `held` prescribed, dof 1 at `one` and the others at 0.
    fn reduce(rows: &[Row], held: &[usize], one: f64) -> Result<Reduction> {
        let layout = DofLayout::uniform(0..2, DofSet::TRANSLATIONS);
        let budget = Budget { max_bytes: usize::MAX, pattern_entries: 0, matrices: 1 };
        Reduction::new(
            &layout,
            std::iter::empty(),
            &[0, 1],
            rows,
            |active| Ok((0..active.len()).map(|d| held.contains(&d).then_some(if d == 1 { one } else { 0.0 })).collect()),
            budget,
        )
    }

    /// A tie row of a secondary degree of freedom held by a boundary condition (dof 1), which
    /// after substituting an earlier tie row (dof 0 following dofs 2 and 3) nearly repeats it, is
    /// dropped (a crosspoint); dof 3 stays free.
    #[test]
    fn a_tie_row_on_a_held_degree_of_freedom_is_dropped() {
        let rows = [tie(vec![(0, 1.0), (2, -0.5), (3, -0.5)], 0), tie(vec![(1, 1.0), (0, -1.0), (2, 0.5), (3, 0.5 - 1e-7)], 1)];
        let r = reduce(&rows, &[1, 2], 0.0).unwrap();
        assert_eq!(r.neq(), 3, "dofs 3, 4, 5 free");
    }

    /// Nearly repeated but inconsistent: on a secondary degree of freedom another tie already
    /// eliminated (a node secondary in two ties) it is refused (#389), and on a held one whose
    /// prescribed values it contradicts too.
    #[test]
    fn a_near_repeat_on_an_eliminated_degree_of_freedom_is_refused() {
        let rows = [tie(vec![(0, 1.0), (2, -0.5), (3, -0.5)], 0), tie(vec![(0, 1.0), (2, -0.5), (3, -0.5 - 1e-7)], 0)];
        let e = reduce(&rows, &[2], 0.0).err().map(|e| e.to_string()).unwrap_or_default();
        assert!(e.contains("almost alike"), "{e}");
        // Dof 1 held at 1 and dof 2 at 0 while the row says dof 1 follows them.
        let rows = [tie(vec![(0, 1.0), (2, -0.5), (3, -0.5)], 0), tie(vec![(1, 1.0), (0, -1.0), (2, 0.5), (3, 0.5 - 1e-7)], 1)];
        let e = reduce(&rows, &[1, 2], 1.0).err().map(|e| e.to_string()).unwrap_or_default();
        assert!(e.contains("contradicts"), "{e}");
    }

    #[test]
    fn reduced_matrix_budget_refuses_byte_counts_that_overflow() {
        for matrices in [1, 2] {
            let budget = Budget { max_bytes: usize::MAX, pattern_entries: 0, matrices };
            let entries = usize::MAX / (8 * (1 + matrices));
            assert!(!budget.exceeded(entries));
            assert!(budget.exceeded(entries + 1));
        }
    }
}
