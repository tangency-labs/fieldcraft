//! Sparse symmetric indefinite factorization P A Pᵀ = L D Lᵀ: L unit lower triangular, D block
//! diagonal with 1 × 1 and 2 × 2 pivots; supernodal, left-looking between supernodes,
//! multithreaded, with a fixed operation order, on the symbolic factorization of the Cholesky
//! factor ([`crate::cholesky`]).
//!
//! **Pivoting** (J. R. Bunch, L. Kaufman, *Some stable methods for calculating inertia and solving
//! symmetric linear systems*, Math. Comp. 31(137):163–179, 1977), in its bounded form (C.
//! Ashcraft, R. G. Grimes, J. G. Lewis, *Accurate symmetric indefinite linear equation solvers*,
//! SIAM J. Matrix Anal. Appl. 20(2):513–561, 1998, the "bounded Bunch–Kaufman" or rook search):
//! with α = (1 + √17)/8 and ω_i the largest off-diagonal magnitude of column i, column k is a
//! 1 × 1 pivot when |a_kk| ≥ α ω_k; otherwise the search walks from column i to the row r of its
//! largest entry: r is a 1 × 1 pivot when |a_rr| ≥ α ω_r, (i, r) a 2 × 2 pivot when ω_r = ω_i
//! (the entry is the largest of both columns), else the search continues from r. Every entry of
//! L is then bounded by 1/(1 − α) ≈ 2.78 for 1 × 1 pivots (and a like constant for 2 × 2 ones),
//! which Bunch–Kaufman's single look does not guarantee; ω grows strictly, so the search ends.
//! A 2 × 2 pivot has |a_ii a_rr| < α² a_ir², a negative determinant: one positive and one
//! negative eigenvalue.
//!
//! **Where pivots come from: within supernodes** (the supernode Bunch–Kaufman of O. Schenk, K.
//! Gärtner, *On fast factorization pivoting methods for sparse symmetric indefinite systems*,
//! Electron. Trans. Numer. Anal. 23:158–179, 2006). The columns of a supernode are fully summed
//! once its descendants' updates are in (left-looking, as the Cholesky factor), so any of them
//! can be the next pivot: a symmetric interchange inside the supernode leaves the row structure
//! of every supernode unchanged, so the static symbolic factorization, the ordering, the
//! supernodal tasks and the fill are the Cholesky factor's exactly. The pivot search looks at
//! the supernode's own rows; entries in the rows below (later supernodes) are not candidates.
//! Multifrontal codes delay a column whose pivot fails into the parent (I. S. Duff, J. K. Reid,
//! *The multifrontal solution of indefinite sparse symmetric linear equations*, ACM TOMS
//! 9(3):302–325, 1983; J. W. H. Liu, *The multifrontal method for sparse matrix solution: theory
//! and practice*, SIAM Review 34(1):82–109, 1992), which changes the structure as it goes; this
//! factor keeps the static structure instead (deterministic, threads as Cholesky, the same
//! memory bound checked up front), and:
//! - a column with no acceptable pivot in its supernode — every remaining candidate negligible,
//!   at most `PIVOT_TOLERANCE` times the largest entry of the matrix's row — is reported as
//!   singular, naming the column; a nonsingular saddle point whose zero diagonal couples only to
//!   rows of a later supernode would also stop there (none of FieldCraft's systems has one: the
//!   hybrid elements' pressures are condensed per element);
//! - L's entries in the rows below the supernode are not bounded by the search; the stiffness,
//!   shifted and mixed tangents solved here have sizeable diagonals, and
//!   [`Symmetric::solve_refined`] adds iterative refinement against A where accuracy matters.
//!
//! **Numeric.** A supernode starts as its columns of A minus L_K D_K L_Kᵀ of every descendant K
//! that reaches it, in supernode order (the L D rows formed on the fly from L and D, so nothing
//! more is stored than Cholesky's L plus the pivot order). It is then factored right-looking by
//! panels of [`PANEL`] columns: inside a panel each column is brought up to date with the panel's
//! columns before it when the search looks at it (the partial factorization of LAPACK's
//! `dsytrf`, Anderson et al., *LAPACK Users' Guide*, SIAM 1999), and after the panel the rest of
//! the supernode receives its update by the tiled, parallel kernel of the Cholesky factor. A
//! 2 × 2 pivot may end a panel one column late.
//!
//! **Inertia** (Sylvester's law): the counts of positive and negative eigenvalues of A are those
//! of D — a 1 × 1 pivot by its sign, a 2 × 2 one by its eigenvalues — so a factorization of
//! K − σ M counts the eigenvalues below σ, and a tangent's negative pivots count its unstable
//! directions. A pivot that would be zero stops the factorization (above), so the count of zero
//! eigenvalues is never reported as anything but a located error.
//!
//! **Threads and determinism.** Every pivot decision and every entry is computed by one task in
//! an order fixed by the matrix, ties in the search broken by the lowest index, so the factor,
//! its pivots and its inertia are the same bit for bit on any number of threads.

use serde::Serialize;

use crate::Result;
use crate::cholesky::{
    Context, Factor, FactorError, PANEL, PIVOT_TOLERANCE, Panel, Pivot, Structure, Supernode, Symbolic, UpperColumns, offset, products,
    solve_supernodes, subtract, times_d,
};
use crate::sparse::Csr;

/// (1 + √17)/8: the growth bound of Bunch and Kaufman.
const ALPHA: f64 = 0.640_388_203_202_207_6;

/// Most steps of one rook search (each strictly raises ω, so the cap only bounds a pathological
/// sequence of near ties); at the cap the last two columns are the 2 × 2 pivot.
const ROOK_STEPS: usize = 64;

/// The counts of positive and negative eigenvalues of a factored symmetric matrix.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub struct Inertia {
    pub positive: usize,
    pub negative: usize,
}

/// One supernode of an LDLᵀ factor.
#[derive(Debug, Clone, Default)]
struct Block {
    /// L below the pivots and D on them, the lower trapezoid by columns in pivot order.
    values: Vec<f64>,
    pivots: Vec<Pivot>,
    /// `order[c]`: the supernode's own column (0..w, in the symbolic order) at pivot position c.
    order: Vec<usize>,
    inertia: Inertia,
    two_by_two: usize,
}

impl Supernode for Block {
    #[inline]
    fn values(&self) -> &[f64] {
        &self.values
    }
    #[inline]
    fn pivots(&self) -> Option<&[Pivot]> {
        Some(&self.pivots)
    }
    #[inline]
    fn order(&self) -> Option<&[usize]> {
        Some(&self.order)
    }
}

/// The pivot a search chose, at supernode-local positions not yet moved to k.
enum Choice {
    One(usize),
    Two(usize, usize),
}

/// Symmetric interchange of local columns p < q (both among the supernode's own columns) in an
/// m-row lower trapezoid: the rows of the columns before p, the diagonals, the strip between them
/// and the rows below q.
fn interchange(v: &mut [f64], m: usize, p: usize, q: usize) {
    let at = |r: usize, c: usize| offset(m, c) + r - c;
    for c in 0..p {
        v.swap(at(p, c), at(q, c));
    }
    v.swap(at(p, p), at(q, q));
    for c in p + 1..q {
        v.swap(at(c, p), at(q, c));
    }
    for r in q + 1..m {
        v.swap(at(r, p), at(r, q));
    }
}

/// The largest magnitude of `v[i]` over `range` except `skip`, with the lowest such index; NaN
/// when an entry is not a number.
fn largest(v: &[f64], range: std::ops::Range<usize>, skip: Option<usize>) -> (usize, f64) {
    let mut best = (range.start, 0.0f64);
    for i in range {
        if Some(i) == skip {
            continue;
        }
        let x = v.get(i).copied().unwrap_or(0.0).abs();
        if x.is_nan() {
            return (i, f64::NAN);
        }
        if x > best.1 {
            best = (i, x);
        }
    }
    best
}

/// The work buffers of one supernode's panels.
struct Scratch {
    /// The current values of the column being searched and of a candidate, by local row.
    a: Vec<f64>,
    b: Vec<f64>,
    /// A row of L D over the panel's factored columns.
    w: Vec<f64>,
}

impl Block {
    /// Column `c` (local, ≥ k) of the part still to factor, rows k..m, into `out`: the trapezoid's
    /// entries less the products with the panel's factored columns c0..k.
    fn current(&self, m: usize, c0: usize, k: usize, c: usize, out: &mut [f64], w: &mut Vec<f64>) {
        let t = Panel { m, data: &self.values };
        for r in k..m {
            if let Some(o) = out.get_mut(r) {
                *o = if r >= c { t.at(r, c) } else { t.at(c, r) };
            }
        }
        w.clear();
        w.extend((c0..k).map(|q| t.at(c, q)));
        times_d(w, &t, c0, self.pivots.get(c0..k).unwrap_or_default());
        for (q, wq) in (c0..k).zip(w.iter()) {
            for (o, l) in out.get_mut(k..m).unwrap_or_default().iter_mut().zip(t.column(q, k, m - k)) {
                *o -= l * wq;
            }
        }
    }

    /// Swap local columns p and q (pivot positions ≥ k) in the trapezoid, the order and the
    /// current columns.
    fn swap(&mut self, m: usize, p: usize, q: usize, scratch: &mut Scratch) {
        if p == q {
            return;
        }
        let (p, q) = (p.min(q), p.max(q));
        interchange(&mut self.values, m, p, q);
        self.order.swap(p, q);
        scratch.a.swap(p, q);
        scratch.b.swap(p, q);
    }

    /// Store a 1 × 1 pivot at k from the current column `col`.
    fn one(&mut self, m: usize, k: usize, col: &[f64]) {
        let d = col.get(k).copied().unwrap_or(0.0);
        let start = offset(m, k);
        if let Some(x) = self.values.get_mut(start) {
            *x = d;
        }
        for r in k + 1..m {
            if let Some(x) = self.values.get_mut(start + r - k) {
                *x = col.get(r).copied().unwrap_or(0.0) / d;
            }
        }
        if let Some(p) = self.pivots.get_mut(k) {
            *p = Pivot::One;
        }
        if d > 0.0 {
            self.inertia.positive += 1;
        } else {
            self.inertia.negative += 1;
        }
    }

    /// Store a 2 × 2 pivot at k, k + 1 from the current columns `a` (k) and `b` (k + 1).
    fn two(&mut self, m: usize, k: usize, a: &[f64], b: &[f64]) {
        let get = |v: &[f64], r: usize| v.get(r).copied().unwrap_or(0.0);
        let (d11, d21, d22) = (get(a, k), get(a, k + 1), get(b, k + 1));
        let det = d11 * d22 - d21 * d21;
        let (s0, s1) = (offset(m, k), offset(m, k + 1));
        for (at, v) in [(s0, d11), (s0 + 1, d21), (s1, d22)] {
            if let Some(x) = self.values.get_mut(at) {
                *x = v;
            }
        }
        for r in k + 2..m {
            let (x, y) = (get(a, r), get(b, r));
            if let Some(l) = self.values.get_mut(s0 + r - k) {
                *l = (x * d22 - y * d21) / det;
            }
            if let Some(l) = self.values.get_mut(s1 + r - k - 1) {
                *l = (y * d11 - x * d21) / det;
            }
        }
        if let Some(p) = self.pivots.get_mut(k) {
            *p = Pivot::First;
        }
        if let Some(p) = self.pivots.get_mut(k + 1) {
            *p = Pivot::Second;
        }
        // The eigenvalues' signs: one of each when det < 0, else both the sign of the trace.
        if det < 0.0 {
            self.inertia.positive += 1;
            self.inertia.negative += 1;
        } else if d11 + d22 > 0.0 {
            self.inertia.positive += 2;
        } else {
            self.inertia.negative += 2;
        }
        self.two_by_two += 1;
    }

    /// Factor the panel from column c0 of a supernode with m rows and w columns: pivots chosen
    /// among the columns c0..w, each column brought up to date with the panel's columns before
    /// it. Returns where the panel ended. `scale(local)` is the largest entry of the matrix in
    /// that column; `column(local)` its postorder column.
    fn panel(
        &mut self,
        m: usize,
        w: usize,
        c0: usize,
        scale: &dyn Fn(usize) -> f64,
        column: &dyn Fn(usize) -> usize,
        scratch: &mut Scratch,
    ) -> std::result::Result<usize, FactorError> {
        let end = (c0 + PANEL).min(w);
        let mut k = c0;
        while k < end {
            let singular = |b: &Block, at: usize| FactorError::Singular { column: column(b.order.get(at).copied().unwrap_or(at)) };
            let negligible = |b: &Block, at: usize, x: f64| !(x > PIVOT_TOLERANCE * scale(b.order.get(at).copied().unwrap_or(at)));
            let mut a = std::mem::take(&mut scratch.a);
            self.current(m, c0, k, k, &mut a, &mut scratch.w);
            scratch.a = a;
            let akk = scratch.a.get(k).copied().unwrap_or(0.0).abs();
            let (r0, omega_k) = largest(&scratch.a, k + 1..w, None);
            if !akk.is_finite() || !omega_k.is_finite() {
                return Err(singular(self, k));
            }
            let choice = if akk >= ALPHA * omega_k {
                Choice::One(k)
            } else {
                // The rook search: `a` holds column i, `b` the candidate r.
                let (mut i, mut omega_i, mut r) = (k, omega_k, r0);
                let mut step = 0;
                loop {
                    let mut b = std::mem::take(&mut scratch.b);
                    self.current(m, c0, k, r, &mut b, &mut scratch.w);
                    scratch.b = b;
                    let (s, omega_r) = largest(&scratch.b, k..w, Some(r));
                    let arr = scratch.b.get(r).copied().unwrap_or(0.0).abs();
                    if !arr.is_finite() || !omega_r.is_finite() {
                        return Err(singular(self, k));
                    }
                    step += 1;
                    if arr >= ALPHA * omega_r {
                        std::mem::swap(&mut scratch.a, &mut scratch.b);
                        break Choice::One(r);
                    }
                    if omega_r <= omega_i || step >= ROOK_STEPS {
                        break Choice::Two(i, r);
                    }
                    std::mem::swap(&mut scratch.a, &mut scratch.b);
                    (i, omega_i, r) = (r, omega_r, s);
                }
            };
            match choice {
                Choice::One(p) => {
                    self.swap(m, k, p, scratch);
                    if negligible(self, k, scratch.a.get(k).copied().unwrap_or(0.0).abs()) {
                        return Err(singular(self, k));
                    }
                    let a = std::mem::take(&mut scratch.a);
                    self.one(m, k, &a);
                    scratch.a = a;
                    k += 1;
                }
                Choice::Two(i, r) => {
                    // i to k, then r to k + 1 (r was at k if i took its place).
                    self.swap(m, k, i, scratch);
                    let r = if r == k { i } else { r };
                    self.swap(m, k + 1, r, scratch);
                    let offdiagonal = scratch.a.get(k + 1).copied().unwrap_or(0.0).abs();
                    let size = scale(self.order.get(k).copied().unwrap_or(k)).max(scale(self.order.get(k + 1).copied().unwrap_or(k + 1)));
                    if !(offdiagonal > PIVOT_TOLERANCE * size) {
                        return Err(singular(self, k));
                    }
                    let (a, b) = (std::mem::take(&mut scratch.a), std::mem::take(&mut scratch.b));
                    self.two(m, k, &a, &b);
                    (scratch.a, scratch.b) = (a, b);
                    k += 2;
                }
            }
        }
        Ok(k)
    }
}

/// Factor supernode `j` into `own`; its descendants are `done[k - base]`.
fn supernode(context: &Context, j: usize, done: &[Block], base: usize, own: &mut Block) -> std::result::Result<(), FactorError> {
    context.assemble(j, done, base, &mut own.values)?;
    let (m, cols) = (context.s.rows(j).len(), context.s.cols(j));
    let (w, f) = (cols.len(), cols.start);
    own.pivots = vec![Pivot::One; w];
    own.order = (0..w).collect();
    own.inertia = Inertia::default();
    own.two_by_two = 0;
    let scale = |local: usize| context.a.scale.get(f + local).copied().unwrap_or(0.0);
    let column = |local: usize| f + local;
    let mut scratch = Scratch { a: vec![0.0; m], b: vec![0.0; m], w: Vec::with_capacity(PANEL + 1) };
    let mut c0 = 0;
    while c0 < w {
        let c1 = own.panel(m, w, c0, &scale, &column, &mut scratch)?;
        if c1 < w {
            // The rest of the supernode minus the panel's L D Lᵀ.
            let tiles = products(&Panel { m, data: &own.values }, c0..c1, c1, m - c1, w - c1, own.pivots.get(c0..c1));
            subtract(&mut own.values, m, &tiles, m - c1, w - c1, |i| c1 + i);
        }
        c0 = c1;
    }
    Ok(())
}

/// A symmetric indefinite factor P A Pᵀ = L D Lᵀ in supernodes.
#[derive(Debug, Clone)]
pub struct LdltFactor {
    pub n: usize,
    perm: Vec<usize>,
    structure: Structure,
    blocks: Vec<Block>,
    entries: usize,
    inertia: Inertia,
    two_by_two: usize,
}

impl LdltFactor {
    pub fn new(a: &UpperColumns, max_entries: usize) -> std::result::Result<LdltFactor, FactorError> {
        LdltFactor::from_symbolic(&Symbolic::new(a, max_entries)?)
    }

    pub(crate) fn from_symbolic(symbolic: &Symbolic) -> std::result::Result<LdltFactor, FactorError> {
        let mut blocks = symbolic.allocate()?.into_iter().map(|values| Block { values, ..Block::default() }).collect::<Vec<_>>();
        symbolic.factor(&mut blocks, supernode)?;
        let mut inertia = Inertia::default();
        let mut two_by_two = 0;
        for b in &blocks {
            inertia.positive += b.inertia.positive;
            inertia.negative += b.inertia.negative;
            two_by_two += b.two_by_two;
        }
        Ok(LdltFactor {
            n: symbolic.n,
            perm: symbolic.perm.clone(),
            structure: symbolic.structure.clone(),
            blocks,
            entries: symbolic.entries,
            inertia,
            two_by_two,
        })
    }

    /// Entries in L (as the Cholesky factor of the same pattern).
    pub fn entries(&self) -> usize {
        self.entries
    }

    pub fn inertia(&self) -> Inertia {
        self.inertia
    }

    /// The number of 2 × 2 pivots.
    pub fn two_by_two(&self) -> usize {
        self.two_by_two
    }

    /// Solve A x = b in place.
    pub fn solve(&self, b: &mut [f64]) -> Result<()> {
        solve_supernodes(self.n, &self.perm, &self.structure, &self.blocks, b)
    }
}

/// Which factorization a symmetric solve uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Kind {
    /// Cholesky; a matrix that is not positive definite is refused.
    #[default]
    Definite,
    /// LDLᵀ with pivoting.
    Indefinite,
    /// Cholesky, and LDLᵀ when the matrix is not positive definite (one symbolic factorization
    /// serves both).
    Auto,
}

/// A factored symmetric matrix: Cholesky or LDLᵀ.
#[derive(Debug, Clone)]
pub enum Symmetric {
    Cholesky(Factor),
    Indefinite(LdltFactor),
}

impl Symmetric {
    /// Real symmetric factorization under a factor-only byte budget. The caller reserves
    /// model, CSR matrices, results and iteration vectors separately; numeric order is unchanged.
    pub fn with_memory_limit(
        a: &UpperColumns,
        max_entries: usize,
        kind: Kind,
        max_bytes: usize,
    ) -> std::result::Result<(Symmetric, crate::cholesky::FactorMemory), FactorError> {
        let symbolic = Symbolic::new(a, max_entries)?;
        let estimate = symbolic.memory_estimate(a);
        estimate.check(max_bytes)?;
        let factor = match kind {
            Kind::Definite => Factor::from_symbolic(&symbolic).map(Symmetric::Cholesky),
            Kind::Indefinite => LdltFactor::from_symbolic(&symbolic).map(Symmetric::Indefinite),
            Kind::Auto => match Factor::from_symbolic(&symbolic) {
                Ok(f) => Ok(Symmetric::Cholesky(f)),
                Err(FactorError::NotPositiveDefinite { .. }) => LdltFactor::from_symbolic(&symbolic).map(Symmetric::Indefinite),
                Err(e) => Err(e),
            },
        }?;
        Ok((factor, estimate))
    }

    pub fn new(a: &UpperColumns, max_entries: usize, kind: Kind) -> std::result::Result<Symmetric, FactorError> {
        match kind {
            Kind::Definite => Factor::new(a, max_entries).map(Symmetric::Cholesky),
            Kind::Indefinite => LdltFactor::new(a, max_entries).map(Symmetric::Indefinite),
            Kind::Auto => {
                let symbolic = Symbolic::new(a, max_entries)?;
                match Factor::from_symbolic(&symbolic) {
                    Ok(f) => Ok(Symmetric::Cholesky(f)),
                    Err(FactorError::NotPositiveDefinite { .. }) => LdltFactor::from_symbolic(&symbolic).map(Symmetric::Indefinite),
                    Err(e) => Err(e),
                }
            }
        }
    }

    pub fn n(&self) -> usize {
        match self {
            Symmetric::Cholesky(f) => f.n,
            Symmetric::Indefinite(f) => f.n,
        }
    }

    pub fn entries(&self) -> usize {
        match self {
            Symmetric::Cholesky(f) => f.entries(),
            Symmetric::Indefinite(f) => f.entries(),
        }
    }

    /// The inertia: a Cholesky factor's matrix is positive definite.
    pub fn inertia(&self) -> Inertia {
        match self {
            Symmetric::Cholesky(f) => Inertia { positive: f.n, negative: 0 },
            Symmetric::Indefinite(f) => f.inertia(),
        }
    }

    pub fn solve(&self, b: &mut [f64]) -> Result<()> {
        match self {
            Symmetric::Cholesky(f) => f.solve(b),
            Symmetric::Indefinite(f) => f.solve(b),
        }
    }

    /// Solve A x = b in place, then refine with the residual against `a` (the matrix that was
    /// factored, both triangles) while it shrinks, at most `steps` times (iterative refinement
    /// in working precision: N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd
    /// ed., SIAM 2002, §12.2). Cholesky solutions are returned as solved.
    pub fn solve_refined(&self, a: &Csr, b: &mut [f64], steps: usize) -> Result<()> {
        let rhs = b.to_vec();
        self.solve(b)?;
        if matches!(self, Symmetric::Cholesky(_)) {
            return Ok(());
        }
        let norm = |v: &[f64]| v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let mut ax = vec![0.0; a.n];
        let residual = |x: &[f64], ax: &mut Vec<f64>| -> Vec<f64> {
            a.mul(x, ax);
            rhs.iter().zip(ax.iter()).map(|(p, q)| p - q).collect()
        };
        let mut r = residual(b, &mut ax);
        let mut size = norm(&r);
        for _ in 0..steps {
            if !(size > 0.0) {
                break;
            }
            let mut dx = r.clone();
            self.solve(&mut dx)?;
            let trial: Vec<f64> = b.iter().zip(&dx).map(|(x, d)| x + d).collect();
            let next = residual(&trial, &mut ax);
            let next_size = norm(&next);
            if !(next_size < size) {
                break;
            }
            b.copy_from_slice(&trial);
            (r, size) = (next, next_size);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Error;
    use fieldcraft_parallel as parallel;

    fn csr(rows: &[Vec<(usize, f64)>]) -> Csr {
        let mut ptr = vec![0];
        let mut idx = Vec::new();
        let mut val = Vec::new();
        for row in rows {
            let mut row = row.clone();
            row.sort_by_key(|e| e.0);
            for (c, v) in row {
                idx.push(c);
                val.push(v);
            }
            ptr.push(idx.len());
        }
        Csr { n: rows.len(), ptr, idx, val }
    }

    fn residual(a: &Csr, f: &impl Fn(&mut [f64]), seed: f64) -> f64 {
        let n = a.n;
        let x: Vec<f64> = (0..n).map(|i| ((i as f64 + seed) * 0.37).sin()).collect();
        let mut b = vec![0.0; n];
        a.mul(&x, &mut b);
        f(&mut b);
        b.iter().zip(&x).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max)
    }

    /// A 3D grid Laplacian (3 coupled unknowns per point) shifted by σ: indefinite for σ inside
    /// its spectrum, with a known count of eigenvalues below σ.
    fn shifted_grid(n: [usize; 3], sigma: f64) -> (Csr, usize) {
        let [nx, ny, nz] = n;
        let point = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
        let size = 3 * nx * ny * nz;
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); size];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let p = 3 * point(i, j, k);
                    for a in 0..3 {
                        rows[p + a].push((p + a, 6.0 - sigma));
                    }
                    let next = [
                        (i + 1 < nx).then(|| point(i + 1, j, k)),
                        (j + 1 < ny).then(|| point(i, j + 1, k)),
                        (k + 1 < nz).then(|| point(i, j, k + 1)),
                    ];
                    for q in next.into_iter().flatten() {
                        for a in 0..3 {
                            rows[p + a].push((3 * q + a, -1.0));
                            rows[3 * q + a].push((p + a, -1.0));
                        }
                    }
                }
            }
        }
        // Eigenvalues of the Dirichlet grid Laplacian: Σ 2 − 2 cos(π m/(n + 1)) over the axes.
        let axis = |n: usize| (1..=n).map(move |m| 2.0 - 2.0 * (std::f64::consts::PI * m as f64 / (n as f64 + 1.0)).cos());
        let mut below = 0;
        for a in axis(nx) {
            for b in axis(ny) {
                for c in axis(nz) {
                    if a + b + c < sigma {
                        below += 3;
                    }
                }
            }
        }
        (csr(&rows), below)
    }

    #[test]
    fn indefinite_systems_solve_with_their_inertia() {
        for sigma in [0.0, 0.9, 3.3, 7.1] {
            let (a, below) = shifted_grid([7, 6, 5], sigma);
            let upper = UpperColumns::from_symmetric(&a);
            let f = LdltFactor::new(&upper, usize::MAX).unwrap();
            assert_eq!(f.inertia(), Inertia { positive: a.n - below, negative: below }, "σ = {sigma}");
            let err = residual(&a, &|b: &mut [f64]| f.solve(b).unwrap(), sigma);
            assert!(err < 1e-10, "σ = {sigma}: {err}");
        }
        // A positive definite matrix: the same solution as Cholesky, no negative pivot.
        let (a, _) = shifted_grid([5, 5, 5], -0.5);
        let upper = UpperColumns::from_symmetric(&a);
        let f = LdltFactor::new(&upper, usize::MAX).unwrap();
        assert_eq!(f.inertia().negative, 0);
        assert_eq!(f.entries(), Factor::new(&upper, usize::MAX).unwrap().entries(), "the same fill");
        assert!(residual(&a, &|b: &mut [f64]| f.solve(b).unwrap(), 1.0) < 1e-12);
    }

    #[test]
    fn two_by_two_pivots_for_zero_diagonals() {
        // A saddle point [[K, B], [Bᵀ, 0]] inside one dense supernode: zero diagonals need 2 × 2
        // pivots, and the inertia is (positive: K's size, negative: the constraints').
        let (nk, nb) = (8, 3);
        let n = nk + nb;
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for i in 0..nk {
            for j in 0..nk {
                // Small next to the constraints' entries, so the search pairs them.
                let v = if i == j { 0.5 } else { 0.05 / (1.0 + (i + j) as f64) };
                rows[i].push((j, v));
            }
        }
        for c in 0..nb {
            for i in 0..nk {
                let v = ((i * 7 + c * 3) % 5) as f64 - 2.0 + if i == c { 3.0 } else { 0.0 };
                rows[i].push((nk + c, v));
                rows[nk + c].push((i, v));
            }
            // Stored zeros: the constraints share one supernode with K (a zero diagonal whose
            // couplings all lie in later supernodes is the restriction the module documents).
            for d in 0..nb {
                rows[nk + c].push((nk + d, 0.0));
            }
        }
        // The zero-diagonal columns first, so the search must pivot past them.
        let order: Vec<usize> = (nk..n).chain(0..nk).collect();
        let mut permuted: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let inverse: Vec<usize> = (0..n).map(|i| order.iter().position(|&o| o == i).unwrap()).collect();
        for (i, row) in rows.iter().enumerate() {
            permuted[inverse[i]] = row.iter().map(|(c, v)| (inverse[*c], *v)).collect();
        }
        let a = csr(&permuted);
        let f = LdltFactor::new(&UpperColumns::from_symmetric(&a), usize::MAX).unwrap();
        assert_eq!(f.inertia(), Inertia { positive: nk, negative: nb });
        assert!(f.two_by_two() > 0);
        assert!(residual(&a, &|b: &mut [f64]| f.solve(b).unwrap(), 0.5) < 1e-12);
    }

    #[test]
    fn singular_and_hostile_matrices_are_located_errors() {
        // A Laplacian with the constant in its null space.
        let n = 12;
        let rows: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|i: usize| {
                let mut r = vec![];
                let mut d = 0.0;
                for j in [i.wrapping_sub(1), i + 1] {
                    if j < n {
                        r.push((j, -1.0));
                        d += 1.0;
                    }
                }
                r.push((i, d));
                r
            })
            .collect();
        let e = LdltFactor::new(&UpperColumns::from_symmetric(&csr(&rows)), usize::MAX).unwrap_err();
        assert!(matches!(e, FactorError::Singular { column } if column < n), "{e:?}");
        // NaN.
        let mut bad = rows.clone();
        bad[5][0].1 = f64::NAN;
        bad[4].iter_mut().filter(|e| e.0 == 5).for_each(|e| e.1 = f64::NAN);
        assert!(matches!(LdltFactor::new(&UpperColumns::from_symmetric(&csr(&bad)), usize::MAX), Err(FactorError::Singular { .. })));
        // The cap is checked before the numeric factorization; empty and 1 × 1 matrices.
        let (a, _) = shifted_grid([3, 3, 3], 2.0);
        assert!(matches!(LdltFactor::new(&UpperColumns::from_symmetric(&a), 10), Err(FactorError::TooLarge { .. })));
        let empty = Csr { n: 0, ptr: vec![0], idx: vec![], val: vec![] };
        LdltFactor::new(&UpperColumns::from_symmetric(&empty), 10).unwrap().solve(&mut []).unwrap();
        let f = LdltFactor::new(&UpperColumns::from_symmetric(&csr(&[vec![(0, -4.0)]])), 10).unwrap();
        let mut b = [2.0];
        f.solve(&mut b).unwrap();
        assert_eq!((b, f.inertia()), ([-0.5], Inertia { positive: 0, negative: 1 }));
    }

    #[test]
    fn same_factor_on_any_number_of_threads() {
        // Large enough for parallel subtrees and tiles, with wide supernodes spanning panels.
        let (a, below) = shifted_grid([14, 12, 10], 2.7);
        let upper = UpperColumns::from_symmetric(&a);
        let one = parallel::with_threads(1, || LdltFactor::new(&upper, usize::MAX)).unwrap();
        assert_eq!(one.inertia().negative, below);
        let widest = one.blocks.iter().map(|b| b.order.len()).max().unwrap();
        assert!(widest > 2 * PANEL, "a supernode spans panels: {widest}");
        for threads in [2, 7] {
            let other = parallel::with_threads(threads, || LdltFactor::new(&upper, usize::MAX)).unwrap();
            for (x, y) in one.blocks.iter().zip(&other.blocks) {
                assert!(x.values.iter().zip(&y.values).all(|(p, q)| p.to_bits() == q.to_bits()) && x.order == y.order, "{threads} threads");
            }
        }
        let err = residual(&a, &|b: &mut [f64]| one.solve(b).unwrap(), 0.1);
        assert!(err < 1e-9, "{err}");
    }

    #[test]
    fn solve_across_pivoted_supernodes_and_rhs_errors() {
        // Two dense blocks sharing a separator; interleaving their columns forces a global
        // postorder, while the first block needs a local interchange and a 2 × 2 pivot.
        let mut rows = vec![Vec::new(); 7];
        for (i, row) in [[0.0, 0.1, 3.0], [0.1, 5.0, 0.2], [3.0, 0.2, 0.5]].iter().enumerate() {
            rows[i].extend(row.iter().enumerate().map(|(j, v)| (j, *v)));
        }
        for i in 3..6 {
            rows[i].extend((3..6).map(|j| (j, if i == j { 4.0 } else { 0.25 })));
        }
        for i in 0..6 {
            rows[i].push((6, 0.15));
            rows[6].push((i, 0.15));
        }
        rows[6].push((6, 8.0));
        let order = [3, 0, 4, 1, 5, 2, 6];
        let inverse: Vec<_> = (0..7).map(|i| order.iter().position(|&j| j == i).unwrap()).collect();
        let mixed = csr(&order.iter().map(|&i| rows[i].iter().map(|&(j, v)| (inverse[j], v)).collect()).collect::<Vec<_>>());
        let f = LdltFactor::new(&UpperColumns::from_symmetric(&mixed), usize::MAX).unwrap();
        assert!(f.perm.iter().enumerate().any(|(i, &j)| i != j));
        assert!(f.blocks.len() > 1);
        assert!(f.blocks.iter().enumerate().any(|(k, b)| {
            b.pivots.contains(&Pivot::One)
                && b.pivots.contains(&Pivot::First)
                && b.order.iter().enumerate().any(|(i, &j)| i != j)
                && f.structure.rows(k).len() > b.order.len()
        }));
        let cases = [
            csr(&[]),
            csr(&[vec![(0, 4.0)]]),
            csr(&[vec![(0, -4.0)]]),
            mixed,
            shifted_grid([4, 3, 3], -0.5).0,
            shifted_grid([4, 3, 3], 2.7).0,
            shifted_grid([4, 3, 3], 7.1).0,
        ];
        for (case, a) in cases.iter().enumerate() {
            let upper = UpperColumns::from_symmetric(a);
            let ldlt = LdltFactor::new(&upper, usize::MAX).unwrap();
            let chol = Factor::new(&upper, usize::MAX).ok();
            for seed in [0.0, 0.5, 3.25] {
                let expected: Vec<_> = (0..a.n).map(|i| ((i as f64 + seed) * 0.37).sin()).collect();
                let mut rhs = vec![0.0; a.n];
                a.mul(&expected, &mut rhs);
                let mut actual = rhs.clone();
                ldlt.solve(&mut actual).unwrap();
                assert!(actual.iter().zip(&expected).all(|(x, y)| (x - y).abs() < 1e-10), "case {case}");
                let mut repeated = rhs.clone();
                ldlt.solve(&mut repeated).unwrap();
                assert!(actual.iter().zip(&repeated).all(|(x, y)| x.to_bits() == y.to_bits()));
                if let Some(chol) = &chol {
                    let mut actual = rhs;
                    chol.solve(&mut actual).unwrap();
                    assert!(actual.iter().zip(&expected).all(|(x, y)| (x - y).abs() < 1e-10), "case {case}");
                }
            }
            for len in [a.n + 1, a.n.saturating_sub(1)] {
                if len == a.n {
                    continue;
                }
                let before = vec![0.25; len];
                let mut rhs = before.clone();
                assert!(matches!(ldlt.solve(&mut rhs), Err(Error::Internal(_))));
                assert_eq!(rhs, before);
                if let Some(chol) = &chol {
                    assert!(matches!(chol.solve(&mut rhs), Err(Error::Internal(_))));
                    assert_eq!(rhs, before);
                }
            }
        }
    }

    #[test]
    fn auto_selection_and_refinement() {
        let (spd, _) = shifted_grid([4, 4, 4], -1.0);
        let s = Symmetric::new(&UpperColumns::from_symmetric(&spd), usize::MAX, Kind::Auto).unwrap();
        assert!(matches!(s, Symmetric::Cholesky(_)));
        let (ind, below) = shifted_grid([4, 4, 4], 4.4);
        let s = Symmetric::new(&UpperColumns::from_symmetric(&ind), usize::MAX, Kind::Auto).unwrap();
        assert!(matches!(s, Symmetric::Indefinite(_)));
        assert_eq!(s.inertia().negative, below);
        let x: Vec<f64> = (0..ind.n).map(|i| (i as f64).cos()).collect();
        let mut b = vec![0.0; ind.n];
        ind.mul(&x, &mut b);
        s.solve_refined(&ind, &mut b, 3).unwrap();
        assert!(b.iter().zip(&x).all(|(p, q)| (p - q).abs() < 1e-12));
        assert!(matches!(
            Symmetric::new(&UpperColumns::from_symmetric(&ind), usize::MAX, Kind::Definite),
            Err(FactorError::NotPositiveDefinite { .. })
        ));
    }

    /// The stiffness pattern of a hex8 mesh (27 neighbours, 3 unknowns per node) on an n³ grid of
    /// nodes, ordered by nested dissection, diagonally dominant and shifted by σ.
    fn hex_pattern(n: usize, sigma: f64) -> Csr {
        let id = |i: usize, j: usize, k: usize| i + n * (j + n * k);
        let nodes = n * n * n;
        let xyz: Vec<[f64; 3]> = (0..nodes).map(|v| [(v % n) as f64, ((v / n) % n) as f64, (v / (n * n)) as f64]).collect();
        let mut adjacency = vec![Vec::new(); nodes];
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    for (a, b, c) in (0..27).map(|t| (t % 3, (t / 3) % 3, t / 9)) {
                        let (x, y, z) = (i + a, j + b, k + c);
                        if (1..=n).contains(&x) && (1..=n).contains(&y) && (1..=n).contains(&z) {
                            adjacency[id(i, j, k)].push(id(x - 1, y - 1, z - 1));
                        }
                    }
                }
            }
        }
        let order = crate::ordering::nested_dissection(&xyz, &adjacency);
        let mut place = vec![0; nodes];
        for (p, v) in order.iter().enumerate() {
            place[*v] = p;
        }
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); 3 * nodes];
        for (v, adj) in adjacency.iter().enumerate() {
            for w in adj {
                for a in 0..3 {
                    for b in 0..3 {
                        let value = if v == *w {
                            if a == b { 30.0 - sigma } else { 0.5 }
                        } else if a == b {
                            -1.0
                        } else {
                            0.01 * ((v + w) % 3) as f64
                        };
                        rows[3 * place[v] + a].push((3 * place[*w] + b, value));
                    }
                }
            }
        }
        csr(&rows)
    }

    /// Factor time and fill of Cholesky and LDLᵀ on a large stiffness pattern (`--nocapture`).
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn speed_against_cholesky() {
        for n in [24, 36] {
            let a = hex_pattern(n, 0.0);
            let upper = UpperColumns::from_symmetric(&a);
            let t = std::time::Instant::now();
            let chol = Factor::new(&upper, usize::MAX).unwrap();
            let tc = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let ldlt = LdltFactor::new(&upper, usize::MAX).unwrap();
            let tl = t.elapsed().as_secs_f64();
            assert_eq!((ldlt.entries(), ldlt.inertia().negative), (chol.entries(), 0));
            let shifted = hex_pattern(n, 22.0);
            let t = std::time::Instant::now();
            let indefinite = LdltFactor::new(&UpperColumns::from_symmetric(&shifted), usize::MAX).unwrap();
            let ti = t.elapsed().as_secs_f64();
            assert!(residual(&shifted, &|b: &mut [f64]| indefinite.solve(b).unwrap(), 0.3) < 1e-8);
            eprintln!(
                "{} equations, {} entries: Cholesky {tc:.2} s, LDLᵀ {tl:.2} s (×{:.2}); indefinite (σ = 22: {} negative, {} 2 × 2 pivots) {ti:.2} s",
                a.n,
                chol.entries(),
                tl / tc,
                indefinite.inertia().negative,
                indefinite.two_by_two()
            );
        }
    }
}
