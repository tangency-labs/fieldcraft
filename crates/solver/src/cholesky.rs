//! Sparse Cholesky factorization A = L Lᵀ of a symmetric positive definite matrix: supernodal,
//! left-looking, multithreaded, with a fixed operation order.
//!
//! **Symbolic.** The elimination tree (J. W. H. Liu, *The role of elimination trees in sparse
//! factorization*, SIAM J. Matrix Anal. Appl. 11(1):134–172, 1990) and the column counts of L
//! come from the row subtrees (T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM 2006,
//! §4.1–4.4), so the factor's size is checked against a cap before any memory is committed. The
//! columns are renumbered in a postorder of the tree (an equivalent ordering: the same fill,
//! Liu 1990), so every subtree is a contiguous range of columns, and consecutive columns with
//! nested structures are merged into **fundamental supernodes** (J. W. H. Liu, E. G. Ng, B. W.
//! Peyton, *On finding supernodes for sparse matrix computations*, SIAM J. Matrix Anal. Appl.
//! 14(1):242–252, 1993): a supernode's columns share one row structure and are stored as one
//! dense lower trapezoid, so the factor holds exactly the entries of L.
//!
//! **Numeric** (left-looking supernodal Cholesky: E. Ng, B. W. Peyton, *Block sparse Cholesky
//! algorithms on advanced uniprocessor computers*, SIAM J. Sci. Comput. 14(5):1034–1056, 1993;
//! E. Rothberg, A. Gupta, *An efficient block-oriented approach to parallel sparse Cholesky
//! factorization*, SIAM J. Sci. Comput. 15(6):1413–1439, 1994): a supernode starts as its
//! columns of A, receives the update L_K L_Kᵀ of every descendant supernode K whose rows reach
//! it (K in increasing order), and is then factored densely by panels of [`PANEL`] columns: a
//! panel is updated with the supernode's columns before it, its diagonal block is factored column
//! by column, and the rows below are solved against that block.
//!
//! **Threads and determinism.** Each entry of L is computed by one task, by a sequence of
//! operations fixed by the matrix and the constants [`TILE`] and [`PANEL`]: sums over a
//! supernode's columns run in column order from zero, and updates are subtracted in supernode
//! order. Disjoint subtrees of the supernodal tree are factored concurrently (a supernode reads
//! only its descendants), and the dense work of a large supernode is split into row tiles;
//! neither changes any entry's sequence of operations, so the factor is the same bit for bit on
//! any number of threads and on every platform (multiplications and additions are never fused).
//!
//! A pivot at or below `PIVOT_TOLERANCE` times the matrix's own diagonal entry means the matrix
//! is singular to working precision (an unconstrained rigid-body motion, a disconnected part, a
//! zero-stiffness direction); the factorization stops and names the column: the first such column
//! in the postorder, whichever thread meets it.

use fieldcraft_parallel as parallel;

#[cfg(test)]
mod diagnostics;
#[cfg(test)]
mod path_experiment;

use crate::{Error, Result};

/// Relative pivot size below which the matrix counts as singular.
pub const PIVOT_TOLERANCE: f64 = 1e-10;

/// Rows per tile of the dense kernels (a tile's partial sums stay in registers and L1).
pub const TILE: usize = 64;

/// Columns per panel of a supernode's dense factorization.
pub const PANEL: usize = 64;

/// Rows of a micro-tile of [`products`] (a multiple dividing [`TILE`]).
const MR: usize = 8;

/// Columns of a micro-tile of [`products`]: MR × NR sums stay in registers.
const NR: usize = 4;

/// Dense work (multiply–adds) from which a kernel's tiles run in parallel.
const PARALLEL_WORK: usize = 1 << 21;

/// Subtrees are handed to threads in pieces of about this fraction of the total work.
const SUBTREE_SHARE: usize = 256;

/// The upper triangle of a symmetric matrix by columns: column k holds A(i, k) for i ≤ k, rows
/// sorted.
#[derive(Debug, Clone)]
pub struct UpperColumns {
    pub n: usize,
    pub ptr: Vec<usize>,
    pub idx: Vec<usize>,
    pub val: Vec<f64>,
}

impl UpperColumns {
    /// From a symmetric matrix stored by full rows (row k's entries with column ≤ k are column k
    /// of the upper triangle).
    pub fn from_symmetric(a: &crate::sparse::Csr) -> UpperColumns {
        let mut ptr = Vec::with_capacity(a.n + 1);
        ptr.push(0);
        let mut idx = Vec::new();
        let mut val = Vec::new();
        for k in 0..a.n {
            let (cols, vals) = a.row(k);
            for (c, v) in cols.iter().zip(vals) {
                if *c <= k {
                    idx.push(*c);
                    val.push(*v);
                }
            }
            ptr.push(idx.len());
        }
        UpperColumns { n: a.n, ptr, idx, val }
    }

    fn column(&self, k: usize) -> (&[usize], &[f64]) {
        let (a, b) = (self.ptr.get(k).copied().unwrap_or(0), self.ptr.get(k + 1).copied().unwrap_or(0));
        (self.idx.get(a..b).unwrap_or_default(), self.val.get(a..b).unwrap_or_default())
    }
}

const NONE: usize = usize::MAX;

/// The elimination tree: parent of each column, `NONE` at roots.
fn etree(a: &UpperColumns) -> Vec<usize> {
    let mut parent = vec![NONE; a.n];
    let mut ancestor = vec![NONE; a.n];
    for k in 0..a.n {
        let (rows, _) = a.column(k);
        for &i0 in rows {
            let mut i = i0;
            // Walk from i up to the root of its current subtree, compressing the path to k.
            while i != NONE && i < k {
                let next = ancestor.get(i).copied().unwrap_or(NONE);
                if let Some(x) = ancestor.get_mut(i) {
                    *x = k;
                }
                if next == NONE
                    && let Some(p) = parent.get_mut(i)
                {
                    *p = k;
                }
                i = next;
            }
        }
    }
    parent
}

/// The pattern of row k of L (excluding the diagonal): the reach of column k's entries in the
/// elimination tree. Fills `stack[top..]` and returns `top`; `mark` holds `k` for visited
/// columns.
fn ereach(a: &UpperColumns, k: usize, parent: &[usize], stack: &mut [usize], mark: &mut [usize], path: &mut Vec<usize>) -> usize {
    let n = a.n;
    let mut top = n;
    if let Some(m) = mark.get_mut(k) {
        *m = k;
    }
    let (rows, _) = a.column(k);
    for &i0 in rows {
        let mut i = i0;
        if i > k {
            continue;
        }
        path.clear();
        while i != NONE && mark.get(i).copied() != Some(k) {
            path.push(i);
            if let Some(m) = mark.get_mut(i) {
                *m = k;
            }
            i = parent.get(i).copied().unwrap_or(NONE);
        }
        while let Some(v) = path.pop() {
            top -= 1;
            if let Some(s) = stack.get_mut(top) {
                *s = v;
            }
        }
    }
    top
}

/// The entries of the Cholesky factor of `a`'s pattern, in its own order, without values.
pub(crate) fn factor_entries(a: &UpperColumns) -> usize {
    symbolic(a).1.iter().sum()
}

/// The elimination tree and the entry count of each column of L, without computing values.
fn symbolic(a: &UpperColumns) -> (Vec<usize>, Vec<usize>) {
    let parent = etree(a);
    let mut counts = vec![1usize; a.n];
    let mut stack = vec![0; a.n];
    let mut mark = vec![NONE; a.n];
    let mut path = Vec::new();
    for k in 0..a.n {
        let top = ereach(a, k, &parent, &mut stack, &mut mark, &mut path);
        for &j in stack.get(top..).unwrap_or_default() {
            if let Some(c) = counts.get_mut(j) {
                *c += 1;
            }
        }
    }
    (parent, counts)
}

/// A postorder of a forest given by parents, children in increasing order: `post[k]` is the
/// k-th node.
fn postorder(parent: &[usize]) -> Vec<usize> {
    let n = parent.len();
    // Children as linked lists in increasing order (built from the largest child down).
    let mut head = vec![NONE; n];
    let mut next = vec![NONE; n];
    for j in (0..n).rev() {
        let Some(&p) = parent.get(j) else { continue };
        if let (Some(h), Some(x)) = (head.get(p).copied(), next.get_mut(j)) {
            *x = h;
            if let Some(h) = head.get_mut(p) {
                *h = j;
            }
        }
    }
    let mut post = Vec::with_capacity(n);
    let mut stack = Vec::new();
    for root in (0..n).filter(|j| parent.get(*j) == Some(&NONE)) {
        stack.push(root);
        while let Some(&top) = stack.last() {
            let child = head.get(top).copied().unwrap_or(NONE);
            if child == NONE {
                post.push(top);
                stack.pop();
            } else {
                // Detach the first remaining child and descend into it.
                if let Some(h) = head.get_mut(top) {
                    *h = next.get(child).copied().unwrap_or(NONE);
                }
                stack.push(child);
            }
        }
    }
    post
}

/// Why a factorization stopped.
#[derive(Debug, Clone, PartialEq)]
pub enum FactorError {
    /// The pivot of this column is not positive (relative to the diagonal).
    NotPositiveDefinite { column: usize },
    /// The symmetric indefinite factorization ([`crate::ldlt`]) found no pivot in this column's
    /// supernode that is not negligible next to the column's entries: the matrix is singular to
    /// working precision there.
    Singular { column: usize },
    /// The factor would hold more entries than allowed.
    TooLarge { entries: usize, limit: usize },
    /// A symbolic storage estimate exceeds the caller's factor-only byte budget.
    MemoryLimit { entries: usize, estimated_bytes: usize, limit_bytes: usize },
    /// The memory for the factor could not be had.
    OutOfMemory { entries: usize },
    /// A broken internal invariant (a bug), reported instead of a panic.
    Internal(String),
    /// The command was cancelled ([`fieldcraft_model::progress`]) between supernodes.
    Cancelled,
}

/// Storage estimate for a real supernodal factorization, derived from its actual symbolic
/// structure (Liu–Ng–Peyton 1993, as cited above) and this implementation's buffers.
/// Includes the supplied upper triangle, reordered matrix, indices, numeric factor and
/// concurrent subtree work. Excludes the caller's CSR matrices, model, results and Lanczos
/// basis: callers must reserve those separately before assigning a factor budget.
/// This is an allocation estimate, not a measurement of process RSS or allocator overhead.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FactorMemory {
    pub entries: usize,
    pub values_bytes: usize,
    pub indices_bytes: usize,
    pub matrix_bytes: usize,
    pub workspace_bytes: usize,
}

impl FactorMemory {
    /// Saturation is deliberately a refusal, never a wrapped small allocation estimate.
    pub fn total_bytes(&self) -> usize {
        self.values_bytes.saturating_add(self.indices_bytes).saturating_add(self.matrix_bytes).saturating_add(self.workspace_bytes)
    }

    pub(crate) fn check(&self, limit_bytes: usize) -> std::result::Result<(), FactorError> {
        let estimated_bytes = self.total_bytes();
        if estimated_bytes == usize::MAX || estimated_bytes > limit_bytes {
            return Err(FactorError::MemoryLimit { entries: self.entries, estimated_bytes, limit_bytes });
        }
        Ok(())
    }
}

fn buffer_bytes<T>(v: &Vec<T>) -> usize {
    v.capacity().saturating_mul(std::mem::size_of::<T>())
}

/// Upper bound on simultaneous products buffers, including every tile's packed input.
/// Counting all tiles also covers nested parallel regions without assuming a thread count.
fn products_bytes(mu: usize, wu: usize, kw: usize) -> usize {
    let tiles = mu.div_ceil(TILE);
    let y = wu.saturating_add(wu.div_ceil(NR).saturating_mul(NR)).saturating_mul(kw);
    let tiles_data = tiles.saturating_mul(TILE).saturating_mul(kw.saturating_add(wu));
    y.saturating_add(tiles_data).saturating_mul(8).saturating_add(tiles.saturating_mul(std::mem::size_of::<Vec<f64>>()))
}

impl FactorError {
    /// The column a numeric failure names.
    fn column(&self) -> Option<usize> {
        match self {
            FactorError::NotPositiveDefinite { column } | FactorError::Singular { column } => Some(*column),
            _ => None,
        }
    }

    /// The same failure at another column.
    fn at(self, column: usize) -> FactorError {
        match self {
            FactorError::NotPositiveDefinite { .. } => FactorError::NotPositiveDefinite { column },
            FactorError::Singular { .. } => FactorError::Singular { column },
            e => e,
        }
    }
}

/// One supernode's update to a later one: supernode `k`'s rows from position `p1` on reach the
/// later one, whose columns they meet at positions `p1..p2`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Update {
    pub(crate) k: usize,
    pub(crate) p1: usize,
    pub(crate) p2: usize,
}

/// The supernodal structure (postorder numbering).
#[derive(Debug, Clone)]
pub(crate) struct Structure {
    /// First column of each supernode, and `n` at the end.
    first: Vec<usize>,
    /// Row structure of each supernode (its own columns first, then the rows below), sorted.
    rptr: Vec<usize>,
    rows: Vec<usize>,
}

impl Structure {
    pub(crate) fn cols(&self, s: usize) -> std::ops::Range<usize> {
        self.first.get(s).copied().unwrap_or(0)..self.first.get(s + 1).copied().unwrap_or(0)
    }

    pub(crate) fn rows(&self, s: usize) -> &[usize] {
        let (a, b) = (self.rptr.get(s).copied().unwrap_or(0), self.rptr.get(s + 1).copied().unwrap_or(0));
        self.rows.get(a..b).unwrap_or_default()
    }

    /// Entries of the supernode's lower trapezoid.
    pub(crate) fn size(&self, s: usize) -> usize {
        trapezoid(self.rows(s).len(), self.cols(s).len())
    }
}

/// Entries of an m-row, w-column lower trapezoid (column c holds rows c..m).
pub(crate) fn trapezoid(m: usize, w: usize) -> usize {
    (0..w).map(|c| m.saturating_sub(c)).sum()
}

/// Where column c of an m-row trapezoid starts.
pub(crate) fn offset(m: usize, c: usize) -> usize {
    c * m - c * c.saturating_sub(1) / 2
}

/// A dense lower trapezoid of m rows: one supernode's values, column by column, column c
/// holding rows c..m of the supernode's row structure.
pub(crate) struct Panel<'a> {
    pub(crate) m: usize,
    pub(crate) data: &'a [f64],
}

impl Panel<'_> {
    /// L(r, c) for r ≥ c.
    pub(crate) fn at(&self, r: usize, c: usize) -> f64 {
        self.data.get(offset(self.m, c) + r - c).copied().unwrap_or(0.0)
    }

    /// Rows r..r+len of column c (r ≥ c).
    pub(crate) fn column(&self, c: usize, r: usize, len: usize) -> &[f64] {
        let start = offset(self.m, c) + r - c;
        self.data.get(start..start + len).unwrap_or_default()
    }
}

/// The pivots of a symmetric indefinite factor's block diagonal D, per column of a supernode
/// ([`crate::ldlt`]): a 1 × 1 pivot, or the first or second column of a 2 × 2 one. D's entries
/// sit on the trapezoid's diagonal, and a 2 × 2 pivot's off-diagonal entry below its first
/// column's diagonal (where L's unit block has a zero).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Pivot {
    One,
    First,
    Second,
}

/// `row` (entries for the columns k0.. of a trapezoid) ← `row` D, with D the block diagonal the
/// trapezoid holds for those columns and `pivots` its pivots from column k0.
pub(crate) fn times_d(row: &mut [f64], l: &Panel, k0: usize, pivots: &[Pivot]) {
    let mut k = 0;
    while k < row.len() {
        let d = |r: usize, c: usize| l.at(k0 + r, k0 + c);
        if pivots.get(k) == Some(&Pivot::First) {
            let (a, b) = (row.get(k).copied().unwrap_or(0.0), row.get(k + 1).copied().unwrap_or(0.0));
            let (d11, d21, d22) = (d(k, k), d(k + 1, k), d(k + 1, k + 1));
            if let Some(x) = row.get_mut(k) {
                *x = a * d11 + b * d21;
            }
            if let Some(x) = row.get_mut(k + 1) {
                *x = a * d21 + b * d22;
            }
            k += 2;
        } else {
            if let Some(x) = row.get_mut(k) {
                *x *= d(k, k);
            }
            k += 1;
        }
    }
}

/// S(i, c) = Σ_k L(p + i, k) Y(p + c, k) over the columns k in `cols` of a trapezoid, for rows i
/// in 0..mu and columns c in 0..wu (meaningful for i ≥ c), each sum in column order. Y is L
/// (Cholesky), or L D with D the block diagonal the trapezoid holds when `pivots` (of the
/// columns `cols`) are given. Returned by row tiles: tile t holds rows t·TILE.. as
/// `tile[c * TILE + i]`.
pub(crate) fn products(l: &Panel, cols: std::ops::Range<usize>, p: usize, mu: usize, wu: usize, pivots: Option<&[Pivot]>) -> Vec<Vec<f64>> {
    let (k0, kw) = (cols.start, cols.len());
    let eligible = pivots.is_none() && kw >= 256 && wu >= 128 && mu.saturating_mul(wu).saturating_mul(kw) >= PARALLEL_WORK;
    #[cfg(test)]
    if eligible && std::env::var("FIELDCRAFT_DENSE_UPDATE").is_ok_and(|v| v == "cached") {
        let _timing = crate::timing::Span::new("dense_update_cached");
        crate::timing::dense_work(true, mu, wu, kw, TILE, NR);
        return products_cached(l, cols, p, mu, wu);
    }
    let _timing = eligible.then(|| crate::timing::Span::new("dense_update_reference"));
    if eligible {
        crate::timing::dense_work(false, mu, wu, kw, TILE, NR);
    }
    // Rows p..p+wu of the columns, transposed so that column c's coefficients are contiguous.
    let mut yt = vec![0.0; wu * kw];
    for k in 0..kw {
        for (c, y) in l.column(k0 + k, p, wu).iter().enumerate() {
            if let Some(slot) = yt.get_mut(c * kw + k) {
                *slot = *y;
            }
        }
    }
    if let Some(pivots) = pivots {
        for row in yt.chunks_mut(kw.max(1)) {
            times_d(row, l, k0, pivots);
        }
    }
    // Y by blocks of NR columns, k-major (zero past wu): what a micro-tile reads at each k.
    let mut yp = vec![[0.0f64; NR]; wu.div_ceil(NR) * kw];
    for (c, row) in yt.chunks(kw.max(1)).enumerate().take(wu) {
        for (k, y) in row.iter().enumerate() {
            if let Some(slot) = yp.get_mut((c / NR) * kw + k).and_then(|b| b.get_mut(c % NR)) {
                *slot = *y;
            }
        }
    }
    let tiles = mu.div_ceil(TILE);
    let tile = |t: usize| -> Vec<f64> {
        let i0 = t * TILE;
        let h = TILE.min(mu - i0);
        // The tile's rows of every column by blocks of MR rows, k-major, padded to TILE rows.
        let mut xp = vec![[0.0f64; MR]; TILE / MR * kw];
        for k in 0..kw {
            for (i, x) in l.column(k0 + k, p + i0, h).iter().enumerate() {
                if let Some(slot) = xp.get_mut((i / MR) * kw + k).and_then(|b| b.get_mut(i % MR)) {
                    *slot = *x;
                }
            }
        }
        let cols = wu.min(i0 + h);
        let mut out = vec![0.0; cols * TILE];
        // Micro-tiles of MR rows and NR columns: each S(i, c) is summed in a register over k in
        // order from 0, the same operations as one column at a time, with each X and Y value
        // loaded once per micro-tile.
        for (b, ys) in yp.chunks(kw.max(1)).enumerate().take(cols.div_ceil(NR)) {
            for (r, xs) in xp.chunks(kw.max(1)).enumerate() {
                let mut acc = [[0.0f64; MR]; NR];
                for (x, y) in xs.iter().zip(ys) {
                    for (a, yj) in acc.iter_mut().zip(y) {
                        for (ai, xi) in a.iter_mut().zip(x) {
                            *ai += xi * yj;
                        }
                    }
                }
                for (j, a) in acc.iter().enumerate() {
                    let c = b * NR + j;
                    if c < cols
                        && let Some(dst) = out.get_mut(c * TILE + r * MR..c * TILE + (r + 1) * MR)
                    {
                        dst.copy_from_slice(a);
                    }
                }
            }
        }
        out
    };
    if mu.saturating_mul(wu).saturating_mul(kw) >= PARALLEL_WORK { parallel::map(tiles, tile) } else { (0..tiles).map(tile).collect() }
}

/// Cache-blocked rank update (Golub–Van Loan, matrix multiplication blocking).
/// Row tiles own their outputs throughout all fixed K panels. A panel barrier
/// completes every partial sum before the next panel can update those outputs.
#[cfg(test)]
fn products_cached(l: &Panel, cols: std::ops::Range<usize>, p: usize, mu: usize, wu: usize) -> Vec<Vec<f64>> {
    const KC: usize = 256;
    let mut tiles: Vec<_> = (0..mu.div_ceil(TILE)).map(|t| vec![0.0; wu.min((t + 1) * TILE).min(mu) * TILE]).collect();
    for k0 in (cols.start..cols.end).step_by(KC) {
        let kw = KC.min(cols.end - k0);
        let mut yp = vec![[0.0; NR]; wu.div_ceil(NR) * kw];
        for k in 0..kw {
            for (c, &y) in l.column(k0 + k, p, wu).iter().enumerate() {
                if let Some(slot) = yp.get_mut((c / NR) * kw + k).and_then(|row| row.get_mut(c % NR)) {
                    *slot = y;
                }
            }
        }
        parallel::for_each_mut(&mut tiles, |t, out| {
            let i0 = t * TILE;
            let h = TILE.min(mu - i0);
            let mut xp = vec![[0.0; MR]; TILE / MR * kw];
            for k in 0..kw {
                for (i, &x) in l.column(k0 + k, p + i0, h).iter().enumerate() {
                    if let Some(slot) = xp.get_mut((i / MR) * kw + k).and_then(|row| row.get_mut(i % MR)) {
                        *slot = x;
                    }
                }
            }
            let columns = wu.min(i0 + h);
            for (b, ys) in yp.chunks(kw).enumerate().take(columns.div_ceil(NR)) {
                for (r, xs) in xp.chunks(kw).enumerate() {
                    let mut acc = [[0.0; MR]; NR];
                    for (x, y) in xs.iter().zip(ys) {
                        for (a, yj) in acc.iter_mut().zip(y) {
                            for (ai, xi) in a.iter_mut().zip(x) {
                                *ai += xi * yj;
                            }
                        }
                    }
                    for (j, a) in acc.iter().enumerate() {
                        let c = b * NR + j;
                        if c < columns
                            && let Some(dst) = out.get_mut(c * TILE + r * MR..c * TILE + (r + 1) * MR)
                        {
                            for (v, partial) in dst.iter_mut().zip(a) {
                                *v += partial;
                            }
                        }
                    }
                }
            }
        });
    }
    tiles
}

/// own(rel(i), rel(c)) −= S(i, c) for i ≥ c, with S in the tiles of [`products`].
pub(crate) fn subtract(own: &mut [f64], m: usize, tiles: &[Vec<f64>], mu: usize, wu: usize, rel: impl Fn(usize) -> usize) {
    for (t, tile) in tiles.iter().enumerate() {
        let i0 = t * TILE;
        let h = TILE.min(mu - i0);
        for c in 0..wu.min(i0 + h) {
            let cc = rel(c);
            for i in i0.max(c)..i0 + h {
                let rr = rel(i);
                if let (Some(x), Some(s)) = (own.get_mut(offset(m, cc) + rr - cc), tile.get(c * TILE + i - i0)) {
                    *x -= s;
                }
            }
        }
    }
}

/// Rows c1..m of the panel c0..c1 solved against its factored diagonal block, by tiles of
/// [`TILE`] rows (`tile[c * TILE + i]`): L(r, c) −= L(r, k) L(c, k) for the panel's k < c in
/// order, then ÷ L(c, c).
fn solve_below(l: &Panel, c0: usize, c1: usize) -> Vec<Vec<f64>> {
    let (m, pw) = (l.m, c1 - c0);
    let tiles = (m - c1).div_ceil(TILE);
    let tile = |t: usize| -> Vec<f64> {
        let i0 = c1 + t * TILE;
        let h = TILE.min(m - i0);
        let mut buf = vec![0.0; pw * TILE];
        for c in 0..pw {
            if let Some(dst) = buf.get_mut(c * TILE..c * TILE + h) {
                dst.copy_from_slice(l.column(c0 + c, i0, h));
            }
        }
        for c in 0..pw {
            let (done, rest) = buf.split_at_mut(c * TILE);
            let Some(target) = rest.get_mut(..TILE) else { continue };
            for (k, src) in done.as_chunks::<TILE>().0.iter().enumerate() {
                let lck = l.at(c0 + c, c0 + k);
                for (y, x) in target.iter_mut().zip(src) {
                    *y -= x * lck;
                }
            }
            let d = l.at(c0 + c, c0 + c);
            for y in target.iter_mut() {
                *y /= d;
            }
        }
        buf
    };
    if (m - c1).saturating_mul(pw).saturating_mul(pw) >= PARALLEL_WORK { parallel::map(tiles, tile) } else { (0..tiles).map(tile).collect() }
}

/// The matrix in postorder by lower columns: column j holds A(i, j) for i ≥ j, rows sorted.
pub(crate) struct Lower {
    ptr: Vec<usize>,
    idx: Vec<usize>,
    val: Vec<f64>,
    diag: Vec<f64>,
    /// The largest magnitude in each row (= column) of the whole symmetric matrix.
    pub(crate) scale: Vec<f64>,
}

impl Lower {
    pub(crate) fn new(a: &UpperColumns, inverse: &[usize]) -> Lower {
        let n = a.n;
        let new = |i: usize| inverse.get(i).copied().unwrap_or(0);
        let mut ptr = vec![0usize; n + 1];
        for k in 0..n {
            for &i in a.column(k).0 {
                if let Some(c) = ptr.get_mut(new(i).min(new(k)) + 1) {
                    *c += 1;
                }
            }
        }
        for j in 0..n {
            let before = ptr.get(j).copied().unwrap_or(0);
            if let Some(c) = ptr.get_mut(j + 1) {
                *c += before;
            }
        }
        let mut pairs = vec![(0usize, 0.0f64); ptr.last().copied().unwrap_or(0)];
        let mut next = ptr.clone();
        let mut diag = vec![0.0; n];
        let mut scale = vec![0.0f64; n];
        for k in 0..n {
            let (rows, vals) = a.column(k);
            for (&i, &v) in rows.iter().zip(vals) {
                let (r, c) = (new(i).max(new(k)), new(i).min(new(k)));
                if r == c
                    && let Some(d) = diag.get_mut(c)
                {
                    *d = v;
                }
                for x in [r, c] {
                    if let Some(s) = scale.get_mut(x) {
                        *s = s.max(v.abs());
                    }
                }
                if let Some(slot) = next.get_mut(c) {
                    if let Some(p) = pairs.get_mut(*slot) {
                        *p = (r, v);
                    }
                    *slot += 1;
                }
            }
        }
        for j in 0..n {
            let (s, e) = (ptr.get(j).copied().unwrap_or(0), ptr.get(j + 1).copied().unwrap_or(0));
            if let Some(col) = pairs.get_mut(s..e) {
                col.sort_by_key(|p| p.0);
            }
        }
        let (idx, val) = pairs.into_iter().unzip();
        Lower { ptr, idx, val, diag, scale }
    }

    pub(crate) fn column(&self, j: usize) -> (&[usize], &[f64]) {
        let (a, b) = (self.ptr.get(j).copied().unwrap_or(0), self.ptr.get(j + 1).copied().unwrap_or(0));
        (self.idx.get(a..b).unwrap_or_default(), self.val.get(a..b).unwrap_or_default())
    }
}

/// What a supernode's factorization reads: the matrix, the structure and the updates each
/// supernode receives.
pub(crate) struct Context<'a> {
    pub(crate) a: &'a Lower,
    pub(crate) s: &'a Structure,
    pub(crate) updates: &'a [Vec<Update>],
}

/// What the numeric factorization of one kind keeps per supernode.
pub(crate) trait Supernode: Send + Sync {
    /// The supernode's values: its lower trapezoid by columns.
    fn values(&self) -> &[f64];
    /// The pivots of D (symmetric indefinite factors), `None` for Cholesky.
    fn pivots(&self) -> Option<&[Pivot]>;
    /// Local pivot position to symbolic column within this supernode; identity for Cholesky.
    #[inline]
    fn order(&self) -> Option<&[usize]> {
        None
    }
}

impl Supernode for Vec<f64> {
    #[inline]
    fn values(&self) -> &[f64] {
        self
    }
    #[inline]
    fn pivots(&self) -> Option<&[Pivot]> {
        None
    }
}

/// Shared triangular solve, with Cholesky's diagonal or LDLᵀ's unit L and block diagonal D.
pub(crate) fn solve_supernodes<T: Supernode>(n: usize, perm: &[usize], s: &Structure, blocks: &[T], b: &mut [f64]) -> Result<()> {
    if b.len() != n {
        return Err(Error::Internal("right-hand side length".into()));
    }
    let mut x: Vec<f64> = perm.iter().map(|&j| b.get(j).copied().unwrap_or(0.0)).collect();
    let get = |v: &[f64], i: usize| v.get(i).copied().unwrap_or(0.0);
    let mut scratch = Vec::new();
    // L z = P b, then D⁻¹, by supernodes; the supernode's own unknowns are left in pivot
    // order for the backward pass.
    for (k, blk) in blocks.iter().enumerate() {
        let rows = s.rows(k);
        let (m, f, w) = (rows.len(), s.cols(k).start, s.cols(k).len());
        let t = Panel { m, data: blk.values() };
        let (own, tail) = x.get_mut(f..).and_then(|v| v.split_at_mut_checked(w)).ok_or_else(|| Error::Internal("solve supernode range".into()))?;
        let pivots = blk.pivots();
        let order = blk.order();
        scratch.clear();
        let y = if let Some(order) = order {
            scratch.extend(order.iter().map(|&c| get(own, c)));
            scratch.as_mut_slice()
        } else {
            &mut *own
        };
        for q in 0..w {
            let end = q + if pivots.and_then(|p| p.get(q)) == Some(&Pivot::First) { 2 } else { 1 };
            let mut yq = get(y, q);
            if pivots.is_none() {
                yq /= t.data.get(offset(m, q)).copied().unwrap_or(1.0);
                if let Some(v) = y.get_mut(q) {
                    *v = yq;
                }
            }
            for (v, l) in y.get_mut(end..w).unwrap_or_default().iter_mut().zip(t.column(q, end, w - end)) {
                *v -= l * yq;
            }
            for (&r, l) in rows.get(w..).unwrap_or_default().iter().zip(t.column(q, w, m - w)) {
                if let Some(v) = tail.get_mut(r.wrapping_sub(f + w)) {
                    *v -= l * yq;
                }
            }
        }
        if pivots.is_some() {
            let mut c = 0;
            while c < w {
                if pivots.and_then(|p| p.get(c)) == Some(&Pivot::First) {
                    let (d11, d21, d22) = (t.at(c, c), t.at(c + 1, c), t.at(c + 1, c + 1));
                    let det = d11 * d22 - d21 * d21;
                    let (p, q) = (get(y, c), get(y, c + 1));
                    if let Some(v) = y.get_mut(c) {
                        *v = (d22 * p - d21 * q) / det;
                    }
                    if let Some(v) = y.get_mut(c + 1) {
                        *v = (d11 * q - d21 * p) / det;
                    }
                    c += 2;
                } else {
                    if let Some(v) = y.get_mut(c) {
                        *v /= t.at(c, c);
                    }
                    c += 1;
                }
            }
        }
        if order.is_some() {
            own.copy_from_slice(&scratch);
        }
    }
    // Lᵀ x = z, supernodes in reverse; each scatters its unknowns back from pivot order.
    for (k, blk) in blocks.iter().enumerate().rev() {
        let rows = s.rows(k);
        let (m, f, w) = (rows.len(), s.cols(k).start, s.cols(k).len());
        let t = Panel { m, data: blk.values() };
        let (own, tail) = x.get_mut(f..).and_then(|v| v.split_at_mut_checked(w)).ok_or_else(|| Error::Internal("solve supernode range".into()))?;
        let pivots = blk.pivots();
        let order = blk.order();
        scratch.clear();
        let y = if order.is_some() {
            scratch.extend_from_slice(own);
            scratch.as_mut_slice()
        } else {
            &mut *own
        };
        for c in (0..w).rev() {
            // Reverse pivot groups, keeping ascending columns within each 2 × 2 group.
            let (q, end) = match pivots.and_then(|p| p.get(c)) {
                Some(Pivot::Second) if c > 0 => (c - 1, c + 1),
                Some(Pivot::First) => (c + 1, c + 2),
                _ => (c, c + 1),
            };
            let mut sum = get(y, q);
            for (v, l) in y.get(end..w).unwrap_or_default().iter().zip(t.column(q, end, w - end)) {
                sum -= l * v;
            }
            for (&r, l) in rows.get(w..).unwrap_or_default().iter().zip(t.column(q, w, m - w)) {
                sum -= l * get(tail, r.wrapping_sub(f + w));
            }
            if pivots.is_none() {
                sum /= t.data.get(offset(m, q)).copied().unwrap_or(1.0);
            }
            if let Some(v) = y.get_mut(q) {
                *v = sum;
            }
        }
        if let Some(order) = order {
            for (pos, &c) in order.iter().enumerate() {
                if let Some(v) = own.get_mut(c) {
                    *v = get(&scratch, pos);
                }
            }
        }
    }
    for (&j, v) in perm.iter().zip(&x) {
        if let Some(t) = b.get_mut(j) {
            *t = *v;
        }
    }
    Ok(())
}

fn internal(what: &str) -> FactorError {
    FactorError::Internal(format!("supernodal factorization: {what}"))
}

impl Context<'_> {
    /// Supernode `j` assembled into `own`: its columns of A minus the update L_K Y_Kᵀ of every
    /// descendant K that reaches it (Y = L, or L D), in supernode order; the descendants are
    /// `done[k - base]`.
    pub(crate) fn assemble<T: Supernode>(&self, j: usize, done: &[T], base: usize, own: &mut [f64]) -> std::result::Result<(), FactorError> {
        self.assemble_profiled(j, done, base, own, None)
    }

    fn assemble_profiled<T: Supernode>(
        &self,
        j: usize,
        done: &[T],
        base: usize,
        own: &mut [f64],
        mut timing: Option<&mut crate::timing::Supernode>,
    ) -> std::result::Result<(), FactorError> {
        let rows = self.s.rows(j);
        let cols = self.s.cols(j);
        let (m, w, f) = (rows.len(), cols.len(), cols.start);
        if own.len() != trapezoid(m, w) {
            return Err(internal("supernode size"));
        }
        own.fill(0.0);
        // The columns of A.
        for c in 0..w {
            let (ai, av) = self.a.column(f + c);
            for (r, v) in ai.iter().zip(av) {
                let at = rows.binary_search(r).map_err(|_| internal("an entry of A outside the structure"))?;
                if at < c {
                    return Err(internal("an entry above the diagonal"));
                }
                if let Some(x) = own.get_mut(offset(m, c) + at - c) {
                    *x = *v;
                }
            }
        }
        // The descendants' updates, in supernode order.
        let mut relative = Vec::new();
        for u in self.updates.get(j).map(Vec::as_slice).unwrap_or_default() {
            let done_k = done.get(u.k.wrapping_sub(base)).ok_or_else(|| internal("an update from outside the subtree"))?;
            let lk = done_k.values();
            let (kr, kw) = (self.s.rows(u.k), self.s.cols(u.k).len());
            let below = kr.get(u.p1..).unwrap_or_default();
            let (mu, wu) = (below.len(), u.p2 - u.p1);
            // Positions of K's rows in J's structure (both sorted).
            relative.clear();
            let mut at = 0;
            for r in below {
                while rows.get(at).is_some_and(|x| x < r) {
                    at += 1;
                }
                if rows.get(at) != Some(r) {
                    return Err(internal("an update row outside the structure"));
                }
                relative.push(at);
            }
            #[cfg(test)]
            if path_experiment::enabled(mu, wu, kw) && done_k.pivots().is_none() {
                path_experiment::update(&Panel { m: kr.len(), data: lk }, 0..kw, u.p1, mu, wu, m, 0, &relative, own)?;
                continue;
            }
            let clock = timing.as_ref().and_then(|t| t.clock());
            let tiles = products(&Panel { m: kr.len(), data: lk }, 0..kw, u.p1, mu, wu, done_k.pivots());
            if let Some(timing) = timing.as_mut() {
                timing.product(clock, mu, wu, kw, TILE, NR);
            }
            subtract(own, m, &tiles, mu, wu, |i| relative.get(i).copied().unwrap_or(0));
        }
        Ok(())
    }

    /// Factor supernode `j` into `own`; its descendants' values are `done[k - base]`. Fails with
    /// the first column (postorder) whose pivot is not positive.
    fn supernode(&self, j: usize, done: &[Vec<f64>], base: usize, own: &mut [f64]) -> std::result::Result<(), FactorError> {
        let mut timing = crate::timing::Supernode::new(j, self.s.rows(j).len(), self.s.cols(j).len());
        let clock = timing.clock();
        self.assemble_profiled(j, done, base, own, Some(&mut timing))?;
        timing.assembly_seconds = crate::timing::Supernode::elapsed(clock);
        timing.finish_assembly();
        let rows = self.s.rows(j);
        let cols = self.s.cols(j);
        let (m, w, f) = (rows.len(), cols.len(), cols.start);
        // Dense factorization by panels.
        for c0 in (0..w).step_by(PANEL) {
            let c1 = (c0 + PANEL).min(w);
            let pw = c1 - c0;
            if c0 > 0 {
                // The panel minus its products with the columns before it.
                #[cfg(test)]
                let tiled = path_experiment::enabled(m - c0, pw, c0);
                #[cfg(not(test))]
                let tiled = false;
                #[cfg(test)]
                if tiled {
                    let (done, target) = own.split_at_mut(offset(m, c0));
                    let relative: Vec<_> = (c0..m).collect();
                    path_experiment::update(&Panel { m, data: done }, 0..c0, c0, m - c0, pw, m, c0, &relative, target)?;
                }
                if !tiled {
                    let clock = timing.clock();
                    let tiles = products(&Panel { m, data: own }, 0..c0, c0, m - c0, pw, None);
                    timing.product(clock, m - c0, pw, c0, TILE, NR);
                    subtract(own, m, &tiles, m - c0, pw, |i| c0 + i);
                }
            }
            // The diagonal block, column by column.
            for c in c0..c1 {
                let dst = offset(m, c);
                for k in c0..c {
                    let src = offset(m, k) + c - k;
                    let lck = own.get(src).copied().unwrap_or(0.0);
                    for r in 0..c1 - c {
                        let x = own.get(src + r).copied().unwrap_or(0.0);
                        if let Some(y) = own.get_mut(dst + r) {
                            *y -= x * lck;
                        }
                    }
                }
                let d = own.get(dst).copied().unwrap_or(0.0);
                let diag = self.a.diag.get(f + c).copied().unwrap_or(0.0);
                if !(d > PIVOT_TOLERANCE * diag.abs()) || !d.is_finite() {
                    return Err(FactorError::NotPositiveDefinite { column: f + c });
                }
                let d = d.sqrt();
                if let Some(x) = own.get_mut(dst) {
                    *x = d;
                }
                for x in own.get_mut(dst + 1..dst + c1 - c).unwrap_or_default() {
                    *x /= d;
                }
            }
            // The rows below the panel, by tiles.
            if c1 < m {
                let clock = timing.clock();
                let solved = solve_below(&Panel { m, data: own }, c0, c1);
                timing.below(clock);
                for (t, tile) in solved.iter().enumerate() {
                    let i0 = c1 + t * TILE;
                    let h = TILE.min(m - i0);
                    for c in 0..pw {
                        let start = offset(m, c0 + c) + i0 - (c0 + c);
                        if let (Some(dst), Some(src)) = (own.get_mut(start..start + h), tile.get(c * TILE..c * TILE + h)) {
                            dst.copy_from_slice(src);
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// The symbolic factorization shared by the Cholesky and the symmetric indefinite (LDLᵀ,
/// [`crate::ldlt`]) factors: the postorder, the supernodes and their row structures, the updates
/// each supernode receives, the matrix in that order, and the subtree tasks of the numeric
/// phase.
pub(crate) struct Symbolic {
    pub(crate) n: usize,
    /// `perm[k]`: the original column that is column k of the factor.
    pub(crate) perm: Vec<usize>,
    pub(crate) structure: Structure,
    a: Lower,
    updates: Vec<Vec<Update>>,
    /// Subtrees factored as parallel tasks, and the supernodes above them, in order.
    tasks: Vec<std::ops::Range<usize>>,
    top: Vec<usize>,
    /// Entries of L.
    pub(crate) entries: usize,
}

impl Symbolic {
    /// The structure of the factor of `a`, refused before any numeric memory is taken when the
    /// factor would hold more than `max_entries` entries.
    pub(crate) fn new(a: &UpperColumns, max_entries: usize) -> std::result::Result<Symbolic, FactorError> {
        let _timing = crate::timing::Span::new("symbolic");
        let n = a.n;
        let internal = |what: &str| FactorError::Internal(format!("supernodal structure: {what}"));
        let (parent0, counts0) = symbolic(a);
        let total: usize = counts0.iter().sum();
        if total > max_entries {
            return Err(FactorError::TooLarge { entries: total, limit: max_entries });
        }
        // Postorder: the same tree and counts, renumbered.
        let perm = postorder(&parent0);
        if perm.len() != n {
            return Err(internal("postorder"));
        }
        let mut inverse = vec![0usize; n];
        for (k, &j) in perm.iter().enumerate() {
            if let Some(x) = inverse.get_mut(j) {
                *x = k;
            }
        }
        let parent: Vec<usize> =
            perm.iter().map(|&j| parent0.get(j).copied().filter(|p| *p != NONE).and_then(|p| inverse.get(p).copied()).unwrap_or(NONE)).collect();
        let counts: Vec<usize> = perm.iter().map(|&j| counts0.get(j).copied().unwrap_or(1)).collect();
        let a = Lower::new(a, &inverse);

        // Fundamental supernodes: column j joins column j − 1's supernode when it is j − 1's
        // parent, j − 1 is its only child, and their structures nest.
        let mut children = vec![0usize; n];
        for &p in &parent {
            if let Some(c) = children.get_mut(p) {
                *c += 1;
            }
        }
        let mut first = vec![0usize];
        for j in 1..n {
            let joins = parent.get(j - 1) == Some(&j) && children.get(j) == Some(&1) && counts.get(j - 1).copied() == counts.get(j).map(|c| c + 1);
            if !joins {
                first.push(j);
            }
        }
        if n > 0 {
            first.push(n);
        }
        let ns = first.len().saturating_sub(1);
        let mut snode = vec![0usize; n];
        for s in 0..ns {
            for j in first.get(s).copied().unwrap_or(0)..first.get(s + 1).copied().unwrap_or(0) {
                if let Some(x) = snode.get_mut(j) {
                    *x = s;
                }
            }
        }
        // The supernodal tree and the row structures: a supernode's rows are its columns, the
        // rows of A below them and its children's rows past their own columns.
        let sparent: Vec<usize> = (0..ns)
            .map(|s| {
                let last = first.get(s + 1).copied().unwrap_or(1).saturating_sub(1);
                parent.get(last).copied().filter(|p| *p != NONE).and_then(|p| snode.get(p).copied()).unwrap_or(NONE)
            })
            .collect();
        let mut kids: Vec<Vec<usize>> = vec![Vec::new(); ns];
        for (s, &p) in sparent.iter().enumerate() {
            if let Some(k) = kids.get_mut(p) {
                k.push(s);
            }
        }
        let mut rptr = vec![0usize];
        let mut rows: Vec<usize> = Vec::new();
        let mut mark = vec![NONE; n];
        let mut set = Vec::new();
        for s in 0..ns {
            let (f, l) = (first.get(s).copied().unwrap_or(0), first.get(s + 1).copied().unwrap_or(0));
            set.clear();
            let mut add = |r: usize, set: &mut Vec<usize>| {
                if let Some(x) = mark.get_mut(r)
                    && *x != s
                {
                    *x = s;
                    set.push(r);
                }
            };
            for j in f..l {
                add(j, &mut set);
                for &r in a.column(j).0 {
                    add(r, &mut set);
                }
            }
            for &c in kids.get(s).map(Vec::as_slice).unwrap_or_default() {
                let (cs, ce) = (rptr.get(c).copied().unwrap_or(0), rptr.get(c + 1).copied().unwrap_or(0));
                let cw = first.get(c + 1).copied().unwrap_or(0) - first.get(c).copied().unwrap_or(0);
                for &r in rows.get(cs + cw..ce).unwrap_or_default() {
                    add(r, &mut set);
                }
            }
            set.sort_unstable();
            if set.len() != counts.get(f).copied().unwrap_or(0) || set.first() != Some(&f) {
                return Err(internal("a supernode's rows differ from its column count"));
            }
            rows.extend_from_slice(&set);
            rptr.push(rows.len());
        }
        let structure = Structure { first, rptr, rows };

        // The updates each supernode receives, from its descendants in increasing order.
        let mut updates: Vec<Vec<Update>> = vec![Vec::new(); ns];
        for k in 0..ns {
            let kr = structure.rows(k);
            let mut p = structure.cols(k).len();
            while let Some(target) = kr.get(p).and_then(|r| snode.get(*r)).copied() {
                let end = structure.cols(target).end;
                let mut q = p;
                while kr.get(q).is_some_and(|r| *r < end) {
                    q += 1;
                }
                if let Some(list) = updates.get_mut(target) {
                    list.push(Update { k, p1: p, p2: q });
                }
                p = q;
            }
        }

        // Tasks: subtrees (contiguous ranges of supernodes in postorder) of at most a share of the
        // work each; the supernodes above them are factored afterwards, in order.
        let mut subtree: Vec<usize> = (0..ns)
            .map(|s| {
                let (m, w) = (structure.rows(s).len(), structure.cols(s).len());
                (0..w).map(|c| (m - c).saturating_mul(m - c)).fold(0usize, usize::saturating_add)
            })
            .collect();
        let mut start: Vec<usize> = (0..ns).collect();
        for s in 0..ns {
            if let Some(&p) = sparent.get(s)
                && p != NONE
            {
                let (ws, ss) = (subtree.get(s).copied().unwrap_or(0), start.get(s).copied().unwrap_or(s));
                if let Some(x) = subtree.get_mut(p) {
                    *x = x.saturating_add(ws);
                }
                if let Some(x) = start.get_mut(p) {
                    *x = (*x).min(ss);
                }
            }
        }
        let all = (0..ns).filter(|s| sparent.get(*s) == Some(&NONE)).filter_map(|s| subtree.get(s).copied()).fold(0, usize::saturating_add);
        let share = (all / SUBTREE_SHARE).max(PARALLEL_WORK);
        let fits = |s: usize| subtree.get(s).is_some_and(|w| *w <= share);
        let mut tasks: Vec<std::ops::Range<usize>> = Vec::new();
        let mut top: Vec<usize> = Vec::new();
        let mut s = 0;
        while s < ns {
            // The largest fitting subtree whose first supernode is s.
            let mut root = NONE;
            let mut r = s;
            while r != NONE && start.get(r) == Some(&s) && fits(r) {
                root = r;
                r = sparent.get(r).copied().unwrap_or(NONE);
            }
            if root == NONE {
                top.push(s);
                s += 1;
            } else {
                tasks.push(s..root + 1);
                s = root + 1;
            }
        }

        Ok(Symbolic { n, perm, structure, a, updates, tasks, top, entries: total })
    }

    /// Estimate from allocated symbolic capacities and the largest live buffers of each
    /// scheduled subtree. Subtrees may overlap in time; top supernodes run afterwards.
    pub(crate) fn memory_estimate(&self, upper: &UpperColumns) -> FactorMemory {
        let sum = |xs: &[usize]| xs.iter().copied().fold(0usize, usize::saturating_add);
        let structure = sum(&[buffer_bytes(&self.structure.first), buffer_bytes(&self.structure.rptr), buffer_bytes(&self.structure.rows)]);
        let indices_bytes = sum(&[
            structure.saturating_mul(2),
            buffer_bytes(&self.perm).saturating_mul(2),
            buffer_bytes(&self.updates),
            self.updates.iter().map(buffer_bytes).fold(0usize, usize::saturating_add),
            buffer_bytes(&self.tasks),
            buffer_bytes(&self.top),
            // LDLT pivot order and pivot tags, plus block/vector headers and task bookkeeping.
            self.n.saturating_mul(std::mem::size_of::<usize>() + std::mem::size_of::<Pivot>()),
            self.supernodes().saturating_mul(256),
        ]);
        let matrix_bytes = sum(&[
            buffer_bytes(&upper.ptr),
            buffer_bytes(&upper.idx),
            buffer_bytes(&upper.val),
            buffer_bytes(&self.a.ptr),
            buffer_bytes(&self.a.idx),
            buffer_bytes(&self.a.val),
            buffer_bytes(&self.a.diag),
            buffer_bytes(&self.a.scale),
        ]);
        let node_work = |j: usize| {
            let (m, w) = (self.structure.rows(j).len(), self.structure.cols(j).len());
            let mut work = products_bytes(m, PANEL.min(w), w);
            for u in self.updates.get(j).map(Vec::as_slice).unwrap_or_default() {
                let mu = self.structure.rows(u.k).len().saturating_sub(u.p1);
                let wu = u.p2.saturating_sub(u.p1);
                let kw = self.structure.cols(u.k).len();
                work = work.max(products_bytes(mu, wu, kw));
            }
            // Relative-row mapping and LDLT's three reusable scratch vectors, including growth.
            work.saturating_add(m.saturating_mul(64))
        };
        let parallel_work = self.tasks.iter().map(|r| r.clone().map(node_work).max().unwrap_or(0)).fold(0usize, usize::saturating_add);
        let top_work = self.top.iter().copied().map(node_work).max().unwrap_or(0);
        // The symbolic phase has column/tree arrays and temporary (index,value) pairs. Count
        // them conservatively alongside retained storage instead of assuming perfect reuse.
        let symbolic_work = self.n.saturating_mul(192).saturating_add(upper.val.len().saturating_mul(16));
        FactorMemory {
            entries: self.entries,
            values_bytes: self.entries.saturating_mul(8),
            indices_bytes,
            matrix_bytes,
            workspace_bytes: parallel_work.max(top_work).max(symbolic_work),
        }
    }

    /// Supernodes.
    pub(crate) fn supernodes(&self) -> usize {
        self.structure.first.len().saturating_sub(1)
    }

    /// The original column of factor column `column`.
    pub(crate) fn original(&self, column: usize) -> usize {
        self.perm.get(column).copied().unwrap_or(column)
    }

    /// Zeroed memory for every supernode's trapezoid, taken without aborting.
    pub(crate) fn allocate(&self) -> std::result::Result<Vec<Vec<f64>>, FactorError> {
        let ns = self.supernodes();
        let mut values: Vec<Vec<f64>> = Vec::new();
        values.try_reserve_exact(ns).map_err(|_| FactorError::OutOfMemory { entries: self.entries })?;
        for s in 0..ns {
            let mut v = Vec::new();
            v.try_reserve_exact(self.structure.size(s)).map_err(|_| FactorError::OutOfMemory { entries: self.entries })?;
            v.resize(self.structure.size(s), 0.0);
            values.push(v);
        }
        Ok(values)
    }

    /// Run `supernode(context, j, done, base, own)` over every supernode, children before
    /// parents: the subtree tasks in parallel, each through its supernodes in order, then the
    /// supernodes above them in order. A numeric failure is the first failing column in
    /// postorder (the one a pass in order would meet first, whichever thread meets it), named by
    /// its original column.
    pub(crate) fn factor<T: Send + Sync>(
        &self,
        values: &mut [T],
        supernode: impl Fn(&Context, usize, &[T], usize, &mut T) -> std::result::Result<(), FactorError> + Sync,
    ) -> std::result::Result<(), FactorError> {
        let context = Context { a: &self.a, s: &self.structure, updates: &self.updates };
        // The tasks run in parallel, each through its supernodes in order.
        let mut pieces: Vec<(std::ops::Range<usize>, &mut [T], Option<FactorError>)> = Vec::with_capacity(self.tasks.len());
        let mut rest: &mut [T] = values;
        let mut at = 0;
        for range in &self.tasks {
            let (_, tail) = std::mem::take(&mut rest).split_at_mut(range.start - at);
            let (piece, tail) = tail.split_at_mut(range.len());
            pieces.push((range.clone(), piece, None));
            rest = tail;
            at = range.end;
        }
        let broken = std::sync::Mutex::new(None::<FactorError>);
        // The pool's threads do not report: they ask the caller's progress whether to stop.
        let progress = fieldcraft_model::progress::current();
        let cancelled = || progress.as_ref().is_some_and(|p| p.is_cancelled());
        parallel::for_each_mut(&mut pieces, |_, (range, piece, failed)| {
            for j in range.clone() {
                if cancelled() {
                    if let Ok(mut b) = broken.lock() {
                        b.get_or_insert(FactorError::Cancelled);
                    }
                    break;
                }
                let (done, own) = piece.split_at_mut(j - range.start);
                let Some(own) = own.first_mut() else { break };
                match supernode(&context, j, done, range.start, own) {
                    Ok(()) => {}
                    Err(e) if e.column().is_some() => {
                        *failed = Some(e);
                        break;
                    }
                    Err(e) => {
                        if let Ok(mut b) = broken.lock() {
                            b.get_or_insert(e);
                        }
                        break;
                    }
                }
            }
        });
        if let Some(e) = broken.into_inner().unwrap_or_else(std::sync::PoisonError::into_inner) {
            return Err(e);
        }
        // The first failing column in postorder is the one a pass in order would meet first:
        // the supernodes above the tasks run up to it.
        let earliest = pieces.iter_mut().filter_map(|p| p.2.take()).min_by_key(|e| e.column());
        drop(pieces);
        let named = |e: FactorError| match e.column() {
            Some(c) => e.at(self.original(c)),
            None => e,
        };
        for &j in &self.top {
            if cancelled() {
                return Err(FactorError::Cancelled);
            }
            if earliest.as_ref().and_then(FactorError::column).is_some_and(|c| self.structure.cols(j).start > c) {
                break;
            }
            let (done, own) = values.split_at_mut(j);
            let Some(own) = own.first_mut() else { break };
            supernode(&context, j, done, 0, own).map_err(named)?;
        }
        match earliest {
            Some(e) => Err(named(e)),
            None => Ok(()),
        }
    }
}

/// A Cholesky factor: L in supernodes, for the matrix with its columns in postorder.
#[derive(Debug, Clone)]
pub struct Factor {
    pub n: usize,
    /// `perm[k]`: the original column that is column k of the factor.
    perm: Vec<usize>,
    structure: Structure,
    /// The values of each supernode (a lower trapezoid by columns).
    values: Vec<Vec<f64>>,
    entries: usize,
}

impl Factor {
    pub fn new(a: &UpperColumns, max_entries: usize) -> std::result::Result<Factor, FactorError> {
        Factor::from_symbolic(&Symbolic::new(a, max_entries)?)
    }

    /// Estimate factor storage without allocating the numeric factor. Symbolic storage is
    /// still allocated; `max_entries` remains an independent cap on the symbolic count.
    pub fn estimate_memory(a: &UpperColumns, max_entries: usize) -> std::result::Result<FactorMemory, FactorError> {
        Ok(Symbolic::new(a, max_entries)?.memory_estimate(a))
    }

    /// Factor within a factor-only byte budget, checked before allocating numeric values.
    /// Callers must deduct their other live allocations from the process budget first.
    pub fn with_memory_limit(a: &UpperColumns, max_entries: usize, max_bytes: usize) -> std::result::Result<(Factor, FactorMemory), FactorError> {
        let symbolic = Symbolic::new(a, max_entries)?;
        let estimate = symbolic.memory_estimate(a);
        estimate.check(max_bytes)?;
        Ok((Factor::from_symbolic(&symbolic)?, estimate))
    }

    pub(crate) fn from_symbolic(symbolic: &Symbolic) -> std::result::Result<Factor, FactorError> {
        let allocation_timing = crate::timing::Span::new("factor_allocate");
        let mut values = symbolic.allocate()?;
        drop(allocation_timing);
        if std::env::var("FIELDCRAFT_SUPERNODE_TIMING").is_ok_and(|v| v == "1") {
            use std::io::Write;
            let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_SUPERNODE_TOP\t{:?}", symbolic.top);
        }
        let numeric_timing = crate::timing::Span::new("numeric_factor");
        symbolic.factor(&mut values, |context, j, done, base, own| context.supernode(j, done, base, own))?;
        drop(numeric_timing);
        Ok(Factor { n: symbolic.n, perm: symbolic.perm.clone(), structure: symbolic.structure.clone(), values, entries: symbolic.entries })
    }

    /// Solve independent right-hand sides with factor entries reused across up to eight
    /// columns. Each lane preserves solve_supernodes' column and row subtraction order.
    /// This is the multiple-RHS triangular solve of Golub–Van Loan §3.1; vectorization
    /// is across independent RHS lanes, never a reassociation of a dot product.
    pub fn solve_block(&self, rhs: &mut [Vec<f64>]) -> Result<()> {
        let _timing = crate::timing::Span::new("backsolve_block");
        self.solve_lanes::<8>(rhs)
    }

    /// Independent pairs share one factor stream, with one writer for each RHS. This
    /// retains column parallelism as well as multiple-RHS reuse; no dot product is split.
    #[cfg(test)]
    fn solve_pairs(&self, rhs: &mut [Vec<f64>]) -> Result<()> {
        let _timing = crate::timing::Span::new("backsolve_pairs");
        if rhs.iter().any(|b| b.len() != self.n) {
            return Err(Error::Internal("paired right-hand side length".into()));
        }
        let mut groups: Vec<_> = rhs.chunks_mut(2).map(|group| (group, Ok(()))).collect();
        parallel::for_each_mut(&mut groups, |_, (group, result)| {
            *result = if group.len() == 1 {
                match group.first_mut() {
                    Some(rhs) => self.solve(rhs),
                    None => Ok(()),
                }
            } else {
                self.solve_lanes::<2>(group)
            };
        });
        for (_, result) in groups {
            result?;
        }
        Ok(())
    }

    fn solve_lanes<const LANES: usize>(&self, rhs: &mut [Vec<f64>]) -> Result<()> {
        let _timing = crate::timing::Span::new("backsolve_lane_group_inclusive");
        if LANES == 0 {
            return Err(Error::Internal("zero right-hand side lanes".into()));
        }
        if rhs.iter().any(|b| b.len() != self.n) {
            return Err(Error::Internal("block right-hand side length".into()));
        }
        for group in rhs.chunks_mut(LANES) {
            let mut x = vec![[0.0; LANES]; self.n];
            for (&p, row) in self.perm.iter().zip(&mut x) {
                for (v, b) in row.iter_mut().zip(group.iter()) {
                    *v = b.get(p).copied().unwrap_or(0.0);
                }
            }
            for (k, values) in self.values.iter().enumerate() {
                let rows = self.structure.rows(k);
                let (m, f, w) = (rows.len(), self.structure.cols(k).start, self.structure.cols(k).len());
                let t = Panel { m, data: values };
                for q in 0..w {
                    let y = x.get_mut(f + q).ok_or_else(|| Error::Internal("block solve column".into()))?;
                    let d = t.at(q, q);
                    for v in y.iter_mut() {
                        *v /= d;
                    }
                    let yq = *y;
                    for (row, &l) in x.get_mut(f + q + 1..f + w).unwrap_or_default().iter_mut().zip(t.column(q, q + 1, w - q - 1)) {
                        for (v, y) in row.iter_mut().zip(yq) {
                            *v -= l * y;
                        }
                    }
                    for (&r, &l) in rows.get(w..).unwrap_or_default().iter().zip(t.column(q, w, m - w)) {
                        if let Some(row) = x.get_mut(r) {
                            for (v, y) in row.iter_mut().zip(yq) {
                                *v -= l * y;
                            }
                        }
                    }
                }
            }
            for (k, values) in self.values.iter().enumerate().rev() {
                let rows = self.structure.rows(k);
                let (m, f, w) = (rows.len(), self.structure.cols(k).start, self.structure.cols(k).len());
                let t = Panel { m, data: values };
                for q in (0..w).rev() {
                    let mut sum = x.get(f + q).copied().ok_or_else(|| Error::Internal("block solve column".into()))?;
                    for (row, &l) in x.get(f + q + 1..f + w).unwrap_or_default().iter().zip(t.column(q, q + 1, w - q - 1)) {
                        for (v, y) in sum.iter_mut().zip(row) {
                            *v -= l * y;
                        }
                    }
                    for (&r, &l) in rows.get(w..).unwrap_or_default().iter().zip(t.column(q, w, m - w)) {
                        if let Some(row) = x.get(r) {
                            for (v, y) in sum.iter_mut().zip(row) {
                                *v -= l * y;
                            }
                        }
                    }
                    let d = t.at(q, q);
                    for v in &mut sum {
                        *v /= d;
                    }
                    if let Some(row) = x.get_mut(f + q) {
                        *row = sum;
                    }
                }
            }
            for (&p, row) in self.perm.iter().zip(&x) {
                for (b, v) in group.iter_mut().zip(row) {
                    if let Some(y) = b.get_mut(p) {
                        *y = *v;
                    }
                }
            }
            crate::timing::triangular_work(group.len(), LANES, self.n, self.entries);
        }
        Ok(())
    }

    /// Entries in L.
    pub fn entries(&self) -> usize {
        self.entries
    }

    /// Solve A x = b in place.
    pub fn solve(&self, b: &mut [f64]) -> Result<()> {
        let _timing = crate::timing::Span::new("backsolve_inclusive");
        let result = solve_supernodes(self.n, &self.perm, &self.structure, &self.values, b);
        if result.is_ok() {
            crate::timing::triangular_work(1, 1, self.n, self.entries);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Csr;

    #[test]
    #[ignore = "explicit bounded cache-panel factorization experiment"]
    fn cache_panel_factor_experiment() {
        let large = std::env::var("FC_FACTOR_SIZE").is_ok_and(|v| v == "large");
        let (child, separator) = if large { (1024, 4911) } else { (384, 768) };
        let n = 2 * child + separator;
        let mut a = UpperColumns { n, ptr: vec![0], idx: vec![], val: vec![] };
        for c in 0..n {
            for r in 0..=c {
                if c >= 2 * child || r / child == c / child {
                    a.idx.push(r);
                    a.val.push(if r == c { n as f64 + 1.0 } else { ((r * 17 + c * 31) % 29) as f64 / 29.0 - 0.5 });
                }
            }
            a.ptr.push(a.idx.len());
        }
        let threads = std::env::var("FC_FACTOR_THREADS").unwrap().parse().unwrap();
        let start = std::time::Instant::now();
        let factor = parallel::with_threads(threads, || Factor::new(&a, usize::MAX)).unwrap();
        let factor_seconds = start.elapsed().as_secs_f64();
        let mut x = vec![1.0; n];
        factor.solve(&mut x).unwrap();
        let solve_seconds = start.elapsed().as_secs_f64();
        let mut ax = vec![0.0; n];
        for c in 0..n {
            let (rows, values) = a.column(c);
            for (&r, &v) in rows.iter().zip(values) {
                ax[r] += v * x[c];
                if r != c {
                    ax[c] += v * x[r];
                }
            }
        }
        let residual = ax.iter().map(|v| (v - 1.0).powi(2)).sum::<f64>().sqrt() / (n as f64).sqrt();
        assert!(residual < 1e-11, "{residual}");
        let hash = factor.values.iter().flatten().fold(0xcbf29ce484222325_u64, |h, v| (h ^ v.to_bits()).wrapping_mul(0x100000001b3));
        let bytes: Vec<_> = x.iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(std::env::var("FC_FACTOR_OUTPUT").unwrap(), bytes).unwrap();
        let largest = (0..factor.values.len()).map(|j| (factor.structure.cols(j).len(), factor.structure.rows(j).len())).max();
        eprintln!(
            "CACHE_FACTOR n={n} child={child} separator={separator} largest={largest:?} factor_seconds={factor_seconds:.9} factor_solve_seconds={solve_seconds:.9} residual={residual:.17e} factor_hash={hash:016x}"
        );
    }

    #[test]
    fn cached_update_fixed_k_panels_and_row_tail() {
        let (m, kw, mu, wu) = (800, 517, 269, 163);
        let values: Vec<f64> = (0..trapezoid(m, kw)).map(|i| ((i * 37 % 101) as f64 - 50.0) / 101.0).collect();
        let panel = Panel { m, data: &values };
        let mut bits = None;
        for threads in [1, 2, 7] {
            let got = parallel::with_threads(threads, || products_cached(&panel, 0..kw, 530, mu, wu));
            let mut max_error = 0.0_f64;
            for (t, tile) in got.iter().enumerate() {
                for c in 0..wu.min(t * TILE + TILE.min(mu - t * TILE)) {
                    for i in 0..TILE.min(mu - t * TILE) {
                        let reference: f64 = (0..kw).map(|k| panel.at(530 + t * TILE + i, k) * panel.at(530 + c, k)).sum();
                        max_error = max_error.max((tile[c * TILE + i] - reference).abs());
                    }
                }
            }
            assert!(max_error < 1e-11, "{max_error}");
            let got: Vec<_> = got.iter().flatten().map(|v| v.to_bits()).collect();
            if let Some(expected) = &bits {
                assert_eq!(&got, expected);
            } else {
                bits = Some(got);
            }
        }
    }

    #[test]
    fn multiple_rhs_preserves_scalar_solve_bits() {
        let n = 65;
        let a: Vec<Vec<f64>> =
            (0..n).map(|i| (0..n).map(|j| if i == j { 100.0 } else { ((i.min(j) * 7 + i.max(j) * 13) % 23) as f64 / 23.0 }).collect()).collect();
        let factor = Factor::new(&UpperColumns::from_symmetric(&dense_to_csr(&a)), usize::MAX).unwrap();
        for width in [0, 1, 6, 9] {
            let rhs: Vec<Vec<f64>> = (0..width).map(|k| (0..n).map(|j| ((j * 7 + k * 11) % 31) as f64 - 15.0).collect()).collect();
            let mut reference = rhs.clone();
            for b in &mut reference {
                factor.solve(b).unwrap();
            }
            for threads in [1, 2, 7] {
                let mut got = rhs.clone();
                parallel::with_threads(threads, || factor.solve_block(&mut got)).unwrap();
                assert!(got.iter().flatten().zip(reference.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
                let mut paired = rhs.clone();
                parallel::with_threads(threads, || factor.solve_pairs(&mut paired)).unwrap();
                assert!(paired.iter().flatten().zip(reference.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
            }
        }
        assert!(factor.solve_block(&mut [vec![1.0; n - 1]]).is_err());
        let mut bad = vec![vec![1.0; n], vec![2.0; n - 1]];
        let original = bad.clone();
        assert!(factor.solve_pairs(&mut bad).is_err());
        assert_eq!(bad, original);
    }

    #[test]
    fn memory_budget_refuses_before_numeric_and_preserves_values() {
        let upper = UpperColumns::from_symmetric(&dense_to_csr(&[vec![4.0, 1.0], vec![1.0, 3.0]]));
        let estimate = Factor::estimate_memory(&upper, 10).unwrap();
        assert_eq!(estimate.entries, 3);
        assert_eq!(estimate.values_bytes, 24);
        assert!(estimate.indices_bytes > 0 && estimate.matrix_bytes > 0 && estimate.workspace_bytes > 0);
        assert!(matches!(Factor::with_memory_limit(&upper, 10, estimate.total_bytes() - 1), Err(FactorError::MemoryLimit { entries: 3, .. })));
        let (bounded, observed) = Factor::with_memory_limit(&upper, 10, estimate.total_bytes()).unwrap();
        let original = Factor::new(&upper, 10).unwrap();
        assert_eq!(observed, estimate);
        assert_eq!(bounded.values, original.values);
        assert_eq!(bounded.perm, original.perm);
        assert_eq!(bounded.values.iter().map(|v| v.len() * 8).sum::<usize>(), estimate.values_bytes);
        let mut x = vec![1.0, 2.0];
        let mut y = x.clone();
        bounded.solve(&mut x).unwrap();
        original.solve(&mut y).unwrap();
        assert_eq!(x, y);
        // With room for symbolic work but no numeric storage, memory wins over the bad pivot.
        let bad = UpperColumns::from_symmetric(&dense_to_csr(&[vec![-1.0]]));
        assert!(matches!(Factor::with_memory_limit(&bad, 10, 0), Err(FactorError::MemoryLimit { .. })));
        assert!(matches!(Factor::with_memory_limit(&upper, 2, usize::MAX), Err(FactorError::TooLarge { .. })));
    }

    #[test]
    fn memory_estimates_fail_closed_on_overflow() {
        let m = FactorMemory { entries: 1, values_bytes: usize::MAX, indices_bytes: 1, matrix_bytes: 0, workspace_bytes: 0 };
        assert_eq!(m.total_bytes(), usize::MAX);
        assert!(m.check(usize::MAX).is_err());
        assert_eq!(products_bytes(usize::MAX, usize::MAX, usize::MAX), usize::MAX);
    }

    #[test]
    fn symmetric_memory_budget_keeps_auto_and_indefinite_results() {
        use crate::ldlt::{Kind, Symmetric};
        let a = UpperColumns::from_symmetric(&dense_to_csr(&[vec![0.0, 1.0], vec![1.0, 0.0]]));
        for kind in [Kind::Auto, Kind::Indefinite] {
            assert!(matches!(Symmetric::with_memory_limit(&a, 10, kind, 0), Err(FactorError::MemoryLimit { .. })));
            let (bounded, estimate) = Symmetric::with_memory_limit(&a, 10, kind, 1_000_000).unwrap();
            let original = Symmetric::new(&a, 10, kind).unwrap();
            assert_eq!(bounded.entries(), estimate.entries);
            let mut x = vec![2.0, 3.0];
            let mut y = x.clone();
            bounded.solve(&mut x).unwrap();
            original.solve(&mut y).unwrap();
            assert_eq!(x, y);
            assert_eq!(x, vec![3.0, 2.0]);
        }
    }

    /// [`products`] one output column at a time, its sum over k in order: the reference the
    /// micro-tiles reorganise.
    fn products_by_columns(l: &Panel, cols: std::ops::Range<usize>, p: usize, mu: usize, wu: usize, pivots: Option<&[Pivot]>) -> Vec<Vec<f64>> {
        let (k0, kw) = (cols.start, cols.len());
        let mut yt = vec![0.0; wu * kw];
        for k in 0..kw {
            for (c, y) in l.column(k0 + k, p, wu).iter().enumerate() {
                yt[c * kw + k] = *y;
            }
        }
        if let Some(pivots) = pivots {
            for row in yt.chunks_mut(kw.max(1)) {
                times_d(row, l, k0, pivots);
            }
        }
        (0..mu.div_ceil(TILE))
            .map(|t| {
                let i0 = t * TILE;
                let h = TILE.min(mu - i0);
                let mut x = vec![0.0; kw * TILE];
                for k in 0..kw {
                    x[k * TILE..k * TILE + h].copy_from_slice(l.column(k0 + k, p + i0, h));
                }
                let cols = wu.min(i0 + h);
                let mut out = vec![0.0; cols * TILE];
                for c in 0..cols {
                    let mut acc = [0.0f64; TILE];
                    for (xs, b) in x.as_chunks::<TILE>().0.iter().zip(&yt[c * kw..(c + 1) * kw]) {
                        for (a, xi) in acc.iter_mut().zip(xs) {
                            *a += xi * b;
                        }
                    }
                    out[c * TILE..(c + 1) * TILE].copy_from_slice(&acc);
                }
                out
            })
            .collect()
    }

    /// The micro-tiled products equal the column-at-a-time ones bit for bit: panels narrower and
    /// wider than a micro-tile, row counts across tile edges, with and without a block diagonal
    /// of 1 × 1 and 2 × 2 pivots.
    #[test]
    fn micro_tiled_products_are_the_column_sums_bit_for_bit() {
        let mut state = 11u64;
        let mut random = move || {
            state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for (m, w, k0, kw, p, mu, wu) in
            [(5, 3, 0, 3, 3, 2, 1), (40, 12, 0, 7, 7, 33, 5), (90, 20, 2, 9, 11, 79, 13), (200, 70, 0, 64, 64, 136, 6), (300, 90, 5, 30, 40, 260, 50)]
        {
            let data: Vec<f64> = (0..trapezoid(m, w)).map(|_| random()).collect();
            let l = Panel { m, data: &data };
            let pivots: Vec<Pivot> = (0..kw)
                .map(|k| {
                    [Pivot::One, Pivot::First, Pivot::Second][if k % 5 == 1 && k + 1 < kw {
                        1
                    } else if k % 5 == 2 {
                        2
                    } else {
                        0
                    }]
                })
                .collect();
            for piv in [None, Some(&pivots[..])] {
                let got = products(&l, k0..k0 + kw, p, mu, wu, piv);
                let expected = products_by_columns(&l, k0..k0 + kw, p, mu, wu, piv);
                let bits = |t: &Vec<Vec<f64>>| t.iter().map(|v| v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>()).collect::<Vec<_>>();
                assert_eq!(bits(&got), bits(&expected), "m {m}, cols {k0}..{}, p {p}, mu {mu}, wu {wu}, pivots {}", k0 + kw, piv.is_some());
            }
        }
    }

    fn dense_to_csr(a: &[Vec<f64>]) -> Csr {
        let n = a.len();
        let mut ptr = vec![0];
        let mut idx = Vec::new();
        let mut val = Vec::new();
        for row in a {
            for (j, v) in row.iter().enumerate() {
                if *v != 0.0 {
                    idx.push(j);
                    val.push(*v);
                }
            }
            ptr.push(idx.len());
        }
        Csr { n, ptr, idx, val }
    }

    fn check_solve(csr: &Csr, f: &Factor, tol: f64) {
        let n = csr.n;
        let x: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let mut b = vec![0.0; n];
        csr.mul(&x, &mut b);
        f.solve(&mut b).unwrap();
        let err = b.iter().zip(&x).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
        assert!(err < tol, "{err}");
    }

    #[test]
    fn factor_and_solve_a_sparse_spd_matrix() {
        // A 2D Laplacian-like matrix on a 7 × 5 grid plus a diagonal shift: sparse with fill.
        let (w, h) = (7, 5);
        let n = w * h;
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..w {
            for j in 0..h {
                let k = i * h + j;
                a[k][k] = 4.5;
                if i + 1 < w {
                    a[k][k + h] = -1.0;
                    a[k + h][k] = -1.0;
                }
                if j + 1 < h {
                    a[k][k + 1] = -1.0;
                    a[k + 1][k] = -1.0;
                }
            }
        }
        let csr = dense_to_csr(&a);
        let f = Factor::new(&UpperColumns::from_symmetric(&csr), usize::MAX).unwrap();
        check_solve(&csr, &f, 1e-13);
        // The cap is checked before the numeric factorization.
        assert!(matches!(Factor::new(&UpperColumns::from_symmetric(&csr), 10), Err(FactorError::TooLarge { .. })));
        // A singular matrix (a pure Laplacian has the constant vector in its null space).
        let mut sing = a.clone();
        for (k, row) in sing.iter_mut().enumerate() {
            let off: f64 = row.iter().enumerate().filter(|(j, _)| *j != k).map(|(_, v)| -v).sum();
            row[k] = off;
        }
        let e = Factor::new(&UpperColumns::from_symmetric(&dense_to_csr(&sing)), usize::MAX).unwrap_err();
        assert_eq!(e, FactorError::NotPositiveDefinite { column: n - 1 });
        // A cancelled command stops between supernodes, in Cholesky and LDLᵀ alike.
        let p = fieldcraft_model::progress::Progress::new();
        p.cancel();
        fieldcraft_model::progress::scoped(Some(p), || {
            let upper = UpperColumns::from_symmetric(&csr);
            assert_eq!(Factor::new(&upper, usize::MAX).unwrap_err(), FactorError::Cancelled);
            assert_eq!(crate::ldlt::LdltFactor::new(&upper, usize::MAX).unwrap_err(), FactorError::Cancelled);
        });
        // Empty and 1 × 1 matrices.
        let empty = Csr { n: 0, ptr: vec![0], idx: vec![], val: vec![] };
        let f = Factor::new(&UpperColumns::from_symmetric(&empty), 10).unwrap();
        assert_eq!(f.entries(), 0);
        f.solve(&mut []).unwrap();
        let f = Factor::new(&UpperColumns::from_symmetric(&dense_to_csr(&[vec![4.0]])), 10).unwrap();
        let mut b = [2.0];
        f.solve(&mut b).unwrap();
        assert_eq!(b, [0.5]);
    }

    /// Two 3D grids (3 coupled unknowns per point, in grid order) joined by a dense separator
    /// ordered last, as nested dissection orders them, and a separate 5 × 5 block: subtrees run
    /// as parallel tasks, the separator is a supernode of several panels whose kernels run in
    /// parallel tiles. The factor solves the system and is the same bit for bit on any number of
    /// threads.
    #[test]
    fn supernodes_panels_and_threads() {
        let (nx, ny, nz) = (8, 8, 6);
        let point = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
        let per = 3 * nx * ny * nz;
        let sep = 400;
        let n = 2 * per + sep + 5;
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let pair = |rows: &mut Vec<Vec<(usize, f64)>>, a: usize, b: usize, v: f64| {
            rows[a].push((b, v));
            rows[b].push((a, v));
        };
        for base in [0, per] {
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let p = base + 3 * point(i, j, k);
                        for a in 0..3 {
                            rows[p + a].push((p + a, 7.0));
                            for b in a + 1..3 {
                                pair(&mut rows, p + a, p + b, 0.25);
                            }
                        }
                        let next = [
                            (i + 1 < nx).then(|| point(i + 1, j, k)),
                            (j + 1 < ny).then(|| point(i, j + 1, k)),
                            (k + 1 < nz).then(|| point(i, j, k + 1)),
                        ];
                        for q in next.into_iter().flatten() {
                            for a in 0..3 {
                                pair(&mut rows, p + a, base + 3 * q + a, -1.0);
                            }
                        }
                    }
                }
            }
        }
        let s0 = 2 * per;
        for a in 0..sep {
            rows[s0 + a].push((s0 + a, 20.0));
            for b in a + 1..sep {
                pair(&mut rows, s0 + a, s0 + b, 0.5 / (1.0 + (b - a) as f64).powi(2));
            }
            pair(&mut rows, s0 + a, (7 * a) % per, -0.1);
            pair(&mut rows, s0 + a, per + (11 * a) % per, -0.1);
        }
        let d0 = s0 + sep;
        for a in 0..5 {
            rows[d0 + a].push((d0 + a, 10.0));
            for b in a + 1..5 {
                pair(&mut rows, d0 + a, d0 + b, 1.0 / (1.0 + (a + b) as f64));
            }
        }
        let mut ptr = vec![0];
        let mut idx = Vec::new();
        let mut val = Vec::new();
        for row in &mut rows {
            row.sort_by_key(|e| e.0);
            for (c, v) in row.iter() {
                idx.push(*c);
                val.push(*v);
            }
            ptr.push(idx.len());
        }
        let csr = Csr { n, ptr, idx, val };
        let upper = UpperColumns::from_symmetric(&csr);
        let one = parallel::with_threads(1, || Factor::new(&upper, usize::MAX)).unwrap();
        let supernodes = one.structure.first.len() - 1;
        let widest = (0..supernodes).map(|s| one.structure.cols(s).len()).max().unwrap();
        assert!(widest > 2 * PANEL, "a supernode spans panels: {widest}");
        check_solve(&csr, &one, 1e-12);
        let estimate = Factor::estimate_memory(&upper, usize::MAX).unwrap();
        for threads in [2, 7] {
            let (other, measured) =
                parallel::with_threads(threads, || Factor::with_memory_limit(&upper, usize::MAX, estimate.total_bytes())).unwrap();
            assert_eq!(measured, estimate, "the allocation estimate is independent of scheduling");
            assert!(one.values.iter().flatten().zip(other.values.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()), "{threads} threads");
        }
    }

    #[test]
    fn postorder_lists_children_first() {
        assert_eq!(postorder(&[2, 2, 4, 4, NONE, NONE]), vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(postorder(&[4, 3, 3, 4, NONE]), vec![0, 1, 2, 3, 4]);
        assert_eq!(postorder(&[NONE, 0]), vec![1, 0]);
        assert_eq!(postorder(&[3, NONE, 1, NONE]), vec![2, 1, 0, 3]);
    }
}
