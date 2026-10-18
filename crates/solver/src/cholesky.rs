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

use crate::{Error, Result};

/// Relative pivot size below which the matrix counts as singular.
pub const PIVOT_TOLERANCE: f64 = 1e-10;

/// Rows per tile of the dense kernels (a tile's partial sums stay in registers and L1).
pub const TILE: usize = 64;

/// Columns per panel of a supernode's dense factorization.
pub const PANEL: usize = 64;

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
    /// The memory for the factor could not be had.
    OutOfMemory { entries: usize },
    /// A broken internal invariant (a bug), reported instead of a panic.
    Internal(String),
    /// The command was cancelled ([`fieldcraft_model::progress`]) between supernodes.
    Cancelled,
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
struct Update {
    k: usize,
    p1: usize,
    p2: usize,
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
    let tiles = mu.div_ceil(TILE);
    let tile = |t: usize| -> Vec<f64> {
        let i0 = t * TILE;
        let h = TILE.min(mu - i0);
        // The tile's rows of every column, padded to TILE rows.
        let mut x = vec![0.0; kw * TILE];
        for k in 0..kw {
            if let Some(dst) = x.get_mut(k * TILE..k * TILE + h) {
                dst.copy_from_slice(l.column(k0 + k, p + i0, h));
            }
        }
        let cols = wu.min(i0 + h);
        let mut out = vec![0.0; cols * TILE];
        for c in 0..cols {
            let mut acc = [0.0f64; TILE];
            let ys = yt.get(c * kw..(c + 1) * kw).unwrap_or_default();
            for (xs, b) in x.as_chunks::<TILE>().0.iter().zip(ys) {
                for (a, xi) in acc.iter_mut().zip(xs) {
                    *a += xi * b;
                }
            }
            if let Some(dst) = out.get_mut(c * TILE..(c + 1) * TILE) {
                dst.copy_from_slice(&acc);
            }
        }
        out
    };
    if mu.saturating_mul(wu).saturating_mul(kw) >= PARALLEL_WORK { parallel::map(tiles, tile) } else { (0..tiles).map(tile).collect() }
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
    fn new(a: &UpperColumns, inverse: &[usize]) -> Lower {
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

    fn column(&self, j: usize) -> (&[usize], &[f64]) {
        let (a, b) = (self.ptr.get(j).copied().unwrap_or(0), self.ptr.get(j + 1).copied().unwrap_or(0));
        (self.idx.get(a..b).unwrap_or_default(), self.val.get(a..b).unwrap_or_default())
    }
}

/// What a supernode's factorization reads: the matrix, the structure and the updates each
/// supernode receives.
pub(crate) struct Context<'a> {
    pub(crate) a: &'a Lower,
    pub(crate) s: &'a Structure,
    updates: &'a [Vec<Update>],
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
            let tiles = products(&Panel { m: kr.len(), data: lk }, 0..kw, u.p1, mu, wu, done_k.pivots());
            subtract(own, m, &tiles, mu, wu, |i| relative.get(i).copied().unwrap_or(0));
        }
        Ok(())
    }

    /// Factor supernode `j` into `own`; its descendants' values are `done[k - base]`. Fails with
    /// the first column (postorder) whose pivot is not positive.
    fn supernode(&self, j: usize, done: &[Vec<f64>], base: usize, own: &mut [f64]) -> std::result::Result<(), FactorError> {
        self.assemble(j, done, base, own)?;
        let rows = self.s.rows(j);
        let cols = self.s.cols(j);
        let (m, w, f) = (rows.len(), cols.len(), cols.start);
        // Dense factorization by panels.
        for c0 in (0..w).step_by(PANEL) {
            let c1 = (c0 + PANEL).min(w);
            let pw = c1 - c0;
            if c0 > 0 {
                // The panel minus its products with the columns before it.
                let tiles = products(&Panel { m, data: own }, 0..c0, c0, m - c0, pw, None);
                subtract(own, m, &tiles, m - c0, pw, |i| c0 + i);
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
                let solved = solve_below(&Panel { m, data: own }, c0, c1);
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
    pub(crate) fn factor<T: Supernode>(
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

    pub(crate) fn from_symbolic(symbolic: &Symbolic) -> std::result::Result<Factor, FactorError> {
        let mut values = symbolic.allocate()?;
        symbolic.factor(&mut values, |context, j, done, base, own| context.supernode(j, done, base, own))?;
        Ok(Factor { n: symbolic.n, perm: symbolic.perm.clone(), structure: symbolic.structure.clone(), values, entries: symbolic.entries })
    }

    /// Entries in L.
    pub fn entries(&self) -> usize {
        self.entries
    }

    /// Solve A x = b in place.
    pub fn solve(&self, b: &mut [f64]) -> Result<()> {
        solve_supernodes(self.n, &self.perm, &self.structure, &self.values, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Csr;

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
        for threads in [2, 7] {
            let other = parallel::with_threads(threads, || Factor::new(&upper, usize::MAX)).unwrap();
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
