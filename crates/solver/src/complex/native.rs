//! Native complex symmetric P Z Pᵀ = L D Lᵀ, without conjugation.
//!
//! J. R. Bunch, *Stable algorithms for solving symmetric and skew-symmetric
//! systems*, Bull. Austral. Math. Soc. 26 (1982), 107–119, §4,
//! <https://doi.org/10.1017/S0004972700005621>, explicitly includes complex
//! symmetric matrices and 1×1/2×2 diagonal pivot blocks. The rook search of
//! Ashcraft, Grimes and Lewis (1998), DOI 10.1137/S0895479896296921, compares
//! complex moduli here: for a mutual column maximum z, |d11*d22| < α²|z|²,
//! so |det(D)| >= (1−α²)|z|². No real inertia is inferred from complex pivots.
//!
//! Reuses FieldCraft's Cholesky ordering, symbolic structure, subtree scheduler
//! and cancellation. Numeric panel/update loops specialize the existing real
//! LDLᵀ kernels to complex arithmetic. Pivots stay within each supernode;
//! a failed local pivot is a located refusal, not proof of a physical pole.
//! There is no positive-definiteness assumption about K−ω²M near resonance.
//! Complex values occupy two real slots; the existing byte-based factor cap
//! counts both. Each entry's sum order is independent of the thread count.

use super::Vector;
use crate::cholesky::{Context, FactorError, Lower, PANEL, PIVOT_TOLERANCE, Pivot, Structure, Symbolic, TILE, UpperColumns, offset, trapezoid};
use crate::sparse::Csr;
use crate::{Error, Result};
use fieldcraft_parallel as parallel;

const ALPHA: f64 = 0.640_388_203_202_207_6;
const ROOK_STEPS: usize = 64;
const MR: usize = 4;
const NR: usize = 2;
const PARALLEL_WORK: usize = 1 << 20;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct C {
    re: f64,
    im: f64,
}
impl C {
    const ZERO: Self = Self { re: 0.0, im: 0.0 };
    const ONE: Self = Self { re: 1.0, im: 0.0 };
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
}
impl std::ops::Add for C {
    type Output = Self;
    #[inline]
    fn add(self, b: Self) -> Self {
        Self { re: self.re + b.re, im: self.im + b.im }
    }
}
impl std::ops::Sub for C {
    type Output = Self;
    #[inline]
    fn sub(self, b: Self) -> Self {
        Self { re: self.re - b.re, im: self.im - b.im }
    }
}
impl std::ops::Mul for C {
    type Output = Self;
    #[inline]
    fn mul(self, b: Self) -> Self {
        Self { re: self.re * b.re - self.im * b.im, im: self.re * b.im + self.im * b.re }
    }
}
impl std::ops::Mul<&C> for &C {
    type Output = C;
    #[inline]
    fn mul(self, b: &C) -> C {
        *self * *b
    }
}
impl std::ops::Mul<C> for &C {
    type Output = C;
    #[inline]
    fn mul(self, b: C) -> C {
        *self * b
    }
}
impl std::ops::Div<f64> for C {
    type Output = Self;
    #[inline]
    fn div(self, b: f64) -> Self {
        Self { re: self.re / b, im: self.im / b }
    }
}
impl std::ops::Div for C {
    type Output = Self;
    #[inline]
    fn div(self, b: Self) -> Self {
        // Scale the divisor before squaring; no |b|² overflow/underflow.
        let scale = b.re.abs().max(b.im.abs());
        let (a, b) = (self / scale, b / scale);
        let d = b.re * b.re + b.im * b.im;
        Self { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d }
    }
}
impl std::ops::AddAssign for C {
    #[inline]
    fn add_assign(&mut self, b: Self) {
        *self = *self + b;
    }
}
impl std::ops::SubAssign for C {
    #[inline]
    fn sub_assign(&mut self, b: Self) {
        *self = *self - b;
    }
}
impl std::ops::SubAssign<&C> for C {
    #[inline]
    fn sub_assign(&mut self, b: &C) {
        *self = *self - *b;
    }
}
impl std::ops::MulAssign for C {
    #[inline]
    fn mul_assign(&mut self, b: Self) {
        *self = *self * b;
    }
}
impl std::ops::DivAssign for C {
    #[inline]
    fn div_assign(&mut self, b: Self) {
        *self = *self / b;
    }
}

/// Symmetric two-pivot solve, scaled before forming the determinant.
fn two_solve(d11: C, d21: C, d22: C, x: C, y: C) -> (C, C) {
    let scale = d11.abs().max(d21.abs()).max(d22.abs());
    let (a, b, c, x, y) = (d11 / scale, d21 / scale, d22 / scale, x / scale, y / scale);
    let det = a * c - b * b;
    ((c * x - b * y) / det, (a * y - b * x) / det)
}

#[derive(Debug, Clone, Default)]
struct Block {
    values: Vec<C>,
    pivots: Vec<Pivot>,
    order: Vec<usize>,
}
impl Block {
    fn values(&self) -> &[C] {
        &self.values
    }
    fn pivots(&self) -> Option<&[Pivot]> {
        Some(&self.pivots)
    }
    fn order(&self) -> Option<&[usize]> {
        Some(&self.order)
    }
}
enum Choice {
    One(usize),
    Two(usize, usize),
}
struct Scratch {
    a: Vec<C>,
    b: Vec<C>,
    w: Vec<C>,
}
struct Panel<'a> {
    m: usize,
    data: &'a [C],
}

impl Panel<'_> {
    /// L(r, c) for r ≥ c.
    fn at(&self, r: usize, c: usize) -> C {
        self.data.get(offset(self.m, c) + r - c).copied().unwrap_or(C::ZERO)
    }

    /// Rows r..r+len of column c (r ≥ c).
    fn column(&self, c: usize, r: usize, len: usize) -> &[C] {
        let start = offset(self.m, c) + r - c;
        self.data.get(start..start + len).unwrap_or_default()
    }
}

fn times_d(row: &mut [C], l: &Panel, k0: usize, pivots: &[Pivot]) {
    let mut k = 0;
    while k < row.len() {
        let d = |r: usize, c: usize| l.at(k0 + r, k0 + c);
        if pivots.get(k) == Some(&Pivot::First) {
            let (a, b) = (row.get(k).copied().unwrap_or(C::ZERO), row.get(k + 1).copied().unwrap_or(C::ZERO));
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
fn products(l: &Panel, cols: std::ops::Range<usize>, p: usize, mu: usize, wu: usize, pivots: Option<&[Pivot]>) -> Vec<Vec<C>> {
    let (k0, kw) = (cols.start, cols.len());
    // Rows p..p+wu of the columns, transposed so that column c's coefficients are contiguous.
    let mut yt = vec![C::ZERO; wu * kw];
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
    let mut yp = vec![[C::ZERO; NR]; wu.div_ceil(NR) * kw];
    for (c, row) in yt.chunks(kw.max(1)).enumerate().take(wu) {
        for (k, y) in row.iter().enumerate() {
            if let Some(slot) = yp.get_mut((c / NR) * kw + k).and_then(|b| b.get_mut(c % NR)) {
                *slot = *y;
            }
        }
    }
    let tiles = mu.div_ceil(TILE);
    let tile = |t: usize| -> Vec<C> {
        let i0 = t * TILE;
        let h = TILE.min(mu - i0);
        // The tile's rows of every column by blocks of MR rows, k-major, padded to TILE rows.
        let mut xp = vec![[C::ZERO; MR]; TILE / MR * kw];
        for k in 0..kw {
            for (i, x) in l.column(k0 + k, p + i0, h).iter().enumerate() {
                if let Some(slot) = xp.get_mut((i / MR) * kw + k).and_then(|b| b.get_mut(i % MR)) {
                    *slot = *x;
                }
            }
        }
        let cols = wu.min(i0 + h);
        let mut out = vec![C::ZERO; cols * TILE];
        // Micro-tiles of MR rows and NR columns: each S(i, c) is summed in a register over k in
        // order from 0, the same operations as one column at a time, with each X and Y value
        // loaded once per micro-tile.
        for (b, ys) in yp.chunks(kw.max(1)).enumerate().take(cols.div_ceil(NR)) {
            for (r, xs) in xp.chunks(kw.max(1)).enumerate() {
                let mut acc = [[C::ZERO; MR]; NR];
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

/// own(rel(i), rel(c)) −= S(i, c) for i ≥ c, with S in the tiles of [`products`].
fn subtract(own: &mut [C], m: usize, tiles: &[Vec<C>], mu: usize, wu: usize, rel: impl Fn(usize) -> usize) {
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

fn interchange(v: &mut [C], m: usize, p: usize, q: usize) {
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

fn largest(v: &[C], range: std::ops::Range<usize>, skip: Option<usize>) -> (usize, f64) {
    let mut best = (range.start, 0.0f64);
    for i in range {
        if Some(i) == skip {
            continue;
        }
        let x = v.get(i).copied().unwrap_or(C::ZERO).abs();
        if x.is_nan() {
            return (i, f64::NAN);
        }
        if x > best.1 {
            best = (i, x);
        }
    }
    best
}

impl Block {
    /// Column `c` (local, ≥ k) of the part still to factor, rows k..m, into `out`: the trapezoid's
    /// entries less the products with the panel's factored columns c0..k.
    fn current(&self, m: usize, c0: usize, k: usize, c: usize, out: &mut [C], w: &mut Vec<C>) {
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
    fn one(&mut self, m: usize, k: usize, col: &[C]) {
        let d = col.get(k).copied().unwrap_or(C::ZERO);
        let start = offset(m, k);
        if let Some(x) = self.values.get_mut(start) {
            *x = d;
        }
        for r in k + 1..m {
            if let Some(x) = self.values.get_mut(start + r - k) {
                *x = col.get(r).copied().unwrap_or(C::ZERO) / d;
            }
        }
        if let Some(p) = self.pivots.get_mut(k) {
            *p = Pivot::One;
        }
    }

    /// Store a 2 × 2 pivot at k, k + 1 from the current columns `a` (k) and `b` (k + 1).
    fn two(&mut self, m: usize, k: usize, a: &[C], b: &[C]) {
        let get = |v: &[C], r: usize| v.get(r).copied().unwrap_or(C::ZERO);
        let (d11, d21, d22) = (get(a, k), get(a, k + 1), get(b, k + 1));
        let (s0, s1) = (offset(m, k), offset(m, k + 1));
        for (at, v) in [(s0, d11), (s0 + 1, d21), (s1, d22)] {
            if let Some(x) = self.values.get_mut(at) {
                *x = v;
            }
        }
        for r in k + 2..m {
            let (x, y) = two_solve(d11, d21, d22, get(a, r), get(b, r));
            if let Some(l) = self.values.get_mut(s0 + r - k) {
                *l = x;
            }
            if let Some(l) = self.values.get_mut(s1 + r - k - 1) {
                *l = y;
            }
        }
        if let Some(p) = self.pivots.get_mut(k) {
            *p = Pivot::First;
        }
        if let Some(p) = self.pivots.get_mut(k + 1) {
            *p = Pivot::Second;
        }
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
            let akk = scratch.a.get(k).copied().unwrap_or(C::ZERO).abs();
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
                    let arr = scratch.b.get(r).copied().unwrap_or(C::ZERO).abs();
                    if !arr.is_finite() || !omega_r.is_finite() {
                        return Err(singular(self, k));
                    }
                    step += 1;
                    if arr >= ALPHA * omega_r {
                        std::mem::swap(&mut scratch.a, &mut scratch.b);
                        break Choice::One(r);
                    }
                    if omega_r <= omega_i {
                        break Choice::Two(i, r);
                    }
                    if step >= ROOK_STEPS {
                        return Err(singular(self, k));
                    }
                    std::mem::swap(&mut scratch.a, &mut scratch.b);
                    (i, omega_i, r) = (r, omega_r, s);
                }
            };
            match choice {
                Choice::One(p) => {
                    self.swap(m, k, p, scratch);
                    if negligible(self, k, scratch.a.get(k).copied().unwrap_or(C::ZERO).abs()) {
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
                    let offdiagonal = scratch.a.get(k + 1).copied().unwrap_or(C::ZERO).abs();
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

fn solve_supernodes(n: usize, perm: &[usize], s: &Structure, blocks: &[Block], b: &mut [C]) -> Result<()> {
    if b.len() != n {
        return Err(Error::Internal("right-hand side length".into()));
    }
    let mut x: Vec<C> = perm.iter().map(|&j| b.get(j).copied().unwrap_or(C::ZERO)).collect();
    let get = |v: &[C], i: usize| v.get(i).copied().unwrap_or(C::ZERO);
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
                yq /= t.data.get(offset(m, q)).copied().unwrap_or(C::ONE);
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
                    let (p, q) = two_solve(d11, d21, d22, get(y, c), get(y, c + 1));
                    if let Some(v) = y.get_mut(c) {
                        *v = p;
                    }
                    if let Some(v) = y.get_mut(c + 1) {
                        *v = q;
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
                sum /= t.data.get(offset(m, q)).copied().unwrap_or(C::ONE);
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
    FactorError::Internal(format!("complex supernodal factor: {what}"))
}

fn assemble(context: &Context, imaginary: &Lower, j: usize, done: &[Block], base: usize, own: &mut [C]) -> std::result::Result<(), FactorError> {
    let rows = context.s.rows(j);
    let cols = context.s.cols(j);
    let (m, w, f) = (rows.len(), cols.len(), cols.start);
    if own.len() != trapezoid(m, w) {
        return Err(internal("supernode size"));
    }
    own.fill(C::ZERO);
    for c in 0..w {
        let (ai, av) = context.a.column(f + c);
        let (bi, bv) = imaginary.column(f + c);
        if ai != bi {
            return Err(internal("real/imaginary pattern"));
        }
        for ((r, re), im) in ai.iter().zip(av).zip(bv) {
            let at = rows.binary_search(r).map_err(|_| internal("entry outside structure"))?;
            if at < c {
                return Err(internal("entry above diagonal"));
            }
            if let Some(x) = own.get_mut(offset(m, c) + at - c) {
                *x = C { re: *re, im: *im };
            }
        }
    }
    let mut relative = Vec::new();
    for u in context.updates.get(j).map(Vec::as_slice).unwrap_or_default() {
        let blk = done.get(u.k.wrapping_sub(base)).ok_or_else(|| internal("update outside subtree"))?;
        let (kr, kw) = (context.s.rows(u.k), context.s.cols(u.k).len());
        let below = kr.get(u.p1..).unwrap_or_default();
        let (mu, wu) = (below.len(), u.p2 - u.p1);
        relative.clear();
        let mut at = 0;
        for r in below {
            while rows.get(at).is_some_and(|x| x < r) {
                at += 1;
            }
            if rows.get(at) != Some(r) {
                return Err(internal("update row outside structure"));
            }
            relative.push(at);
        }
        let tiles = products(&Panel { m: kr.len(), data: &blk.values }, 0..kw, u.p1, mu, wu, Some(&blk.pivots));
        subtract(own, m, &tiles, mu, wu, |i| relative.get(i).copied().unwrap_or(0));
    }
    Ok(())
}

fn supernode(context: &Context, imaginary: &Lower, j: usize, done: &[Block], base: usize, own: &mut Block) -> std::result::Result<(), FactorError> {
    assemble(context, imaginary, j, done, base, &mut own.values)?;
    let (m, cols) = (context.s.rows(j).len(), context.s.cols(j));
    let (w, f) = (cols.len(), cols.start);
    own.pivots = vec![Pivot::One; w];
    own.order = (0..w).collect();
    // The norm of a complex row is at least the maximum of its component norms.
    let scale = |local: usize| context.a.scale.get(f + local).copied().unwrap_or(0.0).max(imaginary.scale.get(f + local).copied().unwrap_or(0.0));
    let column = |local: usize| f + local;
    let mut scratch = Scratch { a: vec![C::ZERO; m], b: vec![C::ZERO; m], w: Vec::with_capacity(PANEL + 1) };
    let mut c0 = 0;
    while c0 < w {
        let c1 = own.panel(m, w, c0, &scale, &column, &mut scratch)?;
        if c1 < w {
            let tiles = products(&Panel { m, data: &own.values }, c0..c1, c1, m - c1, w - c1, own.pivots.get(c0..c1));
            subtract(&mut own.values, m, &tiles, m - c1, w - c1, |i| c1 + i);
        }
        c0 = c1;
    }
    if own.values.iter().any(|v| !v.is_finite()) {
        return Err(FactorError::Singular { column: f });
    }
    Ok(())
}

pub(super) struct Factor {
    n: usize,
    perm: Vec<usize>,
    structure: Structure,
    blocks: Vec<Block>,
    entries: usize,
}
impl Factor {
    pub(super) fn new(a: &Csr, b: &Csr, max_entries: usize) -> std::result::Result<Self, FactorError> {
        // Symbolic checks before numeric allocation, charged in real f64 slots.
        let symbolic = Symbolic::new(&UpperColumns::from_symmetric(a), max_entries / 2).map_err(|e| match e {
            FactorError::TooLarge { entries, .. } => FactorError::TooLarge { entries: entries.saturating_mul(2), limit: max_entries },
            e => e,
        })?;
        let entries = symbolic.entries.checked_mul(2).ok_or_else(|| internal("factor size overflow"))?;
        let mut inverse = vec![0; symbolic.n];
        for (k, &j) in symbolic.perm.iter().enumerate() {
            if let Some(v) = inverse.get_mut(j) {
                *v = k;
            }
        }
        let imaginary = Lower::new(&UpperColumns::from_symmetric(b), &inverse);
        let mut blocks = Vec::new();
        blocks.try_reserve_exact(symbolic.supernodes()).map_err(|_| FactorError::OutOfMemory { entries })?;
        for s in 0..symbolic.supernodes() {
            let mut values = Vec::new();
            values.try_reserve_exact(symbolic.structure.size(s)).map_err(|_| FactorError::OutOfMemory { entries })?;
            values.resize(symbolic.structure.size(s), C::ZERO);
            blocks.push(Block { values, ..Block::default() });
        }
        symbolic.factor(&mut blocks, |context, j, done, base, own| supernode(context, &imaginary, j, done, base, own))?;
        Ok(Self { n: symbolic.n, perm: symbolic.perm.clone(), structure: symbolic.structure.clone(), blocks, entries })
    }
    pub(super) fn entries(&self) -> usize {
        self.entries
    }
    pub(super) fn solve(&self, rhs: &mut Vector) -> Result<()> {
        if rhs.real.len() != self.n || rhs.imaginary.len() != self.n {
            return Err(Error::Invalid("complex factor right-hand side length".into()));
        }
        let mut values: Vec<C> = rhs.real.iter().zip(&rhs.imaginary).map(|(&re, &im)| C { re, im }).collect();
        solve_supernodes(self.n, &self.perm, &self.structure, &self.blocks, &mut values)?;
        for ((r, i), v) in rhs.real.iter_mut().zip(&mut rhs.imaginary).zip(values) {
            if !v.is_finite() {
                return Err(Error::Invalid("non-finite native complex solve".into()));
            }
            *r = v.re;
            *i = v.im;
        }
        Ok(())
    }
}
