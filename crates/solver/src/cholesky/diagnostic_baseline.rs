use super::*;
pub(super) fn baseline_products(
    l: &Panel,
    cols: std::ops::Range<usize>,
    p: usize,
    mu: usize,
    wu: usize,
    pivots: Option<&[Pivot]>,
    stats: &Stats,
    single: bool,
) -> Vec<Vec<f64>> {
    let (k0, kw) = (cols.start, cols.len());
    let packing = Instant::now();
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
    stats.record(0, packing, 1, (wu * kw * 32) as u64);
    let tiles = mu.div_ceil(TILE);
    let tile = |t: usize| -> Vec<f64> {
        let packing = Instant::now();
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
        stats.record(0, packing, 1, (h * kw * 16) as u64);
        let cols = wu.min(i0 + h);
        let mut out = vec![0.0; cols * TILE];
        // Micro-tiles of MR rows and NR columns: each S(i, c) is summed in a register over k in
        // order from 0, the same operations as one column at a time, with each X and Y value
        // loaded once per micro-tile.
        for (b, ys) in yp.chunks(kw.max(1)).enumerate().take(cols.div_ceil(NR)) {
            for (r, xs) in xp.chunks(kw.max(1)).enumerate() {
                let clock = stats.clock();
                let mut acc = [[0.0f64; MR]; NR];
                for (x, y) in xs.iter().zip(ys) {
                    for (a, yj) in acc.iter_mut().zip(y) {
                        for (ai, xi) in a.iter_mut().zip(x) {
                            *ai += xi * yj;
                        }
                    }
                }
                stats.record_clock(1, clock, 1, (kw * (MR + NR) * 8) as u64);
                let clock = stats.clock();
                for (j, a) in acc.iter().enumerate() {
                    let c = b * NR + j;
                    if c < cols
                        && let Some(dst) = out.get_mut(c * TILE + r * MR..c * TILE + (r + 1) * MR)
                    {
                        dst.copy_from_slice(a);
                    }
                }
                stats.record_clock(3, clock, 1, (MR * NR * 8) as u64);
            }
        }
        out
    };
    if mu.saturating_mul(wu).saturating_mul(kw) >= PARALLEL_WORK {
        dispatch_map(tiles, single, stats, tile)
    } else {
        dispatch_map(tiles, true, stats, tile)
    }
}
