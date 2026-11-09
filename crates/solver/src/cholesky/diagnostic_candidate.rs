use super::*;
pub(super) struct Target<'a> {
    pub m: usize,
    /// `values` starts at this column of the packed lower trapezoid.
    pub first: usize,
    pub relative: &'a [usize],
    pub values: &'a mut [f64],
}

struct Column<'a> {
    row: usize,
    values: &'a mut [f64],
}

struct Tile<'a> {
    row: usize,
    col: usize,
    columns: Vec<Column<'a>>,
}

pub(super) fn update(
    l: &Panel<'_>,
    cols: std::ops::Range<usize>,
    p: usize,
    mu: usize,
    wu: usize,
    target: Target<'_>,
    stats: &Stats,
    single: bool,
    isolated_kernel: bool,
) -> std::result::Result<(), FactorError> {
    let (k0, kw) = (cols.start, cols.len());
    let Target { m, first, relative, values } = target;
    if kw == 0 || wu == 0 || mu == 0 {
        return Ok(());
    }
    if wu > mu || relative.len() != mu || relative.windows(2).any(|r| r[0] >= r[1]) || relative.last().is_some_and(|&r| r >= m) {
        return Err(internal("tiled update row map"));
    }
    // Pack once per input row, shared read-only by every output column tile. Both
    // layouts keep k contiguous and use the same MR x NR micro-kernel as products.
    let packing = Instant::now();
    let mut xp = vec![[0.0; MR]; mu.div_ceil(TILE) * (TILE / MR) * kw];
    dispatch_chunks(&mut xp, (TILE / MR) * kw, single, stats, |t, packed| {
        let clock = stats.clock();
        let i0 = t * TILE;
        for k in 0..kw {
            for (i, x) in l.column(k0 + k, p + i0, TILE.min(mu - i0)).iter().enumerate() {
                if let Some(v) = packed.get_mut((i / MR) * kw + k).and_then(|r| r.get_mut(i % MR)) {
                    *v = *x;
                }
            }
        }
        stats.record_clock(0, clock, 1, (TILE.min(mu - i0) * kw * 16) as u64);
    });
    let mut yp = vec![[0.0; NR]; wu.div_ceil(NR) * kw];
    dispatch_chunks(&mut yp, (PANEL / NR) * kw, single, stats, |t, packed| {
        let clock = stats.clock();
        let c0 = t * PANEL;
        for k in 0..kw {
            for (c, y) in l.column(k0 + k, p + c0, PANEL.min(wu - c0)).iter().enumerate() {
                if let Some(v) = packed.get_mut((c / NR) * kw + k).and_then(|r| r.get_mut(c % NR)) {
                    *v = *y;
                }
            }
        }
        stats.record_clock(0, clock, 1, (PANEL.min(wu - c0) * kw * 16) as u64);
    });
    stats.pack_wall.lock().unwrap().push(packing.elapsed().as_secs_f64());
    let setup = stats.clock();
    let nc = wu.div_ceil(PANEL);
    let mut tiles: Vec<_> =
        (0..mu.div_ceil(TILE)).flat_map(|r| (0..nc).map(move |c| Tile { row: r * TILE, col: c * PANEL, columns: Vec::new() })).collect();
    // Split ownership before scheduling. Gaps in the sparse row map belong to one
    // segment but are never touched. Neither locks nor floating-point atomics are needed.
    let mut tail = values;
    let mut cursor = offset(m, first);
    for (c, &cc) in relative.iter().take(wu).enumerate() {
        let start = offset(m, cc);
        let skip = start.checked_sub(cursor).ok_or_else(|| internal("tiled update column order"))?;
        let (_, rest) = tail.split_at_mut_checked(skip).ok_or_else(|| internal("tiled update column offset"))?;
        let (column, rest) = rest.split_at_mut_checked(m - cc).ok_or_else(|| internal("tiled update column length"))?;
        tail = rest;
        cursor = start + m - cc;
        let mut segment = column;
        let mut row = cc;
        for r in c / TILE..mu.div_ceil(TILE) {
            let lo = (r * TILE).max(c);
            let hi = ((r + 1) * TILE).min(mu);
            let begin = relative.get(lo).copied().ok_or_else(|| internal("tiled update row start"))?;
            let end = relative.get(hi - 1).copied().ok_or_else(|| internal("tiled update row end"))? + 1;
            let (_, rest) = segment.split_at_mut_checked(begin - row).ok_or_else(|| internal("tiled update row offset"))?;
            let (values, rest) = rest.split_at_mut_checked(end - begin).ok_or_else(|| internal("tiled update row length"))?;
            segment = rest;
            row = end;
            let tile = tiles.get_mut(r * nc + c / PANEL).ok_or_else(|| internal("tiled update tile index"))?;
            tile.columns.push(Column { row: begin, values });
        }
    }
    stats.record_clock(2, setup, 1, (tiles.len() * std::mem::size_of::<Tile>()) as u64);
    dispatch_each(&mut tiles, single, stats, |_, tile| {
        let h = TILE.min(mu - tile.row);
        for b in (0..tile.columns.len()).step_by(NR) {
            let ys = yp.get((tile.col + b) / NR * kw..((tile.col + b) / NR + 1) * kw).unwrap_or_default();
            for r in (0..TILE).step_by(MR) {
                let xs = xp.get((tile.row + r) / MR * kw..((tile.row + r) / MR + 1) * kw).unwrap_or_default();
                let clock = stats.clock();
                let acc = if isolated_kernel {
                    fixed_kernel(xs, ys)
                } else {
                    let mut acc = [[0.0; MR]; NR];
                    for (x, y) in xs.iter().zip(ys) {
                        for (a, yj) in acc.iter_mut().zip(y) {
                            for (ai, xi) in a.iter_mut().zip(x) {
                                *ai += xi * yj;
                            }
                        }
                    }
                    acc
                };
                stats.record_clock(1, clock, 1, (kw * (MR + NR) * 8) as u64);
                let clock = stats.clock();
                for (j, sums) in acc.iter().enumerate() {
                    let c = tile.col + b + j;
                    if let Some(column) = tile.columns.get_mut(b + j) {
                        for (i, &sum) in sums.iter().enumerate().take(h.saturating_sub(r)) {
                            let ri = tile.row + r + i;
                            if ri >= c
                                && let Some(&mapped) = relative.get(ri)
                                && let Some(value) = column.values.get_mut(mapped.wrapping_sub(column.row))
                            {
                                *value -= sum;
                            }
                        }
                    }
                }
                stats.record_clock(3, clock, 1, (MR * NR * 16) as u64);
            }
        }
    });
    Ok(())
}

// Experimental code-generation boundary only: same MR/NR and ascending k, no FMA.
#[inline(never)]
fn fixed_kernel(xs: &[[f64; MR]], ys: &[[f64; NR]]) -> [[f64; MR]; NR] {
    let mut acc = [[0.0; MR]; NR];
    for (x, y) in xs.iter().zip(ys) {
        for (a, yj) in acc.iter_mut().zip(y) {
            for (ai, xi) in a.iter_mut().zip(x) {
                *ai += xi * yj;
            }
        }
    }
    acc
}
