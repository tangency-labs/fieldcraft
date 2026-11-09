//! Two-dimensional left-looking updates (Ng–Peyton 1993; Rothberg–Gupta 1994, cited
//! by the parent module). Each output tile owns disjoint mutable column segments.
//! Packing precedes the tile tasks; joining those tasks precedes the next descendant
//! update or panel factorization. Thus no tile reads an unfinished factor, and each
//! entry keeps the old ascending-k sum followed by exactly one subtraction per update.

use super::{FactorError, MR, NR, PANEL, Panel, TILE, internal, offset, parallel};

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

pub(super) fn update(l: &Panel<'_>, cols: std::ops::Range<usize>, p: usize, mu: usize, wu: usize, target: Target<'_>) -> Result<(), FactorError> {
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
    let mut xp = vec![[0.0; MR]; mu.div_ceil(TILE) * (TILE / MR) * kw];
    parallel::chunks_mut(&mut xp, (TILE / MR) * kw, |t, packed| {
        let i0 = t * TILE;
        for k in 0..kw {
            for (i, x) in l.column(k0 + k, p + i0, TILE.min(mu - i0)).iter().enumerate() {
                if let Some(v) = packed.get_mut((i / MR) * kw + k).and_then(|r| r.get_mut(i % MR)) {
                    *v = *x;
                }
            }
        }
    });
    let mut yp = vec![[0.0; NR]; wu.div_ceil(NR) * kw];
    parallel::chunks_mut(&mut yp, (PANEL / NR) * kw, |t, packed| {
        let c0 = t * PANEL;
        for k in 0..kw {
            for (c, y) in l.column(k0 + k, p + c0, PANEL.min(wu - c0)).iter().enumerate() {
                if let Some(v) = packed.get_mut((c / NR) * kw + k).and_then(|r| r.get_mut(c % NR)) {
                    *v = *y;
                }
            }
        }
    });
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
    parallel::for_each_mut(&mut tiles, |_, tile| {
        let h = TILE.min(mu - tile.row);
        for b in (0..tile.columns.len()).step_by(NR) {
            let ys = yp.get((tile.col + b) / NR * kw..((tile.col + b) / NR + 1) * kw).unwrap_or_default();
            for r in (0..TILE).step_by(MR) {
                let xs = xp.get((tile.row + r) / MR * kw..((tile.row + r) / MR + 1) * kw).unwrap_or_default();
                let mut acc = [[0.0; MR]; NR];
                for (x, y) in xs.iter().zip(ys) {
                    for (a, yj) in acc.iter_mut().zip(y) {
                        for (ai, xi) in a.iter_mut().zip(x) {
                            *ai += xi * yj;
                        }
                    }
                }
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
            }
        }
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cholesky::{products, subtract, trapezoid};

    fn exercise(mu: usize, wu: usize, kw: usize, gap: usize, threads: usize) -> (f64, f64) {
        let lm = kw + mu;
        let data: Vec<f64> = (0..trapezoid(lm, kw)).map(|i| ((i * 37 % 101) as f64 - 50.0) / 101.0).collect();
        let l = Panel { m: lm, data: &data };
        let relative: Vec<_> = (0..mu).map(|i| 3 + i * gap).collect();
        let m = relative[mu - 1] + 2;
        let mut reference = vec![0.75; trapezoid(m, relative[wu - 1] + 1)];
        let mut candidate = reference.clone();
        let (old, new) = parallel::with_threads(threads, || {
            let clock = std::time::Instant::now();
            let tiles = products(&l, 0..kw, kw, mu, wu, None);
            subtract(&mut reference, m, &tiles, mu, wu, |i| relative[i]);
            let old = clock.elapsed().as_secs_f64();
            drop(tiles);
            let clock = std::time::Instant::now();
            update(&l, 0..kw, kw, mu, wu, Target { m, first: 0, relative: &relative, values: &mut candidate }).unwrap();
            (old, clock.elapsed().as_secs_f64())
        });
        assert!(reference.iter().zip(&candidate).all(|(a, b)| a.to_bits() == b.to_bits()), "{mu}/{wu}/{kw}/{gap}, threads {threads}");
        (old, new)
    }

    #[test]
    fn disjoint_tiles_preserve_sparse_scatter_and_sum_bits() {
        for threads in [1, 2, 7] {
            for (mu, wu, kw, gap) in [(129, 67, 65, 2), (256, 192, 131, 1), (65, 1, 3, 3), (17, 17, 5, 1)] {
                exercise(mu, wu, kw, gap, threads);
            }
        }
    }

    #[test]
    #[ignore = "bounded representative dense-update benchmark; run explicitly in release"]
    fn representative_large_update() {
        for threads in [1, 2, 7, 32] {
            let (baseline, candidate) = exercise(4096, 2048, 512, 1, threads);
            eprintln!("TILED_UPDATE threads={threads} mu=4096 wu=2048 kw=512 baseline={baseline:.9} candidate={candidate:.9} bitwise=true");
        }
    }
}
