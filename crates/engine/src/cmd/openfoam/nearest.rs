//! The OpenFOAM cell nearest each lattice cell. The search is bounded: its reach is derived from the
//! caller's distance, never from the size of the case, and the rings it visits are cut at the grid.

use fieldcraft_lbm::voxel::Grid;
use fieldcraft_model::progress;

use crate::EngineError;

/// Why the nearest cells could not be found: the input is refused, or a fluid lattice cell has no
/// OpenFOAM cell near it. The second carries its numbers in metres, for the caller to show in the
/// display units.
#[derive(Debug)]
pub(super) enum Failure {
    Refused(EngineError),
    Uncovered { position: [f64; 3], reach: f64 },
}

impl From<EngineError> for Failure {
    fn from(e: EngineError) -> Self {
        Self::Refused(e)
    }
}

impl std::fmt::Display for Failure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Refused(e) => write!(f, "{e}"),
            Self::Uncovered { position, reach } => write!(
                f,
                "no OpenFOAM cell centre within {reach:e} m of the lattice cell at {position:?} m: the case does not cover the step's domain there"
            ),
        }
    }
}

type Result<T, E = EngineError> = std::result::Result<T, E>;

/// A refusal of the case's input, with the reason it cannot be sampled.
fn refuse(why: impl Into<String>) -> EngineError {
    EngineError::invalid("path", why)
}

/// The OpenFOAM centres bucketed by lattice cell, in compressed sparse-row form: cell `c` holds
/// `ids[offsets[c]..offsets[c + 1]]`, and those ids are in increasing order. The memory is one `u32`
/// per cell plus one per centre, not a `Vec` per cell.
struct Buckets {
    offsets: Vec<u32>,
    ids: Vec<u32>,
}

impl Buckets {
    /// Buckets the centres that lie in the grid; centres outside it are skipped.
    fn new(grid: &Grid, centres: &[[f64; 3]]) -> Result<Self> {
        let cells = grid.len();
        let mut offsets = vec![0u32; cells + 1];
        // Pass 1: the count of centres per cell.
        for centre in centres {
            if !centre.iter().all(|v| v.is_finite()) {
                return Err(refuse("a cell centre that is not a finite point"));
            }
            if let Some(count) = grid.cell_at(*centre).and_then(|c| offsets.get_mut(c)) {
                // At most `centres.len()` counts, which the caller bounds by `u32::MAX`.
                *count += 1;
            }
        }
        // Inclusive prefix sums: `offsets[c]` becomes the end of cell `c`'s run.
        let mut end = 0u32;
        for count in offsets.iter_mut().take(cells) {
            end += *count;
            *count = end;
        }
        if let Some(total) = offsets.last_mut() {
            *total = end;
        }
        // Pass 2: fill each cell's run from its back, walking the centres in reverse, so each run
        // ends up in increasing id order and `offsets[c]` ends at the run's start.
        let mut ids = vec![0u32; end as usize];
        for (id, centre) in centres.iter().enumerate().rev() {
            let Some(c) = grid.cell_at(*centre) else {
                continue;
            };
            let Some(slot_end) = offsets.get_mut(c) else {
                return Err(refuse("a cell centre falls outside the bucket table"));
            };
            let Some(slot) = slot_end.checked_sub(1) else {
                return Err(refuse("a cell's bucket count is smaller than its centres"));
            };
            *slot_end = slot;
            let id = u32::try_from(id).map_err(|_| refuse("more OpenFOAM cells than a cell index holds"))?;
            let Some(out) = ids.get_mut(slot as usize) else {
                return Err(refuse("a cell centre falls outside the bucket table"));
            };
            *out = id;
        }
        Ok(Self { offsets, ids })
    }

    /// The ids of the centres in lattice cell `cell`, in increasing order.
    fn cell(&self, cell: usize) -> &[u32] {
        let (Some(&start), Some(&end)) = (self.offsets.get(cell), self.offsets.get(cell + 1)) else {
            return &[];
        };
        self.ids.get(start as usize..end as usize).unwrap_or(&[])
    }
}

/// For each lattice cell, the index of the OpenFOAM cell whose centre is nearest (`None` in a solid cell;
/// ties go to the lower index). A fluid cell with no centre within `reach` (metres) is an error, as is a
/// case with no centre inside the grid at all.
pub(super) fn nearest_cells(grid: &Grid, solid: &[bool], centres: &[[f64; 3]], reach: f64) -> Result<Vec<Option<u32>>, Failure> {
    if solid.len() != grid.len() {
        return Err(refuse(format!("the solid mask has {} cells but the lattice has {}", solid.len(), grid.len())).into());
    }
    if !(reach.is_finite() && reach > 0.0) {
        return Err(refuse(format!("the reach must be a finite, positive distance, not {reach:e} m")).into());
    }
    if !(grid.spacing.is_finite() && grid.spacing > 0.0) {
        return Err(refuse(format!("the lattice spacing must be a finite, positive length, not {:e} m", grid.spacing)).into());
    }
    if centres.len() > u32::MAX as usize {
        return Err(refuse("more OpenFOAM cells than a cell index holds").into());
    }
    let buckets = Buckets::new(grid, centres)?;
    if buckets.ids.is_empty() {
        return Err(refuse("no OpenFOAM cell centre lies in the step's domain: is the case another model's?").into());
    }
    // A centre within `reach` lies within `ceil(reach / spacing)` rings of a cell, and no ring reaches
    // past the widest side of the grid.
    let widest = grid.cells.iter().copied().max().unwrap_or(0);
    let rings = ((reach / grid.spacing).ceil() as usize).min(widest);

    let mut out = Vec::with_capacity(grid.len());
    for (cell, &is_solid) in solid.iter().enumerate() {
        if is_solid {
            out.push(None);
            continue;
        }
        match nearest_in_rings(grid, &buckets, centres, cell, rings)? {
            Some((d2, id)) if d2.sqrt() <= reach => out.push(Some(id)),
            _ => {
                return Err(Failure::Uncovered { position: grid.centre(cell), reach });
            }
        }
    }
    Ok(out)
}

/// The nearest centre to lattice cell `cell`, searched ring by ring out to `rings`, as `(d², id)`.
fn nearest_in_rings(grid: &Grid, buckets: &Buckets, centres: &[[f64; 3]], cell: usize, rings: usize) -> Result<Option<(f64, u32)>> {
    let [ci, cj, ck] = grid.coords(cell);
    let p = grid.centre(cell);
    let [nx, ny, nz] = grid.cells;
    let mut best: Option<(f64, u32)> = None;
    for r in 0..=rings {
        progress::check()?;
        // A centre in ring `r` or beyond is at least `(r - 1) · spacing` away, so a best nearer than
        // that cannot be beaten.
        if let Some((d2, _)) = best
            && d2.sqrt() < (r as f64 - 1.0) * grid.spacing
        {
            break;
        }
        let (x0, x1) = span(ci, r, nx);
        let (y0, y1) = span(cj, r, ny);
        let (z0, z1) = span(ck, r, nz);
        for k in z0..=z1 {
            for j in y0..=y1 {
                if k.abs_diff(ck) == r || j.abs_diff(cj) == r {
                    for i in x0..=x1 {
                        scan(grid, buckets, centres, p, [i, j, k], &mut best);
                    }
                } else {
                    // Only the two cells at distance `r` along x lie on the ring in this row.
                    for i in [ci.checked_sub(r), ci.checked_add(r)].into_iter().flatten() {
                        if (x0..=x1).contains(&i) {
                            scan(grid, buckets, centres, p, [i, j, k], &mut best);
                        }
                    }
                }
            }
        }
    }
    Ok(best)
}

/// The clipped range `[c − r, c + r]` within `0..n`.
fn span(c: usize, r: usize, n: usize) -> (usize, usize) {
    (c.saturating_sub(r), c.saturating_add(r).min(n.saturating_sub(1)))
}

/// Offers every centre bucketed in cell `ijk` to the running best.
fn scan(grid: &Grid, buckets: &Buckets, centres: &[[f64; 3]], p: [f64; 3], ijk: [usize; 3], best: &mut Option<(f64, u32)>) {
    for &id in buckets.cell(grid.index(ijk)) {
        let Some(q) = centres.get(id as usize) else {
            continue;
        };
        let d2 = distance_squared(*q, p);
        if beats(d2, id, *best) {
            *best = Some((d2, id));
        }
    }
}

/// Whether `(d2, id)` is nearer than the best so far: a smaller distance, or an equal one from a lower id.
fn beats(d2: f64, id: u32, best: Option<(f64, u32)>) -> bool {
    match best {
        None => true,
        Some((best_d2, best_id)) => d2 < best_d2 || (d2 == best_d2 && id < best_id),
    }
}

fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    (ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2)
}

#[cfg(test)]
mod tests {
    use super::{Failure, nearest_cells};
    use fieldcraft_lbm::voxel::Grid;

    fn row() -> Grid {
        Grid { origin: [0.0; 3], spacing: 1.0, cells: [5, 1, 1] }
    }

    fn refused(result: Result<Vec<Option<u32>>, Failure>, text: &str) {
        match result {
            Err(e) => assert!(e.to_string().contains(text), "wrong refusal: {e}"),
            Ok(v) => panic!("expected a refusal containing {text:?}, got {v:?}"),
        }
    }

    #[test]
    fn finds_the_nearest_centre_beyond_empty_cells() {
        let solid = [false, false, true, false, false];
        let centres = [[0.4, 0.5, 0.5], [4.9, 0.5, 0.5]];
        let got = nearest_cells(&row(), &solid, &centres, 5.0).expect("a covered row");
        assert_eq!(got, vec![Some(0), Some(0), None, Some(1), Some(1)]);
    }

    #[test]
    fn a_centre_in_the_last_ring_is_found() {
        // Cell 0's centre (0.5) is 4 from the only centre (4.5), in ring ceil(4 / 1) = 4, the last
        // one a reach of 4 searches; a reach just short of it refuses, never reads a zero.
        let fluid = [false; 5];
        let centres = [[4.5, 0.5, 0.5]];
        assert_eq!(nearest_cells(&row(), &fluid, &centres, 4.0).unwrap(), vec![Some(0); 5]);
        refused(nearest_cells(&row(), &fluid, &centres, 3.9), "within");
    }

    #[test]
    fn a_case_with_no_centre_in_the_domain_is_refused() {
        let solid = [false; 5];
        let centres = [[100.0, 0.5, 0.5]];
        refused(nearest_cells(&row(), &solid, &centres, 5.0), "no OpenFOAM cell centre lies in the step's domain");
    }

    #[test]
    fn a_fluid_cell_beyond_reach_is_an_error() {
        let solid = [false; 5];
        let centres = [[0.5, 0.5, 0.5]];
        refused(nearest_cells(&row(), &solid, &centres, 1.5), "within");
    }

    #[test]
    fn ties_go_to_the_lower_index() {
        // Cell 0's centre (0.5) is 0.25 from both centres, which lie inside the grid: the lower index
        // wins in either listing order, so it is the index and not the first listed.
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [2, 1, 1] };
        let (a, b) = ([0.75, 0.5, 0.5], [0.25, 0.5, 0.5]);
        for centres in [[a, b], [b, a]] {
            let got = nearest_cells(&grid, &[false, true], &centres, 1.0).expect("a covered cell");
            assert_eq!(got, vec![Some(0), None]);
        }
    }

    #[test]
    fn an_uncovered_cell_carries_its_position_and_reach() {
        let centres = [[0.5, 0.5, 0.5]];
        match nearest_cells(&row(), &[false; 5], &centres, 1.5) {
            Err(Failure::Uncovered { position, reach }) => assert_eq!((position, reach), ([2.5, 0.5, 0.5], 1.5)),
            other => panic!("expected an uncovered cell, got {other:?}"),
        }
    }

    #[test]
    fn hostile_inputs_never_panic() {
        let solid = [false; 5];
        let centres = [[0.5, 0.5, 0.5]];
        refused(nearest_cells(&row(), &solid, &[[f64::NAN, 0.5, 0.5]], 5.0), "not a finite point");
        for reach in [f64::NAN, 0.0, -1.0] {
            refused(nearest_cells(&row(), &solid, &centres, reach), "reach");
        }
        refused(nearest_cells(&row(), &[false; 4], &centres, 5.0), "solid mask");
        let bad = Grid { origin: [0.0; 3], spacing: f64::NAN, cells: [5, 1, 1] };
        refused(nearest_cells(&bad, &solid, &centres, 5.0), "spacing");
    }
}
