//! Fill-reducing ordering: nested dissection by coordinate bisection.
//!
//! Nested dissection (A. George, *Nested dissection of a regular finite element mesh*, SIAM J.
//! Numer. Anal. 10(2):345–363, 1973) orders a separator after the two parts it splits, so the
//! Cholesky factor of a 3D mesh with n nodes fills O(n^{4/3}) instead of the O(n^{5/3}) of a band
//! or envelope ordering. The separator here is geometric: the nodes are split at the median of
//! the longest axis of their bounding box, and the nodes of the smaller side that touch the other
//! side form the separator (the coordinate nested dissection of H. D. Simon, *Partitioning of
//! unstructured problems for parallel processing*, Computing Systems in Engineering 2:135–148,
//! 1991, with vertex separators). Ties are broken by node index, so the ordering is
//! deterministic.
//!
//! A plane cuts a thin curved wall (a shell meshed with solids, a bellows, a revolved part) many
//! times, so its separators are large there. [`graph_dissection`] takes its separators from the
//! graph instead: a middle level of the breadth-first level structure rooted at a
//! pseudo-peripheral node (A. George, J. W. H. Liu, *An automatic nested dissection algorithm for
//! irregular finite element problems*, SIAM J. Numer. Anal. 15(5):1053–1069, 1978; the
//! pseudo-peripheral node of N. E. Gibbs, W. G. Poole, P. K. Stockmeyer, *An algorithm for
//! reducing the bandwidth and profile of a sparse matrix*, SIAM J. Numer. Anal. 13(2):236–250,
//! 1976), which cuts such a wall across its surface. [`fill_reducing_order`] computes both and
//! keeps the one whose Cholesky factor has fewer entries; the symbolic count costs O(|L|) time
//! and O(n + |A|) memory, small next to the numeric factorization.

/// Parts this small are ordered as they come.
const LEAF: usize = 16;

/// A permutation of the nodes: `order[k]` is the node eliminated k-th.
pub fn nested_dissection(xyz: &[[f64; 3]], adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = xyz.len();
    let mut order = Vec::with_capacity(n);
    // 0 = not in the current part, 1 = left, 2 = right.
    let mut side = vec![0u8; n];
    // An explicit stack instead of recursion: (part, separator to emit after it).
    enum Task {
        Split(Vec<usize>),
        Emit(Vec<usize>),
    }
    let mut stack = vec![Task::Split((0..n).collect())];
    while let Some(task) = stack.pop() {
        let mut part = match task {
            Task::Emit(sep) => {
                order.extend(sep);
                continue;
            }
            Task::Split(p) => p,
        };
        if part.len() <= LEAF {
            order.extend(part);
            continue;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &v in &part {
            let p = xyz.get(v).copied().unwrap_or_default();
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let axis = (0..3).fold(0, |best, k| if hi[k] - lo[k] > hi[best] - lo[best] { k } else { best });
        let key = |v: usize| xyz.get(v).map_or(0.0, |p| p[axis]);
        part.sort_by(|&a, &b| key(a).total_cmp(&key(b)).then(a.cmp(&b)));
        let half = part.len() / 2;
        for (i, &v) in part.iter().enumerate() {
            if let Some(s) = side.get_mut(v) {
                *s = if i < half { 1 } else { 2 };
            }
        }
        let touches = |v: usize, other: u8| adjacency.get(v).is_some_and(|row| row.iter().any(|w| side.get(*w) == Some(&other)));
        let (left, right) = part.split_at(half);
        let left_boundary: Vec<usize> = left.iter().copied().filter(|&v| touches(v, 2)).collect();
        let right_boundary: Vec<usize> = right.iter().copied().filter(|&v| touches(v, 1)).collect();
        let separator = if left_boundary.len() <= right_boundary.len() { left_boundary } else { right_boundary };
        for &v in &part {
            if let Some(s) = side.get_mut(v) {
                *s = 0;
            }
        }
        // Mark the separator so the parts leave it out.
        for &v in &separator {
            if let Some(s) = side.get_mut(v) {
                *s = 3;
            }
        }
        let a: Vec<usize> = left.iter().copied().filter(|v| side.get(*v) != Some(&3)).collect();
        let b: Vec<usize> = right.iter().copied().filter(|v| side.get(*v) != Some(&3)).collect();
        for &v in &separator {
            if let Some(s) = side.get_mut(v) {
                *s = 0;
            }
        }
        // Order: a, then b, then the separator (the stack runs in reverse).
        stack.push(Task::Emit(separator));
        stack.push(Task::Split(b));
        stack.push(Task::Split(a));
    }
    order
}

/// The nested dissection order of the nodes that fills the factor less: by coordinate planes or
/// by the graph's level structures. Ties keep the coordinate order.
pub fn fill_reducing_order(xyz: &[[f64; 3]], adjacency: &[Vec<usize>]) -> Vec<usize> {
    let by_planes = nested_dissection(xyz, adjacency);
    if adjacency.len() <= LEAF {
        // Both order a leaf as it comes: there is no choice to make.
        return by_planes;
    }
    let by_graph = graph_dissection(adjacency);
    if fill(adjacency, &by_graph) < fill(adjacency, &by_planes) { by_graph } else { by_planes }
}

/// Searches for a pseudo-peripheral node per part.
const PERIPHERAL_SEARCHES: usize = 8;

/// The breadth-first level structure of the nodes with `owner == token` reached from `root`.
/// Marks each reached node's level in `level` (which must be `usize::MAX` there on entry).
fn level_structure(adjacency: &[Vec<usize>], root: usize, token: usize, owner: &[usize], level: &mut [usize]) -> Vec<Vec<usize>> {
    let mut out = vec![vec![root]];
    if let Some(l) = level.get_mut(root) {
        *l = 0;
    }
    loop {
        let depth = out.len();
        let mut next = Vec::new();
        for &v in out.last().into_iter().flatten() {
            for &w in adjacency.get(v).into_iter().flatten() {
                if owner.get(w) == Some(&token)
                    && let Some(l) = level.get_mut(w)
                    && *l == usize::MAX
                {
                    *l = depth;
                    next.push(w);
                }
            }
        }
        if next.is_empty() {
            return out;
        }
        out.push(next);
    }
}

/// Clears the levels of `structure`'s nodes.
fn clear(level: &mut [usize], structure: &[Vec<usize>]) {
    for &v in structure.iter().flatten() {
        if let Some(l) = level.get_mut(v) {
            *l = usize::MAX;
        }
    }
}

/// A permutation of the nodes by nested dissection on the graph alone (George & Liu 1978, see
/// the module notes): each part's separator is the nodes of its middle level that touch the
/// next level. Disconnected parts are ordered one after the other. Ties go by node index, so
/// the ordering is deterministic.
pub fn graph_dissection(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = adjacency.len();
    let mut order = Vec::with_capacity(n);
    // The part each node is in (a token per part), and its level in the current search.
    let mut owner = vec![0usize; n];
    let mut level = vec![usize::MAX; n];
    let mut token = 0usize;
    enum Task {
        Split(Vec<usize>),
        Emit(Vec<usize>),
    }
    let mut stack = vec![Task::Split((0..n).collect())];
    while let Some(task) = stack.pop() {
        let part = match task {
            Task::Emit(sep) => {
                order.extend(sep);
                continue;
            }
            Task::Split(p) => p,
        };
        let Some(&first) = part.first() else { continue };
        if part.len() <= LEAF {
            order.extend(part);
            continue;
        }
        token += 1;
        for &v in &part {
            if let Some(o) = owner.get_mut(v) {
                *o = token;
            }
        }
        let degree = |v: usize, owner: &[usize]| adjacency.get(v).map_or(0, |row| row.iter().filter(|w| owner.get(**w) == Some(&token)).count());
        // A pseudo-peripheral root: from the first node, move to the least connected node of the
        // last level while that deepens the structure.
        let mut structure = level_structure(adjacency, first, token, &owner, &mut level);
        for _ in 0..PERIPHERAL_SEARCHES {
            let Some(candidate) = structure.last().and_then(|last| last.iter().copied().min_by_key(|&v| (degree(v, &owner), v))) else { break };
            clear(&mut level, &structure);
            let deeper = level_structure(adjacency, candidate, token, &owner, &mut level);
            if deeper.len() <= structure.len() {
                // Keep the structure found: restore its levels.
                clear(&mut level, &deeper);
                for (depth, l) in structure.iter().enumerate() {
                    for &v in l {
                        if let Some(x) = level.get_mut(v) {
                            *x = depth;
                        }
                    }
                }
                break;
            }
            structure = deeper;
        }
        let reached: usize = structure.iter().map(Vec::len).sum();
        if reached < part.len() {
            // Disconnected: the reached component, then the rest.
            let rest: Vec<usize> = part.iter().copied().filter(|&v| level.get(v) == Some(&usize::MAX)).collect();
            clear(&mut level, &structure);
            let mut component: Vec<usize> = structure.into_iter().flatten().collect();
            component.sort_unstable();
            stack.push(Task::Split(rest));
            stack.push(Task::Split(component));
            continue;
        }
        if structure.len() < 3 {
            clear(&mut level, &structure);
            order.extend(part);
            continue;
        }
        // The middle level: the first at which the levels up to it hold half the part, kept
        // off the last level so both sides are non-empty.
        let mut count = 0;
        let mut middle = 1;
        for (i, l) in structure.iter().enumerate().skip(1) {
            count += l.len();
            if 2 * count >= part.len() {
                middle = i.min(structure.len() - 2);
                break;
            }
        }
        let touches_next =
            |v: usize| adjacency.get(v).is_some_and(|row| row.iter().any(|&w| owner.get(w) == Some(&token) && level.get(w) == Some(&(middle + 1))));
        let mut separator = Vec::new();
        let mut before: Vec<usize> = structure.iter().take(middle).flatten().copied().collect();
        for &v in structure.get(middle).into_iter().flatten() {
            if touches_next(v) { separator.push(v) } else { before.push(v) }
        }
        let mut after: Vec<usize> = structure.iter().skip(middle + 1).flatten().copied().collect();
        clear(&mut level, &structure);
        before.sort_unstable();
        separator.sort_unstable();
        after.sort_unstable();
        // Order: before, then after, then the separator (the stack runs in reverse).
        stack.push(Task::Emit(separator));
        stack.push(Task::Split(after));
        stack.push(Task::Split(before));
    }
    order
}

/// The entries of the Cholesky factor of the node graph `adjacency` eliminated in `order` (one
/// per node pair): the fill an ordering gives. `order` must be a permutation of the nodes.
pub fn fill(adjacency: &[Vec<usize>], order: &[usize]) -> usize {
    let n = adjacency.len();
    let mut position = vec![usize::MAX; n];
    for (k, &v) in order.iter().enumerate() {
        if let Some(p) = position.get_mut(v) {
            *p = k;
        }
    }
    let mut ptr = Vec::with_capacity(n + 1);
    ptr.push(0);
    let mut idx = Vec::new();
    let mut rows = Vec::new();
    for (k, &v) in order.iter().enumerate() {
        rows.clear();
        rows.extend(adjacency.get(v).into_iter().flatten().filter_map(|w| position.get(*w).copied()).filter(|&r| r < k));
        rows.push(k);
        rows.sort_unstable();
        rows.dedup();
        idx.extend_from_slice(&rows);
        ptr.push(idx.len());
    }
    let val = vec![0.0; idx.len()];
    crate::cholesky::factor_entries(&crate::cholesky::UpperColumns { n: order.len(), ptr, idx, val })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_permutation_with_separators_last() {
        // A 12 × 12 grid graph.
        let m = 12;
        let id = |i: usize, j: usize| i * m + j;
        let xyz: Vec<[f64; 3]> = (0..m * m).map(|v| [(v / m) as f64, (v % m) as f64, 0.0]).collect();
        let mut adj = vec![Vec::new(); m * m];
        for i in 0..m {
            for j in 0..m {
                for (di, dj) in [(0i64, 1i64), (1, 0), (0, -1), (-1, 0), (0, 0)] {
                    let (a, b) = (i as i64 + di, j as i64 + dj);
                    if (0..m as i64).contains(&a) && (0..m as i64).contains(&b) {
                        adj[id(i, j)].push(id(a as usize, b as usize));
                    }
                }
            }
        }
        let order = nested_dissection(&xyz, &adj);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..m * m).collect::<Vec<_>>());
        // The top separator (last 12 entries) is one grid line.
        let last: Vec<usize> = order[m * m - m..].to_vec();
        assert!(last.iter().all(|v| v / m == last[0] / m) || last.iter().all(|v| v % m == last[0] % m), "{last:?}");
        assert_eq!(order, nested_dissection(&xyz, &adj), "deterministic");
    }

    /// The node graph of a structured hex8 grid of `dims` cells (27-point stencil), with the
    /// node coordinates from `at(i, j, k)`; `wrap` closes the first index around.
    fn hex_grid(dims: [usize; 3], wrap: bool, at: impl Fn(usize, usize, usize) -> [f64; 3]) -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
        let ni = if wrap { dims[0] } else { dims[0] + 1 };
        let (nj, nk) = (dims[1] + 1, dims[2] + 1);
        let id = |i: usize, j: usize, k: usize| (i % ni) + ni * (j + nj * k);
        let mut xyz = vec![[0.0; 3]; ni * nj * nk];
        let mut adj = vec![Vec::new(); ni * nj * nk];
        for k in 0..nk {
            for j in 0..nj {
                for i in 0..ni {
                    xyz[id(i, j, k)] = at(i, j, k);
                    for t in 0..27 {
                        let (di, dj, dk) = (t % 3, (t / 3) % 3, t / 9);
                        let (a, b, c) = (i + ni + di - 1, j + dj, k + dk);
                        let a = if wrap {
                            a % ni
                        } else if a < ni || a >= 2 * ni {
                            continue;
                        } else {
                            a - ni
                        };
                        if (1..=nj).contains(&b) && (1..=nk).contains(&c) {
                            adj[id(i, j, k)].push(id(a, b - 1, c - 1));
                        }
                    }
                }
            }
        }
        for row in &mut adj {
            row.sort_unstable();
            row.dedup();
        }
        (xyz, adj)
    }

    fn is_permutation(order: &[usize], n: usize) -> bool {
        let mut sorted = order.to_vec();
        sorted.sort_unstable();
        sorted == (0..n).collect::<Vec<_>>()
    }

    /// A thin cylinder wall, two solid elements thick, 160 around and 40 high.
    fn cylinder_wall() -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
        hex_grid([160, 2, 40], true, |i, j, k| {
            let t = i as f64 / 160.0 * std::f64::consts::TAU;
            let r = 50.0 + j as f64 * 0.5;
            [r * t.cos(), r * t.sin(), k as f64 * 2.0]
        })
    }

    #[test]
    fn graph_dissection_is_a_deterministic_permutation() {
        let (_, adj) = cylinder_wall();
        let order = graph_dissection(&adj);
        assert!(is_permutation(&order, adj.len()));
        assert_eq!(order, graph_dissection(&adj));
        // Two disconnected blocks and an isolated node.
        let (_, mut two) = hex_grid([6, 6, 6], false, |i, j, k| [i as f64, j as f64, k as f64]);
        let m = two.len();
        two.extend(two.clone().into_iter().map(|row| row.into_iter().map(|v| v + m).collect()));
        two.push(Vec::new());
        assert!(is_permutation(&graph_dissection(&two), 2 * m + 1));
        assert!(graph_dissection(&[]).is_empty());
    }

    #[test]
    fn the_graph_cuts_a_thin_curved_wall_with_less_fill() {
        let (xyz, adj) = cylinder_wall();
        let planes = fill(&adj, &nested_dissection(&xyz, &adj));
        let graph = fill(&adj, &graph_dissection(&adj));
        assert!(graph < planes, "graph {graph} planes {planes}");
        assert_eq!(fill(&adj, &fill_reducing_order(&xyz, &adj)), graph);
    }

    #[test]
    fn a_block_keeps_the_better_ordering() {
        let (xyz, adj) = hex_grid([16, 16, 16], false, |i, j, k| [i as f64, j as f64, k as f64]);
        let planes = fill(&adj, &nested_dissection(&xyz, &adj));
        let chosen = fill(&adj, &fill_reducing_order(&xyz, &adj));
        assert!(chosen <= planes, "chosen {chosen} planes {planes}");
    }

    /// Factor entries and factor time (s) of a diagonally dominant 3-unknowns-per-node matrix on
    /// the graph `adj`, eliminated in `order`.
    fn factor_3dof(adj: &[Vec<usize>], order: &[usize]) -> (usize, f64) {
        let mut place = vec![0; adj.len()];
        for (p, &v) in order.iter().enumerate() {
            place[v] = p;
        }
        let neq = 3 * adj.len();
        let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); neq];
        for (v, row) in adj.iter().enumerate() {
            for &w in row {
                for a in 0..3 {
                    for b in 0..3 {
                        let (r, c) = (3 * place[w] + b, 3 * place[v] + a);
                        if r > c {
                            continue;
                        }
                        let value = if v == w {
                            if a == b { 30.0 } else { 0.5 }
                        } else if a == b {
                            -1.0
                        } else {
                            0.01 * ((v + w) % 3) as f64
                        };
                        cols[c].push((r, value));
                    }
                }
            }
        }
        let mut ptr = vec![0];
        let (mut idx, mut val) = (Vec::new(), Vec::new());
        for mut col in cols {
            col.sort_by_key(|e| e.0);
            idx.extend(col.iter().map(|e| e.0));
            val.extend(col.iter().map(|e| e.1));
            ptr.push(idx.len());
        }
        let upper = crate::cholesky::UpperColumns { n: neq, ptr, idx, val };
        let t = std::time::Instant::now();
        let factor = crate::cholesky::Factor::new(&upper, usize::MAX).unwrap();
        (factor.entries(), t.elapsed().as_secs_f64())
    }

    /// The measurement table of the ordering choice (fill and time per ordering).
    #[test]
    #[ignore = "measurement; run with --ignored --nocapture"]
    fn measure_orderings() {
        let cases = [
            ("cylinder wall 160x2x40", cylinder_wall()),
            (
                "bellows 120x2x120",
                hex_grid([120, 2, 120], true, |i, j, k| {
                    let t = i as f64 / 120.0 * std::f64::consts::TAU;
                    let z = k as f64;
                    // Convolutions of the wall in z, the wall thickness along the normal.
                    let r = 40.0 + 12.0 * (z / 120.0 * 8.0 * std::f64::consts::TAU).sin() + j as f64 * 0.5;
                    [r * t.cos(), r * t.sin(), z * 1.0]
                }),
            ),
            ("block 24^3", hex_grid([24, 24, 24], false, |i, j, k| [i as f64, j as f64, k as f64])),
        ];
        for (name, (xyz, adj)) in cases {
            let t = std::time::Instant::now();
            let planes = nested_dissection(&xyz, &adj);
            let tp = t.elapsed();
            let t = std::time::Instant::now();
            let graph = graph_dissection(&adj);
            let tg = t.elapsed();
            let t = std::time::Instant::now();
            let (fp, fg) = (fill(&adj, &planes), fill(&adj, &graph));
            let tf = t.elapsed();
            let (ep, sp) = factor_3dof(&adj, &planes);
            let (eg, sg) = factor_3dof(&adj, &graph);
            println!(
                "| {name} | {} | {fp} / {fg} | {ep} / {eg} | {sp:.3} / {sg:.3} | {:.1} / {:.1} / {:.1} |",
                3 * adj.len(),
                tp.as_secs_f64() * 1e3,
                tg.as_secs_f64() * 1e3,
                tf.as_secs_f64() * 1e3
            );
        }
    }
}
