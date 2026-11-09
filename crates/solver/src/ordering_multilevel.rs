//! Deterministic multilevel nested dissection, independently implemented from:
//! Karypis & Kumar, SIAM J. Sci. Comput. 20(1), 359–392 (1998), sections 3–5;
//! Fiduccia & Mattheyses, DAC 1982 (locked moves and best-prefix rollback);
//! Amestoy, Davis & Duff, SIAM J. Matrix Anal. Appl. 17(4), 886–905 (1996)
//! (quotient-graph approximate degrees for bounded leaves).
//! No external partitioner is linked. Unlike randomized HEM/GGGP, visits and
//! starts are fixed; refinement uses heaps rather than claiming linear-time buckets.
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};

// Exclusive accumulated serial phase times; clocks never affect decisions.
struct Profile {
    start: Option<std::time::Instant>,
    capacity_prune: bool,
    moves: std::cell::Cell<[usize; 4]>,
    impossible_seconds: std::cell::Cell<f64>,
    seconds: std::cell::Cell<[f64; 7]>,
    calls: std::cell::Cell<[usize; 7]>,
}
struct Phase<'a> {
    profile: &'a Profile,
    index: usize,
    start: Option<std::time::Instant>,
}
impl Profile {
    fn new() -> Self {
        Self {
            start: std::env::var("FIELDCRAFT_ORDERING_PROFILE").is_ok_and(|v| v == "1").then(std::time::Instant::now),
            capacity_prune: !std::env::var("FIELDCRAFT_ORDERING_REFINEMENT").is_ok_and(|v| v == "baseline"),
            moves: [0; 4].into(),
            impossible_seconds: 0.0.into(),
            seconds: [0.0; 7].into(),
            calls: [0; 7].into(),
        }
    }
    fn phase(&self, index: usize) -> Phase<'_> {
        Phase { profile: self, index, start: self.start.map(|_| std::time::Instant::now()) }
    }
    fn move_query(&self, start: Option<std::time::Instant>, pops: usize, blocked: usize, impossible: bool) {
        if let Some(start) = start {
            let old = self.moves.get();
            self.moves.set([old[0] + 1, old[1] + usize::from(impossible), old[2] + pops, old[3] + blocked]);
            if impossible {
                self.impossible_seconds.set(self.impossible_seconds.get() + start.elapsed().as_secs_f64());
            }
        }
    }
    fn finish(&self) {
        if let Some(start) = self.start {
            use std::io::Write;
            let total = start.elapsed().as_secs_f64();
            let times = self.seconds.get();
            let [queries, impossible, pops, blocked] = self.moves.get();
            let _ = writeln!(
                std::io::stderr().lock(),
                "FIELDCRAFT_ORDER_FM_QUERY\t{queries}\t{impossible}\t{pops}\t{blocked}\t{:.9}",
                self.impossible_seconds.get()
            );
            for ((name, seconds), calls) in ["graph_build", "coarsening", "initial_partition", "refinement", "separator", "leaf", "uncoarsening"]
                .iter()
                .zip(times)
                .zip(self.calls.get())
            {
                let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_ORDER_PHASE\t{name}\t{seconds:.9}\t{calls}");
            }
            let _ =
                writeln!(std::io::stderr().lock(), "FIELDCRAFT_ORDER_PHASE\trecursion_other\t{:.9}\t1", (total - times.iter().sum::<f64>()).max(0.0));
            let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_ORDER_PHASE\ttotal\t{total:.9}\t1");
        }
    }
}
impl Drop for Phase<'_> {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            let mut times = self.profile.seconds.get();
            let mut calls = self.profile.calls.get();
            if let (Some(t), Some(c)) = (times.get_mut(self.index), calls.get_mut(self.index)) {
                *t += start.elapsed().as_secs_f64();
                *c += 1;
            }
            self.profile.seconds.set(times);
            self.profile.calls.set(calls);
        }
    }
}

const NONE: usize = usize::MAX;
const LEAF: usize = 64;
const COARSE: usize = 128;

#[derive(Clone)]
struct Graph {
    rows: Vec<Vec<(usize, i64)>>,
    weights: Vec<usize>,
    keys: Vec<usize>,
}

fn value<T: Copy + Default>(a: &[T], i: usize) -> T {
    a.get(i).copied().unwrap_or_default()
}

impl Graph {
    fn row(&self, v: usize) -> &[(usize, i64)] {
        self.rows.get(v).map(Vec::as_slice).unwrap_or_default()
    }

    fn n(&self) -> usize {
        self.rows.len()
    }

    fn induced(adj: &[Vec<usize>], part: &[usize], position: &mut [usize]) -> Self {
        for (i, &v) in part.iter().enumerate() {
            if let Some(p) = position.get_mut(v) {
                *p = i;
            }
        }
        let rows = part
            .iter()
            .map(|&v| {
                adj.get(v)
                    .into_iter()
                    .flatten()
                    .filter_map(|&u| {
                        let j = position.get(u).copied().unwrap_or(NONE);
                        (j != NONE).then_some((j, 1))
                    })
                    .collect()
            })
            .collect();
        for &v in part {
            if let Some(p) = position.get_mut(v) {
                *p = NONE;
            }
        }
        Self { rows, weights: vec![1; part.len()], keys: part.to_vec() }
    }

    fn cut(&self, side: &[usize]) -> i64 {
        self.rows
            .iter()
            .enumerate()
            .map(|(v, row)| row.iter().filter(|&&(u, _)| u > v && value(side, u) != value(side, v)).map(|&(_, w)| w).sum::<i64>())
            .sum()
    }

    fn gain(&self, v: usize, side: &[usize]) -> i64 {
        self.row(v).iter().map(|&(u, w)| if value(side, u) == value(side, v) { -w } else { w }).sum()
    }

    fn coarsen(&self, profile: &Profile) -> Option<(Self, Vec<usize>)> {
        let _phase = profile.phase(1);
        let mut visits: Vec<_> = (0..self.n()).collect();
        visits.sort_unstable_by_key(|&v| (self.row(v).len(), value(&self.keys, v)));
        let mut map = vec![NONE; self.n()];
        let (mut weights, mut keys) = (Vec::new(), Vec::new());
        let limit = (self.weights.iter().sum::<usize>() / 32).max(2);
        for v in visits {
            if map.get(v) != Some(&NONE) {
                continue;
            }
            let partner = self
                .row(v)
                .iter()
                .filter(|&&(u, _)| map.get(u) == Some(&NONE) && value(&self.weights, v) + value(&self.weights, u) <= limit)
                .max_by_key(|&&(u, w)| (w, Reverse(value(&self.keys, u))))
                .map(|&(u, _)| u);
            let k = weights.len();
            let mut weight = value(&self.weights, v);
            let mut key = value(&self.keys, v);
            if let Some(m) = map.get_mut(v) {
                *m = k;
            }
            if let Some(u) = partner {
                if let Some(m) = map.get_mut(u) {
                    *m = k;
                }
                weight += value(&self.weights, u);
                key = key.min(value(&self.keys, u));
            }
            weights.push(weight);
            keys.push(key);
        }
        if weights.len() * 10 >= self.n() * 9 {
            return None;
        }
        let mut rows = vec![Vec::new(); weights.len()];
        for (v, row) in self.rows.iter().enumerate() {
            let a = value(&map, v);
            for &(u, w) in row {
                let b = value(&map, u);
                if a != b
                    && let Some(dst) = rows.get_mut(a)
                {
                    dst.push((b, w));
                }
            }
        }
        for row in &mut rows {
            row.sort_unstable_by_key(|e| e.0);
            let mut merged: Vec<(usize, i64)> = Vec::new();
            for &(u, w) in row.iter() {
                if let Some(last) = merged.last_mut()
                    && last.0 == u
                {
                    last.1 += w;
                    continue;
                }
                merged.push((u, w));
            }
            *row = merged;
        }
        Some((Self { rows, weights, keys }, map))
    }
}

type Heap = BinaryHeap<(i64, Reverse<usize>, usize, usize)>;

fn best_move(
    heap: &mut Heap,
    locked: &[bool],
    versions: &[usize],
    weights: &[usize],
    room: usize,
    minimum_weight: usize,
    profile: &Profile,
) -> Option<(i64, usize)> {
    let start = profile.start.map(|_| std::time::Instant::now());
    let impossible = room < minimum_weight;
    // No unlocked vertex can fit if even the global minimum weight exceeds room.
    // The reference heap scan would return None after discarding stale entries and
    // restoring every live entry. Retaining those stale entries changes no choice.
    if profile.capacity_prune && impossible {
        profile.move_query(start, 0, 0, true);
        return None;
    }
    let mut pops = 0;
    let mut rejected = 0;
    let mut blocked = Vec::new();
    let mut best = None;
    while let Some(entry @ (gain, _, v, revision)) = heap.pop() {
        pops += 1;
        if value(locked, v) || value(versions, v) != revision {
            continue;
        }
        if value(weights, v) > room {
            rejected += 1;
            blocked.push(entry);
            continue;
        }
        best = Some((gain, v));
        // Peek semantics: the selected vertex is locked or its version changes after selection.
        blocked.push(entry);
        break;
    }
    heap.extend(blocked);
    profile.move_query(start, pops, rejected, impossible);
    best
}

/// Weighted graph FM: locked single moves, bounded non-improving search, rollback
/// to the best balanced prefix. Every neighbor-gain update follows sorted adjacency.
fn refine(g: &Graph, side: &mut [usize], profile: &Profile) {
    let _phase = profile.phase(3);
    let minimum_weight = g.weights.iter().copied().min().unwrap_or(1);
    let total: usize = g.weights.iter().sum();
    let cap = (total * 53 / 100 + g.weights.iter().copied().max().unwrap_or(1)).min(total.saturating_sub(1));
    for _ in 0..4 {
        let mut weights = [0usize; 2];
        for (v, &s) in side.iter().enumerate() {
            if let Some(w) = weights.get_mut(s) {
                *w += value(&g.weights, v);
            }
        }
        let mut gains: Vec<_> = (0..g.n()).map(|v| g.gain(v, side)).collect();
        let mut versions = vec![0; g.n()];
        let mut locked = vec![false; g.n()];
        let mut heaps: [Heap; 2] = std::array::from_fn(|_| Heap::new());
        for (v, &gain) in gains.iter().enumerate() {
            if let Some(h) = heaps.get_mut(value(side, v)) {
                h.push((gain, Reverse(value(&g.keys, v)), v, 0));
            }
        }
        let (mut cumulative, mut best_gain, mut prefix) = (0i64, 0i64, 0usize);
        let mut moves = Vec::new();
        let stall_limit = (g.n() / 100 + 32).min(512);
        for _ in 0..g.n() {
            let candidates: Vec<_> = heaps
                .iter_mut()
                .enumerate()
                .filter_map(|(s, h)| {
                    best_move(h, &locked, &versions, &g.weights, cap.saturating_sub(value(&weights, 1 - s)), minimum_weight, profile)
                        .map(|(gain, v)| (gain, v, s))
                })
                .collect();
            let Some(&(gain, v, old)) = candidates.iter().max_by_key(|&&(gain, v, _)| (gain, Reverse(value(&g.keys, v)))) else { break };
            let weight = value(&g.weights, v);
            if let Some(w) = weights.get_mut(old) {
                *w -= weight;
            }
            if let Some(w) = weights.get_mut(1 - old) {
                *w += weight;
            }
            if let Some(s) = side.get_mut(v) {
                *s = 1 - old;
            }
            if let Some(l) = locked.get_mut(v) {
                *l = true;
            }
            moves.push((v, old));
            cumulative += gain;
            if cumulative > best_gain && weights.iter().all(|&w| w <= cap) {
                best_gain = cumulative;
                prefix = moves.len();
            }
            for &(u, edge) in g.row(v) {
                if value(&locked, u) {
                    continue;
                }
                let delta = if value(side, u) == old { 2 * edge } else { -2 * edge };
                if let Some(x) = gains.get_mut(u) {
                    *x += delta;
                }
                if let Some(x) = versions.get_mut(u) {
                    *x += 1;
                }
                if let Some(h) = heaps.get_mut(value(side, u)) {
                    h.push((value(&gains, u), Reverse(value(&g.keys, u)), u, value(&versions, u)));
                }
            }
            if moves.len() - prefix > stall_limit {
                break;
            }
        }
        for &(v, old) in moves.get(prefix..).unwrap_or_default().iter().rev() {
            if let Some(s) = side.get_mut(v) {
                *s = old;
            }
        }
        if prefix == 0 {
            break;
        }
    }
}

fn grow(g: &Graph, seed: usize, profile: &Profile) -> Vec<usize> {
    let phase = profile.phase(2);
    let mut side = vec![1usize; g.n()];
    let target = g.weights.iter().sum::<usize>() / 2;
    let mut gains: Vec<_> = (0..g.n()).map(|v| g.gain(v, &side)).collect();
    let mut versions = vec![0usize; g.n()];
    let mut heap = BinaryHeap::new();
    for (v, &gain) in gains.iter().enumerate() {
        heap.push((false, gain, Reverse(value(&g.keys, v)), v, 0));
    }
    let mut held = 0;
    let mut next = Some(seed);
    for _ in 0..g.n() {
        let Some(v) = next else { break };
        if let Some(s) = side.get_mut(v) {
            *s = 0;
        }
        held += value(&g.weights, v);
        if held >= target {
            break;
        }
        for &(u, w) in g.row(v) {
            if value(&side, u) == 0 {
                continue;
            }
            if let Some(gain) = gains.get_mut(u) {
                *gain += 2 * w;
            }
            if let Some(revision) = versions.get_mut(u) {
                *revision += 1;
            }
            heap.push((true, value(&gains, u), Reverse(value(&g.keys, u)), u, value(&versions, u)));
        }
        next = None;
        while let Some((_, _, _, u, revision)) = heap.pop() {
            if value(&side, u) == 1 && value(&versions, u) == revision {
                next = Some(u);
                break;
            }
        }
    }
    drop(phase);
    refine(g, &mut side, profile);
    side
}

fn partition(g: &Graph, profile: &Profile) -> (Vec<usize>, usize) {
    let mut graphs = vec![g.clone()];
    let mut maps = Vec::new();
    for _ in 0..32 {
        let Some(last) = graphs.last() else { break };
        if last.n() <= COARSE {
            break;
        }
        let Some((coarse, map)) = last.coarsen(profile) else { break };
        maps.push(map);
        graphs.push(coarse);
    }
    let Some(last) = graphs.last() else { return (Vec::new(), 0) };
    // Fixed starts distributed over the stable coarse ids. No RNG or scheduling dependence.
    let mut seeds: Vec<_> = (0..8).map(|i| i * last.n() / 8).filter(|&v| v < last.n()).collect();
    seeds.sort_unstable();
    seeds.dedup();
    let mut side = seeds.into_iter().map(|v| grow(last, v, profile)).min_by_key(|s| (last.cut(s), s.clone())).unwrap_or_default();
    for level in (0..maps.len()).rev() {
        if let (Some(map), Some(fine)) = (maps.get(level), graphs.get(level)) {
            let phase = profile.phase(6);
            side = map.iter().map(|&v| value(&side, v)).collect();
            drop(phase);
            refine(fine, &mut side, profile);
        }
    }
    (side, maps.len())
}

/// Minimum cardinality vertex cover of the bipartite crossing-edge graph.
/// Deterministic augmenting paths, followed by alternating reachability (Kőnig).
fn cover(g: &Graph, side: &[usize], profile: &Profile) -> Vec<bool> {
    let _phase = profile.phase(4);
    let mut mate = vec![NONE; g.n()];
    let mut parent = vec![NONE; g.n()];
    let mut seen = vec![0usize; g.n()];
    let mut queue = VecDeque::new();
    for root in 0..g.n() {
        if value(side, root) != 0 || mate.get(root) != Some(&NONE) {
            continue;
        }
        let token = root + 1;
        queue.clear();
        queue.push_back(root);
        if let Some(s) = seen.get_mut(root) {
            *s = token;
        }
        let mut end = None;
        while let Some(v) = queue.pop_front() {
            for &(u, _) in g.row(v) {
                if value(side, u) != 1 || value(&seen, u) == token {
                    continue;
                }
                if let Some(s) = seen.get_mut(u) {
                    *s = token;
                }
                if let Some(p) = parent.get_mut(u) {
                    *p = v;
                }
                let other = mate.get(u).copied().unwrap_or(NONE);
                if other == NONE {
                    end = Some(u);
                    break;
                }
                if value(&seen, other) != token {
                    if let Some(s) = seen.get_mut(other) {
                        *s = token;
                    }
                    queue.push_back(other);
                }
            }
            if end.is_some() {
                break;
            }
        }
        while let Some(u) = end {
            let v = parent.get(u).copied().unwrap_or(NONE);
            let old = mate.get(v).copied().unwrap_or(NONE);
            if let Some(m) = mate.get_mut(u) {
                *m = v;
            }
            if let Some(m) = mate.get_mut(v) {
                *m = u;
            }
            end = (old != NONE).then_some(old);
        }
    }
    let mut reachable = vec![false; g.n()];
    queue.clear();
    for v in 0..g.n() {
        if value(side, v) == 0 && mate.get(v) == Some(&NONE) {
            queue.push_back(v);
            if let Some(r) = reachable.get_mut(v) {
                *r = true;
            }
        }
    }
    while let Some(v) = queue.pop_front() {
        for &(u, _) in g.row(v) {
            if value(side, u) != 1 || mate.get(v) == Some(&u) || value(&reachable, u) {
                continue;
            }
            if let Some(r) = reachable.get_mut(u) {
                *r = true;
            }
            let other = mate.get(u).copied().unwrap_or(NONE);
            if other != NONE && !value(&reachable, other) {
                if let Some(r) = reachable.get_mut(other) {
                    *r = true;
                }
                queue.push_back(other);
            }
        }
    }
    (0..g.n()).map(|v| if value(side, v) == 0 { !value(&reachable, v) } else { value(&reachable, v) }).collect()
}

/// Quotient-graph approximate minimum degree for bounded leaves. Element degree
/// sums are upper bounds (overlap is allowed); no claim of exact minimum degree or
/// of reproducing another AMD implementation's ordering.
fn amd(g: &Graph, profile: &Profile) -> Vec<usize> {
    let _phase = profile.phase(5);
    let mut direct: Vec<Vec<usize>> = g.rows.iter().map(|r| r.iter().map(|e| e.0).collect()).collect();
    let mut incident = vec![Vec::<usize>::new(); g.n()];
    let mut elements: Vec<Vec<usize>> = Vec::new();
    let mut active = vec![true; g.n()];
    let mut order = Vec::with_capacity(g.n());
    for remaining in (1..=g.n()).rev() {
        let degree = |v: usize| {
            let direct_count = direct.get(v).into_iter().flatten().filter(|&&u| value(&active, u)).count();
            let element_count: usize = incident
                .get(v)
                .into_iter()
                .flatten()
                .map(|&e| elements.get(e).into_iter().flatten().filter(|&&u| value(&active, u)).count().saturating_sub(1))
                .sum();
            (direct_count + element_count).min(remaining - 1)
        };
        let Some(v) = (0..g.n()).filter(|&v| value(&active, v)).min_by_key(|&v| (degree(v), value(&g.keys, v))) else { break };
        let old = incident.get(v).cloned().unwrap_or_default();
        let mut neighbors = direct.get(v).cloned().unwrap_or_default();
        for &e in &old {
            neighbors.extend(elements.get(e).into_iter().flatten().copied());
        }
        neighbors.retain(|&u| u != v && value(&active, u));
        neighbors.sort_unstable();
        neighbors.dedup();
        if let Some(a) = active.get_mut(v) {
            *a = false;
        }
        for &e in &old {
            if let Some(row) = elements.get_mut(e) {
                row.clear();
            }
        }
        let e = elements.len();
        elements.push(neighbors.clone());
        for &u in &neighbors {
            if let Some(row) = incident.get_mut(u) {
                row.retain(|&i| elements.get(i).is_some_and(|x| !x.is_empty()));
                row.push(e);
            }
            if let Some(row) = direct.get_mut(u) {
                row.retain(|&i| value(&active, i) && neighbors.binary_search(&i).is_err());
            }
        }
        order.push(value(&g.keys, v));
    }
    order
}

pub(crate) fn order(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let profile = Profile::new();
    let build = profile.phase(0);
    let n = adjacency.len();
    // Canonical undirected graph; duplicate, self and out-of-range entries cannot
    // affect traversal or create input-dependent unchecked indexing.
    let mut adj = vec![Vec::new(); n];
    for (v, row) in adjacency.iter().enumerate() {
        for &u in row {
            if u < n && u != v {
                if let Some(r) = adj.get_mut(v) {
                    r.push(u);
                }
                if let Some(r) = adj.get_mut(u) {
                    r.push(v);
                }
            }
        }
    }
    for row in &mut adj {
        row.sort_unstable();
        row.dedup();
    }
    drop(build);
    let mut position = vec![NONE; n];
    enum Task {
        Split(Vec<usize>, usize),
        Emit(Vec<usize>),
    }
    let mut stack = vec![Task::Split((0..n).collect(), 0)];
    let mut result = Vec::with_capacity(n);
    while let Some(task) = stack.pop() {
        let (part, depth) = match task {
            Task::Emit(p) => {
                result.extend(p);
                continue;
            }
            Task::Split(p, d) => (p, d),
        };
        if part.is_empty() {
            continue;
        }
        let build = profile.phase(0);
        let g = Graph::induced(&adj, &part, &mut position);
        drop(build);
        if part.len() <= LEAF {
            result.extend(amd(&g, &profile));
            continue;
        }
        // Separate connected components before coarsening, avoiding artificial cuts.
        let mut seen = vec![false; part.len()];
        let mut components = Vec::new();
        for root in 0..part.len() {
            if value(&seen, root) {
                continue;
            }
            let mut q = vec![root];
            if let Some(s) = seen.get_mut(root) {
                *s = true;
            }
            let mut head = 0;
            while let Some(&v) = q.get(head) {
                head += 1;
                for &(u, _) in g.row(v) {
                    if !value(&seen, u) {
                        if let Some(s) = seen.get_mut(u) {
                            *s = true;
                        }
                        q.push(u);
                    }
                }
            }
            let mut component: Vec<_> = q.iter().map(|&v| value(&g.keys, v)).collect();
            component.sort_unstable();
            components.push(component);
        }
        if components.len() > 1 {
            for c in components.into_iter().rev() {
                stack.push(Task::Split(c, depth));
            }
            continue;
        }
        let (side, levels) = partition(&g, &profile);
        let separator = cover(&g, &side, &profile);
        let (mut left, mut right, mut sep) = (Vec::new(), Vec::new(), Vec::new());
        for (v, &key) in g.keys.iter().enumerate() {
            if value(&separator, v) {
                sep.push(key);
            } else if value(&side, v) == 0 {
                left.push(key);
            } else {
                right.push(key);
            }
        }
        if std::env::var("FIELDCRAFT_ORDERING_TRACE").is_ok_and(|v| v == "1") {
            use std::io::Write;
            let _ = writeln!(
                std::io::stderr().lock(),
                "FIELDCRAFT_ORDER_SEPARATOR\t{depth}\t{}\t{}\t{}\t{}\t{levels}\t{}",
                part.len(),
                left.len(),
                right.len(),
                sep.len(),
                g.cut(&side)
            );
        }
        // Bounded fallback for dense/unsplittable parts, with a proper permutation.
        if left.is_empty() || right.is_empty() || depth >= 64 {
            result.extend(part);
            continue;
        }
        stack.push(Task::Emit(sep));
        stack.push(Task::Split(right, depth + 1));
        stack.push(Task::Split(left, depth + 1));
    }
    profile.finish();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cover_is_minimum_on_every_three_by_three_crossing_graph() {
        for bits in 0..512 {
            let mut g = Graph { rows: vec![vec![]; 6], weights: vec![1; 6], keys: (0..6).collect() };
            for a in 0..3 {
                for b in 0..3 {
                    if bits & (1 << (3 * a + b)) != 0 {
                        g.rows[a].push((3 + b, 1));
                        g.rows[3 + b].push((a, 1));
                    }
                }
            }
            let side = [0, 0, 0, 1, 1, 1];
            let got = cover(&g, &side, &Profile::new());
            let valid = |mask: usize| g.rows.iter().enumerate().all(|(v, row)| row.iter().all(|&(u, _)| mask & (1 << v | 1 << u) != 0));
            let want = (0usize..64).filter(|&m| valid(m)).map(usize::count_ones).min().unwrap();
            let mask = got.iter().enumerate().fold(0, |m, (v, &yes)| if yes { m | 1 << v } else { m });
            assert!(valid(mask));
            assert_eq!(got.iter().filter(|&&v| v).count(), want as usize);
        }
    }

    #[test]
    fn capacity_pruning_preserves_weighted_heap_choices_and_stale_entries() {
        for weights in [vec![1, 1, 1, 1], vec![2, 3, 5, 7]] {
            let minimum = *weights.iter().min().unwrap();
            let mut reference = Profile::new();
            reference.capacity_prune = false;
            let mut candidate = Profile::new();
            candidate.capacity_prune = true;
            let mut a: Heap = (0..4).flat_map(|v| [(10 - v as i64, Reverse(v), v, 0), (20 - v as i64, Reverse(v), v, 1)]).collect();
            let mut b = a.clone();
            let mut locked = vec![false; 4];
            let versions = vec![1; 4];
            for room in [0, 1, 0, 2, 0, 3, 5, 0, 7, 1, 9] {
                let x = best_move(&mut a, &locked, &versions, &weights, room, minimum, &reference);
                let y = best_move(&mut b, &locked, &versions, &weights, room, minimum, &candidate);
                assert_eq!(x, y);
                if let Some((_, v)) = x {
                    locked[v] = true;
                }
            }
        }
    }

    #[test]
    fn coarsening_preserves_weights_and_partition_cut() {
        let n = 64;
        let g = Graph {
            rows: (0..n).map(|v| [(v + n - 1) % n, (v + 1) % n].into_iter().map(|u| (u, 1)).collect()).collect(),
            weights: vec![1; n],
            keys: (0..n).collect(),
        };
        let (coarse, map) = g.coarsen(&Profile::new()).unwrap();
        assert_eq!(coarse.weights.iter().sum::<usize>(), n);
        for stride in 1..5 {
            let side: Vec<_> = (0..coarse.n()).map(|v| (v / stride) % 2).collect();
            let fine: Vec<_> = map.iter().map(|&v| side[v]).collect();
            assert_eq!(g.cut(&fine), coarse.cut(&side));
        }
    }

    #[test]
    fn fm_keeps_the_best_prefix_and_balance() {
        let n = 80;
        let g = Graph {
            rows: (0..n).map(|v| [(v + n - 1) % n, (v + 1) % n].into_iter().map(|u| (u, 1)).collect()).collect(),
            weights: vec![1; n],
            keys: (0..n).collect(),
        };
        let mut side: Vec<_> = (0..n).map(|v| v % 2).collect();
        let before = g.cut(&side);
        refine(&g, &mut side, &Profile::new());
        assert!(g.cut(&side) < before);
        assert!((36..=44).contains(&side.iter().filter(|&&s| s == 0).count()));
    }

    #[test]
    fn disconnected_dense_star_and_invalid_entries_are_bounded() {
        for mut adj in [
            vec![],
            vec![vec![]; 200],
            (0..90).map(|v| (0..90).filter(|&u| u != v).collect()).collect(),
            (0..200).map(|v| if v == 0 { (1..200).collect() } else { vec![0] }).collect(),
        ] {
            if let Some(row) = adj.first_mut() {
                row.extend([NONE, 0, 0]);
            }
            let got = order(&adj);
            let mut sorted = got.clone();
            sorted.sort_unstable();
            assert_eq!(sorted, (0..adj.len()).collect::<Vec<_>>());
            for row in &mut adj {
                row.reverse();
            }
            assert_eq!(got, order(&adj));
        }
    }
}
