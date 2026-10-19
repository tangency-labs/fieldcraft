//! The block topology of the projection method: hexahedral blocks over shared vertices, built
//! from a carved grid, from the round-feature templates (butterfly, seven-block sphere) and by the
//! O-grid and extrude operations.
//!
//! A block's corners are in hex8 order (1-2-3-4 counter-clockwise seen from 5-6-7-8's side) and
//! its parameter directions run 1→2, 1→4 and 1→5. Faces are numbered as hex8 faces: 1 = 1-2-3-4,
//! 2 = 5-6-7-8, 3 = 1-2-6-5, 4 = 2-3-7-6, 5 = 3-4-8-7, 6 = 4-1-5-8. Positions are only a first
//! guess: the projection moves every boundary vertex onto the geometry.

use std::collections::BTreeMap;

use crate::vec3::{P, add, det, dist, scale, sub};
use crate::{Error, Result};

/// Vertices this close, relative to the template's size, are one.
const MERGE: f64 = 1e-9; // unit-factors: relative to the template's size

/// Most blocks one topology holds.
pub const MAX_BLOCKS: usize = 20_000;

/// Corner `c` of a block in parameter space (hex8 order).
pub(crate) const CORNERS: [[usize; 3]; 8] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];

/// The direction and side of each hex8 face (1-based face `f` at index `f - 1`).
pub(crate) const FACE_SIDES: [(usize, usize); 6] = [(2, 0), (2, 1), (1, 0), (0, 1), (1, 1), (0, 0)];

/// The corner at parameter `c`.
pub(crate) fn corner_at(c: [usize; 3]) -> usize {
    CORNERS.iter().position(|x| *x == c).unwrap_or(0)
}

/// Hex8 face number (1-based) of direction `d`, side `s`.
pub(crate) fn face_number(d: usize, s: usize) -> u8 {
    FACE_SIDES.iter().position(|x| *x == (d, s)).map_or(1, |i| i as u8 + 1)
}

/// The four edges of a block along direction `d`, each from its side-0 corner to its side-1 corner.
fn edges_along(block: &[u32; 8], d: usize) -> [[u32; 2]; 4] {
    let mut out = [[0; 2]; 4];
    for (slot, (k, bits)) in out.iter_mut().zip(CORNERS.iter().enumerate().filter(|(_, bits)| bits[d] == 0)) {
        let mut e = *bits;
        e[d] = 1;
        *slot = [block[k], block[corner_at(e)]];
    }
    out
}

/// Blocks over shared vertices.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Topology {
    pub vertices: Vec<P>,
    pub blocks: Vec<[u32; 8]>,
}

/// The trilinear map of a block's corners at parameter `u`.
pub(crate) fn trilinear(c: &[P; 8], u: [f64; 3]) -> P {
    let mut x = [0.0; 3];
    for (k, bits) in CORNERS.iter().enumerate() {
        let w: f64 = (0..3).map(|d| if bits[d] == 1 { u[d] } else { 1.0 - u[d] }).product();
        x = add(x, scale(c[k], w));
    }
    x
}

/// The smallest corner Jacobian determinant of a trilinear block (positive for a valid one).
pub(crate) fn min_corner_det(c: &[P; 8]) -> f64 {
    let mut m = f64::INFINITY;
    for (k, bits) in CORNERS.iter().enumerate() {
        // The three edges leaving corner k along increasing parameters, signed to keep the
        // right-handed order.
        let along = |d: usize| {
            let mut b = *bits;
            b[d] = 1 - b[d];
            let e = sub(c[corner_at(b)], c[k]);
            if bits[d] == 1 { scale(e, -1.0) } else { e }
        };
        m = m.min(det(along(0), along(1), along(2)));
    }
    m
}

impl Topology {
    /// The corner positions of block `b`.
    pub fn corners(&self, b: usize) -> Result<[P; 8]> {
        let block = self.blocks.get(b).ok_or_else(|| Error::Invalid(format!("there is no block {b}")))?;
        let mut out = [[0.0; 3]; 8];
        for (o, &v) in out.iter_mut().zip(block) {
            *o = *self.vertices.get(v as usize).ok_or_else(|| Error::Invalid(format!("block {b} uses vertex {v}, which does not exist")))?;
        }
        Ok(out)
    }

    /// Every block uses existing, distinct vertices, is right-handed and caps are kept.
    pub fn check(&self) -> Result<()> {
        if self.blocks.is_empty() {
            return Err(Error::Invalid("the topology has no blocks".into()));
        }
        if self.blocks.len() > MAX_BLOCKS {
            return Err(Error::Limit(format!("{} blocks; at most {MAX_BLOCKS}", self.blocks.len())));
        }
        if let Some(v) = self.vertices.iter().position(|v| !v.iter().all(|x| x.is_finite() && x.abs() <= fieldcraft_model::MAX_COORD)) {
            return Err(Error::Invalid(format!("vertex {v} is not a finite point")));
        }
        for (b, block) in self.blocks.iter().enumerate() {
            let mut sorted = *block;
            sorted.sort_unstable();
            if sorted.windows(2).any(|w| w[0] == w[1]) {
                return Err(Error::Invalid(format!("block {b} uses one vertex twice")));
            }
            if !(min_corner_det(&self.corners(b)?) > 0.0) {
                return Err(Error::Invalid(format!(
                    "block {b} is inside out or folded at a corner: give its corners in hex8 order (1-2-3-4 counter-clockwise seen from 5-6-7-8)"
                )));
            }
        }
        Ok(())
    }

    /// The bounding-box diagonal of the vertices.
    pub fn diagonal(&self) -> f64 {
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for v in &self.vertices {
            for d in 0..3 {
                lo[d] = lo[d].min(v[d]);
                hi[d] = hi[d].max(v[d]);
            }
        }
        if self.vertices.is_empty() { 0.0 } else { dist(lo, hi) }
    }

    /// A grid of blocks between the planes `planes[d]` (increasing coordinates along X, Y and Z),
    /// block `(i, j, k)` at index `(k·nj + j)·ni + i`.
    pub fn grid(planes: [&[f64]; 3]) -> Result<Topology> {
        for (d, p) in planes.iter().enumerate() {
            if p.len() < 2 || !p.iter().all(|x| x.is_finite()) || p.windows(2).any(|w| !(w[1] > w[0])) {
                return Err(Error::Invalid(format!("the grid planes along {} must be at least two increasing numbers", ["X", "Y", "Z"][d])));
            }
        }
        let n = planes.map(|p| p.len() - 1);
        if n.iter().try_fold(1usize, |a, &b| a.checked_mul(b)).is_none_or(|c| c > MAX_BLOCKS) {
            return Err(Error::Limit(format!("a grid of {n:?} blocks; at most {MAX_BLOCKS}")));
        }
        let vid = |i: usize, j: usize, k: usize| ((k * (n[1] + 1) + j) * (n[0] + 1) + i) as u32;
        let mut vertices = Vec::new();
        for z in planes[2] {
            for y in planes[1] {
                for x in planes[0] {
                    vertices.push([*x, *y, *z]);
                }
            }
        }
        let mut blocks = Vec::new();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    blocks.push(CORNERS.map(|c| vid(i + c[0], j + c[1], k + c[2])));
                }
            }
        }
        Ok(Topology { vertices, blocks })
    }

    /// Remove the blocks of a grid of `n` blocks per direction whose index lies in every range
    /// of `ranges` (inclusive, per direction), keeping the others in order; vertices no block
    /// uses are dropped.
    pub fn carve(&mut self, n: [usize; 3], ranges: &[[[usize; 2]; 3]]) -> Result<()> {
        if self.blocks.len() != n.iter().product::<usize>() {
            return Err(Error::Invalid("carving needs the grid's own block numbering".into()));
        }
        let mut keep = Vec::new();
        for (b, block) in self.blocks.iter().enumerate() {
            let at = [b % n[0], (b / n[0]) % n[1], b / (n[0] * n[1])];
            let gone = ranges.iter().any(|r| (0..3).all(|d| r[d][0] <= at[d] && at[d] <= r[d][1]));
            if !gone {
                keep.push(*block);
            }
        }
        self.blocks = keep;
        self.compact();
        if self.blocks.is_empty() { Err(Error::Invalid("carving removed every block".into())) } else { Ok(()) }
    }

    /// Drop unused vertices and renumber.
    fn compact(&mut self) {
        let mut map: BTreeMap<u32, u32> = BTreeMap::new();
        let mut vertices = Vec::new();
        for block in &mut self.blocks {
            for v in block.iter_mut() {
                let next = map.len() as u32;
                let new = *map.entry(*v).or_insert_with(|| {
                    vertices.push(self.vertices.get(*v as usize).copied().unwrap_or([0.0; 3]));
                    next
                });
                *v = new;
            }
        }
        self.vertices = vertices;
    }

    /// Join vertices closer than `tol` (the first one stays).
    pub fn merge(&mut self, tol: f64) {
        let mut map: Vec<u32> = (0..self.vertices.len() as u32).collect();
        // Sorted by X, so only a window is compared.
        let mut order: Vec<usize> = (0..self.vertices.len()).collect();
        order.sort_by(|&a, &b| self.vertices[a][0].total_cmp(&self.vertices[b][0]).then(a.cmp(&b)));
        for (oi, &a) in order.iter().enumerate() {
            if map[a] as usize != a {
                continue;
            }
            for &b in order.iter().skip(oi + 1) {
                if self.vertices[b][0] - self.vertices[a][0] > tol {
                    break;
                }
                if map[b] as usize == b && dist(self.vertices[a], self.vertices[b]) <= tol {
                    map[b] = a as u32;
                }
            }
        }
        for block in &mut self.blocks {
            for v in block.iter_mut() {
                *v = map.get(*v as usize).copied().unwrap_or(*v);
            }
        }
        self.compact();
    }

    /// Replace block `b` by an O-grid: an inner block, scaled by `inset[d]` (0 to 1) about the
    /// block's centre along each of its directions (from the open face, along a direction with one
    /// open face), and one block between each face of the inner
    /// block and the block's own face, except across the `open` faces (hex8 face numbers), where
    /// the inner block reaches the face. Faces of `b` stay faces of the result, so neighbours
    /// still fit, except the open ones, which become the inner block's face ringed by the side
    /// blocks' faces (the butterfly pattern on that face). The inner block takes index `b`; the
    /// others are appended. Returns, per open face, the faces now on it: the inner block's first,
    /// then the ring's.
    pub fn o_grid(&mut self, b: usize, open: &[u8], inset: [f64; 3]) -> Result<Vec<Vec<(usize, u8)>>> {
        if !inset.iter().all(|i| *i > 0.0 && *i < 1.0) {
            return Err(Error::Invalid(format!("the O-grid inset must be between 0 and 1, got {inset:?}")));
        }
        if let Some(f) = open.iter().find(|f| !(1..=6).contains(*f)) {
            return Err(Error::Invalid(format!("face {f} is not a hex8 face (1 to 6)")));
        }
        let c = self.corners(b)?;
        let outer = self.blocks.get(b).copied().ok_or_else(|| Error::Invalid(format!("there is no block {b}")))?;
        let is_open = |d: usize, s: usize| open.contains(&face_number(d, s));
        // Each corner's inner copy: on an open face it stays on the face.
        let mut inner = [0u32; 8];
        for (k, bits) in CORNERS.iter().enumerate() {
            // Along a direction with one open face the inner block reaches `inset` of the way
            // across from it; with none it is centred.
            let u = [0, 1, 2].map(|d| {
                let b = bits[d] as f64;
                match (is_open(d, 0), is_open(d, 1)) {
                    (true, true) => b,
                    (true, false) => {
                        if bits[d] == 0 {
                            0.0
                        } else {
                            inset[d]
                        }
                    }
                    (false, true) => {
                        if bits[d] == 1 {
                            1.0
                        } else {
                            1.0 - inset[d]
                        }
                    }
                    (false, false) => 0.5 + inset[d] * (b - 0.5),
                }
            });
            inner[k] = self.vertices.len() as u32;
            self.vertices.push(trilinear(&c, u));
        }
        if let Some(slot) = self.blocks.get_mut(b) {
            *slot = inner;
        }
        let mut on_open: Vec<Vec<(usize, u8)>> = open.iter().map(|&f| vec![(b, f)]).collect();
        for (d, s) in FACE_SIDES {
            if is_open(d, s) {
                continue;
            }
            // Local directions: the face's two directions, then from the outer face inwards.
            let (e1, e2) = match d {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let pick = |a: usize, bb: usize, w: usize| {
                let mut bits = [0; 3];
                bits[d] = s;
                bits[e1] = a;
                bits[e2] = bb;
                let k = corner_at(bits);
                if w == 0 { outer[k] } else { inner[k] }
            };
            let mut block = CORNERS.map(|q| pick(q[0], q[1], q[2]));
            let pos = |blk: &[u32; 8]| blk.map(|v| self.vertices.get(v as usize).copied().unwrap_or([0.0; 3]));
            let flipped = min_corner_det(&pos(&block)) < 0.0;
            if flipped {
                block = CORNERS.map(|q| pick(q[1], q[0], q[2]));
            }
            let index = self.blocks.len();
            self.blocks.push(block);
            // This block's face on each open face of the old block.
            for (list, &f) in on_open.iter_mut().zip(open) {
                let Some(&(od, os)) = FACE_SIDES.get(usize::from(f) - 1) else { continue };
                let local = if od == e1 {
                    0
                } else if od == e2 {
                    1
                } else {
                    continue;
                };
                let local = if flipped { 1 - local } else { local };
                list.push((index, face_number(local, os)));
            }
        }
        Ok(on_open)
    }

    /// The hex8 face of block `b` whose outward direction is nearest `dir`.
    pub fn face_towards(&self, b: usize, dir: P) -> Result<u8> {
        let c = self.corners(b)?;
        let centre = scale(c.iter().fold([0.0; 3], |a, p| add(a, *p)), 0.125);
        let mut best = (f64::NEG_INFINITY, 1u8);
        for (fi, (d, s)) in FACE_SIDES.iter().enumerate() {
            let fc = CORNERS.iter().zip(&c).filter(|(bits, _)| bits[*d] == *s).fold([0.0; 3], |a, (_, p)| add(a, *p));
            let out = sub(scale(fc, 0.25), centre);
            let l = crate::vec3::norm(out);
            let score = if l > 0.0 { crate::vec3::dot(out, dir) / l } else { f64::NEG_INFINITY };
            if score > best.0 {
                best = (score, fi as u8 + 1);
            }
        }
        Ok(best.1)
    }

    /// The direction of block `b` (0 to 2) most nearly along `dir`, and whether it runs against it.
    pub fn direction_along(&self, b: usize, dir: P) -> Result<(usize, bool)> {
        let c = self.corners(b)?;
        let mut best = (f64::NEG_INFINITY, 0, false);
        for d in 0..3 {
            let mut v = [0.0; 3];
            for [p, q] in edges_along(&[0, 1, 2, 3, 4, 5, 6, 7], d) {
                v = add(v, sub(c[q as usize], c[p as usize]));
            }
            let l = crate::vec3::norm(v);
            let cos = if l > 0.0 { crate::vec3::dot(v, dir) / l } else { 0.0 };
            if cos.abs() > best.0 {
                best = (cos.abs(), d, cos < 0.0);
            }
        }
        Ok((best.1, best.2))
    }

    /// The corners of a block face in cyclic order.
    fn face_quad(&self, b: usize, f: u8) -> Result<[u32; 4]> {
        let (d, s) = FACE_SIDES.get(usize::from(f).wrapping_sub(1)).copied().ok_or_else(|| Error::Invalid(format!("face {f} is not a hex8 face")))?;
        let block = self.blocks.get(b).copied().ok_or_else(|| Error::Invalid(format!("there is no block {b}")))?;
        let (e1, e2) = match d {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        Ok([[0, 0], [1, 0], [1, 1], [0, 1]].map(|[a, bb]| {
            let mut bits = [0; 3];
            bits[d] = s;
            bits[e1] = a;
            bits[e2] = bb;
            block[corner_at(bits)]
        }))
    }

    /// Fill the space between two block faces with `layers` blocks, each face's corner matched
    /// to the nearest of the other's (over the eight ways of laying one quad on the other).
    /// Returns the new blocks.
    pub fn connect(&mut self, from: (usize, u8), to: (usize, u8), layers: usize) -> Result<Vec<usize>> {
        if !(1..=1000).contains(&layers) {
            return Err(Error::Invalid(format!("a connection has 1 to 1000 layers, got {layers}")));
        }
        let a = self.face_quad(from.0, from.1)?;
        let b = self.face_quad(to.0, to.1)?;
        let verts = self.vertices.clone();
        let at = |v: u32| verts.get(v as usize).copied().unwrap_or([0.0; 3]);
        let mut best: Option<(f64, [u32; 4])> = None;
        for rev in [false, true] {
            for shift in 0..4 {
                let m: [u32; 4] = [0, 1, 2, 3].map(|i| {
                    let k = if rev { (4 + shift - i) % 4 } else { (shift + i) % 4 };
                    b[k]
                });
                let cost: f64 = (0..4).map(|i| dist(at(a[i]), at(m[i]))).sum();
                if best.is_none_or(|x| cost < x.0) {
                    best = Some((cost, m));
                }
            }
        }
        let (_, m) = best.ok_or_else(|| Error::Invalid("no way to match the faces".into()))?;
        // Intermediate layers of vertices, then the blocks.
        let mut rows: Vec<[u32; 4]> = vec![a];
        for l in 1..layers {
            let t = l as f64 / layers as f64;
            let pts = [0, 1, 2, 3].map(|i| crate::vec3::lerp(at(a[i]), at(m[i]), t));
            let row = pts.map(|p| {
                self.vertices.push(p);
                self.vertices.len() as u32 - 1
            });
            rows.push(row);
        }
        rows.push(m);
        let mut made = Vec::new();
        for w in rows.windows(2) {
            let (q, t) = (w[0], w[1]);
            let mut block = [q[0], q[1], q[2], q[3], t[0], t[1], t[2], t[3]];
            let pos = |blk: &[u32; 8]| blk.map(|v| self.vertices.get(v as usize).copied().unwrap_or([0.0; 3]));
            if min_corner_det(&pos(&block)) < 0.0 {
                block = [q[0], q[3], q[2], q[1], t[0], t[3], t[2], t[1]];
            }
            made.push(self.blocks.len());
            self.blocks.push(block);
        }
        Ok(made)
    }

    /// Split block `b` along its direction `d` (0 to 2) at the fractions `at` (increasing, between
    /// 0 and 1, measured from its side-0 face), and with it the whole sheet of blocks the cut runs
    /// through: every block holding an edge parallel to the cut edges (chained through shared
    /// edges) is split the same way, so the topology stays conforming. Each new vertex lies on its
    /// edge at the fraction, counted from the end `b`'s side-0 face reaches through the chain. A
    /// split block keeps its index for its first slab (from its own side 0); the other slabs are
    /// appended. Returns the appended blocks.
    pub fn split(&mut self, b: usize, d: usize, at: &[f64]) -> Result<Vec<usize>> {
        if d > 2 {
            return Err(Error::Invalid(format!("a block has directions 0 to 2, not {d}")));
        }
        if at.is_empty() || at.len() > 1000 || !at.iter().all(|t| *t > 0.0 && *t < 1.0) || at.windows(2).any(|w| !(w[1] > w[0])) {
            return Err(Error::Invalid(format!("split fractions must be increasing numbers between 0 and 1, got {at:?}")));
        }
        let seed = self.blocks.get(b).copied().ok_or_else(|| Error::Invalid(format!("there is no block {b}")))?;
        let key = |[p, q]: [u32; 2]| [p.min(q), p.max(q)];
        // The sheet: each cut edge (ascending key) with the vertex its fractions count from, and
        // each block in it with its cut direction, grown through an edge → blocks index.
        let mut by_edge: BTreeMap<[u32; 2], Vec<usize>> = BTreeMap::new();
        for (i, block) in self.blocks.iter().enumerate() {
            for dd in 0..3 {
                for e in edges_along(block, dd) {
                    by_edge.entry(key(e)).or_default().push(i);
                }
            }
        }
        let mut start: BTreeMap<[u32; 2], u32> = BTreeMap::new();
        for e in edges_along(&seed, d) {
            start.insert(key(e), e[0]);
        }
        let mut cut: BTreeMap<usize, usize> = BTreeMap::from([(b, d)]);
        let mut queue = vec![(b, d)];
        while let Some((i, dd)) = queue.pop() {
            let block = self.blocks.get(i).copied().unwrap_or([0; 8]);
            let along = edges_along(&block, dd);
            // The block runs the sheet's way, or against it: read off an edge already in the
            // sheet (the one it was reached through), and orient its other edges to match.
            let against = along.iter().find_map(|e| start.get(&key(*e)).map(|s| *s != e[0])).unwrap_or(false);
            for e in along {
                let k = key(e);
                start.entry(k).or_insert(if against { e[1] } else { e[0] });
                for &n in by_edge.get(&k).map_or(&[][..], Vec::as_slice) {
                    if cut.contains_key(&n) {
                        continue;
                    }
                    let nb = self.blocks.get(n).copied().unwrap_or([0; 8]);
                    if let Some(nd) = (0..3).find(|&x| edges_along(&nb, x).iter().any(|f| key(*f) == k)) {
                        cut.insert(n, nd);
                        queue.push((n, nd));
                    }
                }
            }
        }
        // Each cut block's four edges must all count from its side-0 face or all from its side-1
        // face; a sheet that comes back on itself reversed cannot take uneven fractions.
        let mut flipped: BTreeMap<usize, bool> = BTreeMap::new();
        for (&i, &dd) in &cut {
            let block = self.blocks.get(i).copied().unwrap_or([0; 8]);
            let ways = edges_along(&block, dd).map(|e| start.get(&key(e)) != Some(&e[0]));
            if ways.iter().any(|w| *w != ways[0]) {
                return Err(Error::Invalid(format!(
                    "block {i} would be cut from both ends: the sheet through block {b} comes back on itself reversed; split another block"
                )));
            }
            flipped.insert(i, ways[0]);
        }
        if self.blocks.len() + cut.len() * at.len() > MAX_BLOCKS {
            return Err(Error::Limit(format!("the split makes more than {MAX_BLOCKS} blocks")));
        }
        // The vertices along each cut edge, from its start.
        let mut rows: BTreeMap<[u32; 2], Vec<u32>> = BTreeMap::new();
        for (k, &s) in &start {
            let e = if k[0] == s { k[1] } else { k[0] };
            let (Some(ps), Some(pe)) = (self.vertices.get(s as usize).copied(), self.vertices.get(e as usize).copied()) else {
                return Err(Error::Invalid(format!("vertex {s} or {e} does not exist")));
            };
            let mut row = vec![s];
            for &t in at {
                row.push(self.vertices.len() as u32);
                self.vertices.push(crate::vec3::lerp(ps, pe, t));
            }
            row.push(e);
            rows.insert(*k, row);
        }
        let n = at.len() + 1;
        let mut made = Vec::new();
        for (&i, &dd) in &cut {
            let block = self.blocks.get(i).copied().unwrap_or([0; 8]);
            let rev = flipped.get(&i).copied().unwrap_or(false);
            // Slab `j`, counted from the block's own side 0.
            let slab = |j: usize| -> [u32; 8] {
                let mut out = block;
                for (k, bits) in CORNERS.iter().enumerate() {
                    let (mut low, mut high) = (*bits, *bits);
                    low[dd] = 0;
                    high[dd] = 1;
                    let row = rows.get(&key([block[corner_at(low)], block[corner_at(high)]])).map_or(&[][..], Vec::as_slice);
                    let local = j + bits[dd];
                    out[k] = row.get(if rev { n - local } else { local }).copied().unwrap_or(block[k]);
                }
                out
            };
            if let Some(first) = self.blocks.get_mut(i) {
                *first = slab(0);
            }
            for j in 1..n {
                made.push(self.blocks.len());
                self.blocks.push(slab(j));
            }
        }
        Ok(made)
    }

    /// Remove blocks (the others keep their order and are renumbered).
    pub fn remove(&mut self, blocks: &[usize]) -> Result<()> {
        if let Some(b) = blocks.iter().find(|b| **b >= self.blocks.len()) {
            return Err(Error::Invalid(format!("there is no block {b}")));
        }
        let mut k = 0;
        self.blocks.retain(|_| {
            k += 1;
            !blocks.contains(&(k - 1))
        });
        self.compact();
        if self.blocks.is_empty() { Err(Error::Invalid("every block was removed".into())) } else { Ok(()) }
    }

    /// Extrude block faces `(block, hex8 face)` by `offset` in `layers` blocks each; the faces'
    /// shared vertices stay shared. Returns the new blocks' indices, the last layer's outer face
    /// is face 2 of the blocks made from it.
    pub fn extrude(&mut self, faces: &[(usize, u8)], offset: P, layers: usize) -> Result<Vec<usize>> {
        if !(1..=1000).contains(&layers) {
            return Err(Error::Invalid(format!("an extrusion has 1 to 1000 layers, got {layers}")));
        }
        if !offset.iter().all(|x| x.is_finite()) {
            return Err(Error::Invalid("the extrusion offset must be finite".into()));
        }
        let mut quads = Vec::new();
        for &(b, f) in faces {
            let (d, s) =
                FACE_SIDES.get(usize::from(f).wrapping_sub(1)).copied().ok_or_else(|| Error::Invalid(format!("face {f} is not a hex8 face")))?;
            let block = self.blocks.get(b).copied().ok_or_else(|| Error::Invalid(format!("there is no block {b}")))?;
            let (e1, e2) = match d {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let q = [[0, 0], [1, 0], [1, 1], [0, 1]].map(|[a, bb]| {
                let mut bits = [0; 3];
                bits[d] = s;
                bits[e1] = a;
                bits[e2] = bb;
                block[corner_at(bits)]
            });
            quads.push(q);
        }
        let mut made = Vec::new();
        let mut layer_of: BTreeMap<u32, u32> = BTreeMap::new();
        let mut base: Vec<[u32; 4]> = quads;
        for _ in 0..layers {
            layer_of.clear();
            let mut next = Vec::new();
            for q in &base {
                let mut top = [0u32; 4];
                for (t, &v) in top.iter_mut().zip(q) {
                    *t = match layer_of.get(&v) {
                        Some(&w) => w,
                        None => {
                            let p = self.vertices.get(v as usize).copied().ok_or_else(|| Error::Invalid(format!("vertex {v} does not exist")))?;
                            let w = self.vertices.len() as u32;
                            self.vertices.push(add(p, scale(offset, 1.0 / layers as f64)));
                            layer_of.insert(v, w);
                            w
                        }
                    };
                }
                let mut block = [q[0], q[1], q[2], q[3], top[0], top[1], top[2], top[3]];
                let pos = |blk: &[u32; 8]| blk.map(|v| self.vertices.get(v as usize).copied().unwrap_or([0.0; 3]));
                if min_corner_det(&pos(&block)) < 0.0 {
                    block = [q[0], q[3], q[2], q[1], top[0], top[3], top[2], top[1]];
                }
                made.push(self.blocks.len());
                self.blocks.push(block);
                next.push(top);
            }
            base = next;
        }
        Ok(made)
    }

    /// The seven-block sphere: an inner cube with its corners at `inset · radius` from the
    /// centre and six blocks between the cube's faces and the sphere (the cube-to-sphere O-grid).
    pub fn sphere(center: P, radius: f64, inset: f64) -> Result<Topology> {
        if !(radius > 0.0 && radius.is_finite()) || !center.iter().all(|x| x.is_finite()) {
            return Err(Error::invalid(crate::message![
                "a sphere needs a finite centre and a positive radius, got ",
                crate::message::length(radius)
            ]));
        }
        let h = radius / 3f64.sqrt();
        let mut t = Topology::grid([&[center[0] - h, center[0] + h], &[center[1] - h, center[1] + h], &[center[2] - h, center[2] + h]])?;
        t.o_grid(0, &[], [inset; 3])?;
        Ok(t)
    }

    /// The butterfly (five-block O-grid of a disc) of radius `radius` about `center` in the plane
    /// normal to `axis`, stacked through `levels` (distances along `axis` from `center`,
    /// increasing): a centre square with its corners at `inset · radius` and four blocks around
    /// it out to the circle, per level.
    pub fn butterfly(center: P, axis: P, radius: f64, levels: &[f64], inset: f64) -> Result<Topology> {
        if !(radius > 0.0 && radius.is_finite()) || !center.iter().all(|x| x.is_finite()) {
            return Err(Error::invalid(crate::message![
                "a butterfly needs a finite centre and a positive radius, got ",
                crate::message::length(radius)
            ]));
        }
        let z = crate::vec3::norm(axis);
        if !(z > 0.0 && z.is_finite()) {
            return Err(Error::Invalid("the butterfly's axis must be a nonzero direction".into()));
        }
        let w = scale(axis, 1.0 / z);
        let t0 = if w[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u = {
            let c = crate::vec3::cross(w, t0);
            scale(c, 1.0 / crate::vec3::norm(c))
        };
        let v = crate::vec3::cross(w, u);
        let h = radius / 2f64.sqrt();
        let grid = Topology::grid([&[-h, h], &[-h, h], levels])?;
        let mut t = Topology {
            vertices: grid.vertices.iter().map(|p| add(center, add(add(scale(u, p[0]), scale(v, p[1])), scale(w, p[2])))).collect(),
            blocks: grid.blocks,
        };
        for b in 0..levels.len() - 1 {
            t.o_grid(b, &[1, 2], [inset; 3])?;
        }
        t.merge(MERGE * (radius + levels.iter().fold(0.0f64, |m, x| m.max(x.abs()))));
        Ok(t)
    }
}
