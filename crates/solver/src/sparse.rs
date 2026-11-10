//! Sparse symmetric matrices in compressed rows, both triangles stored.
//!
//! Assembly works on the node graph: a node's row holds every node it shares an element with,
//! and each node pair is a dense block of their degrees of freedom ([`DofLayout`]: a node's are
//! contiguous), so the pattern is built once from the element connectivity and the stiffness and
//! mass share it. Sums run in a fixed order (elements by
//! label, then the element's own node order), so assembled values are the same bit for bit on
//! every run and platform.

use fieldcraft_fem::{DofLayout, DofSet};

use crate::{Error, Result};

/// Rows per parallel block of a product.
const ROWS: usize = 4096;

/// A square sparse matrix in compressed rows with sorted column indices.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    pub n: usize,
    pub ptr: Vec<usize>,
    pub idx: Vec<usize>,
    pub val: Vec<f64>,
}

impl Csr {
    /// The row's column indices and values.
    pub fn row(&self, r: usize) -> (&[usize], &[f64]) {
        let (a, b) = (self.ptr.get(r).copied().unwrap_or(0), self.ptr.get(r + 1).copied().unwrap_or(0));
        (self.idx.get(a..b).unwrap_or_default(), self.val.get(a..b).unwrap_or_default())
    }

    /// y = A x.
    /// Rows are computed in parallel blocks; each row's sum runs in column order.
    pub fn mul(&self, x: &[f64], y: &mut [f64]) {
        self.mul_with_rows(x, y, ROWS);
    }

    /// Independent row tasks; changing their size never changes a row's column-order sum.
    pub(crate) fn mul_with_rows(&self, x: &[f64], y: &mut [f64], rows: usize) {
        let rows = rows.max(1);
        let n = self.n.min(y.len());
        fieldcraft_parallel::chunks_mut(y.get_mut(..n).unwrap_or_default(), rows, |b, block| {
            for (i, yr) in block.iter_mut().enumerate() {
                let (cols, vals) = self.row(b * rows + i);
                *yr = cols.iter().zip(vals).map(|(c, v)| v * x.get(*c).copied().unwrap_or(0.0)).sum();
            }
        });
    }

    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n)
            .map(|r| {
                let (cols, vals) = self.row(r);
                cols.binary_search(&r).ok().and_then(|i| vals.get(i)).copied().unwrap_or(0.0)
            })
            .collect()
    }

    /// Make the matrix exactly symmetric by copying the lower triangle (what the factorization
    /// reads) onto the upper one; the pattern must be structurally symmetric.
    pub fn mirror_lower(&mut self) {
        for r in 0..self.n {
            let (start, end) = (self.ptr.get(r).copied().unwrap_or(0), self.ptr.get(r + 1).copied().unwrap_or(0));
            for at in start..end {
                let Some(&c) = self.idx.get(at) else { continue };
                if c <= r {
                    continue;
                }
                let (cs, ce) = (self.ptr.get(c).copied().unwrap_or(0), self.ptr.get(c + 1).copied().unwrap_or(0));
                let mirror = self.idx.get(cs..ce).and_then(|cols| cols.binary_search(&r).ok()).and_then(|i| self.val.get(cs + i)).copied();
                if let (Some(v), Some(x)) = (mirror, self.val.get_mut(at)) {
                    *x = v;
                }
            }
        }
    }

    /// `self + a·other` on the same pattern.
    pub fn add_scaled(&self, a: f64, other: &Csr) -> Result<Csr> {
        if self.ptr != other.ptr || self.idx != other.idx {
            return Err(Error::Internal("matrices with different patterns".into()));
        }
        let val = self.val.iter().zip(&other.val).map(|(x, y)| x + a * y).collect();
        Ok(Csr { n: self.n, ptr: self.ptr.clone(), idx: self.idx.clone(), val })
    }
}

/// The node-block pattern of a mesh: for each node, the sorted nodes it shares an element with
/// (itself included).
#[derive(Debug, Clone)]
pub struct NodeGraph {
    pub adjacency: Vec<Vec<usize>>,
}

impl NodeGraph {
    pub fn new(nodes: usize, elements: impl Iterator<Item = impl AsRef<[usize]>>) -> NodeGraph {
        let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); nodes];
        for e in elements {
            let e = e.as_ref();
            for &a in e {
                if let Some(row) = adjacency.get_mut(a) {
                    row.extend_from_slice(e);
                }
            }
        }
        for row in &mut adjacency {
            row.sort_unstable();
            row.dedup();
        }
        NodeGraph { adjacency }
    }

    /// The graph with more blocks, each coupling all its slots (contact: the slots a contact
    /// element couples change with the deformation).
    pub fn with_blocks<'a>(&self, blocks: impl Iterator<Item = &'a [usize]>) -> NodeGraph {
        let mut adjacency = self.adjacency.clone();
        let mut touched = Vec::new();
        for b in blocks {
            for &a in b {
                if let Some(row) = adjacency.get_mut(a) {
                    row.extend_from_slice(b);
                    touched.push(a);
                }
            }
        }
        touched.sort_unstable();
        touched.dedup();
        for a in touched {
            if let Some(row) = adjacency.get_mut(a) {
                row.sort_unstable();
                row.dedup();
            }
        }
        NodeGraph { adjacency }
    }

    /// Nonzeros of the dof-level matrix: a dense block per node pair.
    pub fn entries(&self, layout: &DofLayout) -> usize {
        self.adjacency.iter().enumerate().map(|(a, row)| layout.range(a).len() * row.iter().map(|b| layout.range(*b).len()).sum::<usize>()).sum()
    }

    /// A zero matrix over the degrees of freedom of the graph's layout slots (the mesh nodes
    /// `0..nodes`, and reference points when the graph has rows for them) on this pattern.
    pub fn zero_matrix(&self, layout: &DofLayout) -> Csr {
        let n = self.adjacency.len().checked_sub(1).map_or(0, |last| layout.range(last).end);
        let mut ptr = Vec::with_capacity(n + 1);
        ptr.push(0);
        let mut idx = Vec::with_capacity(self.entries(layout));
        for (a, row) in self.adjacency.iter().enumerate() {
            for _ in layout.range(a) {
                for &b in row {
                    idx.extend(layout.range(b));
                }
                ptr.push(idx.len());
            }
        }
        let val = vec![0.0; idx.len()];
        Csr { n, ptr, idx, val }
    }
}

/// Most kinds a node carries ([`fieldcraft_fem::DofKind::ALL`]).
const MAX_KINDS: usize = 8;

/// Where node b's block lies in node a's rows of a matrix built by [`NodeGraph::zero_matrix`],
/// for the degrees of freedom `kinds` of both: the value index of (row r, column c) is
/// `rows[r] + cols[c]`. The rows of a node share one pattern, so this is one search per node
/// pair.
struct Block {
    rows: [usize; MAX_KINDS],
    cols: [usize; MAX_KINDS],
}

impl Block {
    fn new(k: &Csr, layout: &DofLayout, a: usize, b: usize, kinds: DofSet) -> Result<Block> {
        Self::from_pattern(&k.ptr, &k.idx, layout, a, b, kinds)
    }

    fn from_pattern(ptr: &[usize], idx: &[usize], layout: &DofLayout, a: usize, b: usize, kinds: DofSet) -> Result<Block> {
        let first = layout.range(a).start;
        let (start, end) = (ptr.get(first).copied().unwrap_or(0), ptr.get(first + 1).copied().unwrap_or(0));
        let pattern = idx.get(start..end).ok_or_else(|| Error::Internal("pattern".into()))?;
        let at = pattern.binary_search(&layout.range(b).start).map_err(|_| Error::Internal("entry outside the pattern".into()))?;
        let missing = |kind| Error::Internal(format!("an element degree of freedom ({kind:?}) its node does not have"));
        let mut block = Block { rows: [0; MAX_KINDS], cols: [0; MAX_KINDS] };
        for ((row, col), kind) in block.rows.iter_mut().zip(&mut block.cols).zip(kinds.iter()) {
            let r = layout.dof(a, kind).ok_or_else(|| missing(kind))?;
            *row = ptr.get(r).copied().unwrap_or(0) + at;
            *col = layout.kinds(b).position(kind).ok_or_else(|| missing(kind))?;
        }
        Ok(block)
    }
}

fn add_at(k: &mut Csr, at: usize, v: f64) {
    if let Some(x) = k.val.get_mut(at) {
        *x += v;
    }
}

/// Add an element matrix into a matrix built by [`NodeGraph::zero_matrix`]: row-major over the
/// element's degrees of freedom, node by node and `kinds` in order within a node; `nodes` are
/// the element's layout slots.
pub fn add_block(k: &mut Csr, layout: &DofLayout, nodes: &[usize], kinds: DofSet, data: &[f64]) -> Result<()> {
    let per = kinds.len();
    let width = per * nodes.len();
    if data.len() != width * width {
        return Err(Error::Internal("element matrix size".into()));
    }
    for (i, &a) in nodes.iter().enumerate() {
        for (j, &b) in nodes.iter().enumerate() {
            let block = Block::new(k, layout, a, b, kinds)?;
            for r in 0..per {
                for c in 0..per {
                    add_at(k, block.rows[r] + block.cols[c], data.get((per * i + r) * width + per * j + c).copied().unwrap_or(0.0));
                }
            }
        }
    }
    Ok(())
}

/// Add an n × n element matrix to the block of each of `kinds` (a solid's consistent mass on its
/// translations).
pub fn add_scalar(k: &mut Csr, layout: &DofLayout, nodes: &[usize], kinds: DofSet, data: &[f64]) -> Result<()> {
    let n = nodes.len();
    if data.len() != n * n {
        return Err(Error::Internal("element matrix size".into()));
    }
    for (i, &a) in nodes.iter().enumerate() {
        for (j, &b) in nodes.iter().enumerate() {
            let block = Block::new(k, layout, a, b, kinds)?;
            let v = data.get(i * n + j).copied().unwrap_or(0.0);
            for c in 0..kinds.len() {
                add_at(k, block.rows[c] + block.cols[c], v);
            }
        }
    }
    Ok(())
}

/// Add a per-node diagonal to each of `kinds` (a solid's lumped mass on its translations).
pub fn add_diagonal(k: &mut Csr, layout: &DofLayout, nodes: &[usize], kinds: DofSet, diagonal: &[f64]) -> Result<()> {
    for (&a, v) in nodes.iter().zip(diagonal) {
        let block = Block::new(k, layout, a, a, kinds)?;
        for c in 0..kinds.len() {
            add_at(k, block.rows[c] + block.cols[c], *v);
        }
    }
    Ok(())
}

/// Matrix shape in a bounded assembly round.
#[derive(Clone, Copy)]
pub(crate) enum AssemblyShape {
    Block,
    Scalar,
    Diagonal,
}

/// Borrowed element data, in the original element-label order.
pub(crate) struct AssemblyMatrix<'a> {
    pub nodes: &'a [usize],
    pub kinds: DofSet,
    pub data: &'a [f64],
    pub shape: AssemblyShape,
}

/// Disjoint node-row owners scatter a round without atomics or partial floating sums.
/// Each owner's incidences retain element order then local node order. Thus every
/// entry sees exactly the scalar assembly sequence, including repeated nodes.
/// Fixed 256-node ranges are independent of the worker count (Hughes 1987 assembly).
pub(crate) fn add_round(k: &mut Csr, layout: &DofLayout, matrices: &[AssemblyMatrix<'_>]) -> Result<()> {
    const NODES: usize = 256;
    let slots = matrices.iter().flat_map(|m| m.nodes.iter()).copied().max().map_or(0, |n| n + 1);
    let mut incidence = vec![Vec::new(); slots.div_ceil(NODES)];
    for (e, m) in matrices.iter().enumerate() {
        let n = m.nodes.len();
        let expected = match m.shape {
            AssemblyShape::Block => (n * m.kinds.len()).pow(2),
            AssemblyShape::Scalar => n * n,
            AssemblyShape::Diagonal => n,
        };
        if m.data.len() != expected {
            return Err(Error::Internal("element matrix size".into()));
        }
        for (i, &node) in m.nodes.iter().enumerate() {
            if let Some(list) = incidence.get_mut(node / NODES) {
                list.push((e, i));
            }
        }
    }
    let mut remaining = k.val.as_mut_slice();
    let mut offset = 0;
    let mut tasks = Vec::with_capacity(incidence.len());
    for (b, list) in incidence.into_iter().enumerate() {
        let end_row = layout.range(((b + 1) * NODES).min(slots) - 1).end;
        let end = k.ptr.get(end_row).copied().ok_or_else(|| Error::Internal("assembly row range".into()))?;
        let len = end.checked_sub(offset).filter(|len| *len <= remaining.len()).ok_or_else(|| Error::Internal("assembly row range".into()))?;
        let (values, rest) = remaining.split_at_mut(len);
        tasks.push((offset, values, list, None));
        remaining = rest;
        offset = end;
    }
    let ptr = &k.ptr;
    let idx = &k.idx;
    fieldcraft_parallel::for_each_mut(&mut tasks, |_, (offset, values, list, failed)| {
        for &(e, i) in list.iter() {
            let Some(m) = matrices.get(e) else { continue };
            let Some(&a) = m.nodes.get(i) else { continue };
            let per = m.kinds.len();
            for (j, &b) in m.nodes.iter().enumerate() {
                if matches!(m.shape, AssemblyShape::Diagonal) && j != i {
                    continue;
                }
                let block = match Block::from_pattern(ptr, idx, layout, a, b, m.kinds) {
                    Ok(block) => block,
                    Err(error) => {
                        *failed = Some((e, i, j, error));
                        return;
                    }
                };
                for r in 0..per {
                    for c in 0..per {
                        if !matches!(m.shape, AssemblyShape::Block) && r != c {
                            continue;
                        }
                        let at = block.rows[r] + block.cols[c];
                        let source = match m.shape {
                            AssemblyShape::Block => (per * i + r) * (per * m.nodes.len()) + per * j + c,
                            AssemblyShape::Scalar => i * m.nodes.len() + j,
                            AssemblyShape::Diagonal => i,
                        };
                        if let Some(value) = at.checked_sub(*offset).and_then(|p| values.get_mut(p)) {
                            *value += m.data.get(source).copied().unwrap_or(0.0);
                        }
                    }
                }
            }
        }
    });
    // Worker completion order cannot select a different first structural error.
    let failure = tasks.into_iter().filter_map(|(_, _, _, e)| e).min_by_key(|(e, i, j, _)| (*e, *i, *j));
    if let Some((_, _, _, error)) = failure {
        return Err(error);
    }
    Ok(())
}

/// Add `v` to entry (`i`, `j`), which must be in the pattern.
pub fn add_entry(k: &mut Csr, i: usize, j: usize, v: f64) -> Result<()> {
    let (start, end) = (k.ptr.get(i).copied().unwrap_or(0), k.ptr.get(i + 1).copied().unwrap_or(0));
    let at = k.idx.get(start..end).and_then(|cols| cols.binary_search(&j).ok()).ok_or_else(|| Error::Internal("entry outside the pattern".into()))?;
    add_at(k, start + at, v);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembly_row_owners_preserve_scalar_bits() {
        let mut builder = DofLayout::builder();
        for i in 0..777 {
            builder.node(i, DofSet::TRANSLATIONS);
        }
        let layout = builder.build();
        let nodes: Vec<Vec<usize>> = (0..700).map(|i| vec![i, i + 1, i + 65, i]).collect();
        let graph = NodeGraph::new(777, nodes.iter());
        for shape in [AssemblyShape::Block, AssemblyShape::Scalar, AssemblyShape::Diagonal] {
            let size = match shape {
                AssemblyShape::Block => 144,
                AssemblyShape::Scalar => 16,
                AssemblyShape::Diagonal => 4,
            };
            let data: Vec<Vec<f64>> = nodes
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    (0..size)
                        .map(|j| match (i + j) % 3 {
                            0 => 1e16,
                            1 => -1e16,
                            _ => 0.25,
                        })
                        .collect()
                })
                .collect();
            let matrices: Vec<_> =
                nodes.iter().zip(&data).map(|(nodes, data)| AssemblyMatrix { nodes, kinds: DofSet::TRANSLATIONS, data, shape }).collect();
            let mut reference = graph.zero_matrix(&layout);
            for m in &matrices {
                match shape {
                    AssemblyShape::Block => add_block(&mut reference, &layout, m.nodes, m.kinds, m.data),
                    AssemblyShape::Scalar => add_scalar(&mut reference, &layout, m.nodes, m.kinds, m.data),
                    AssemblyShape::Diagonal => add_diagonal(&mut reference, &layout, m.nodes, m.kinds, m.data),
                }
                .unwrap();
            }
            for threads in [1, 2, 7] {
                let mut actual = graph.zero_matrix(&layout);
                fieldcraft_parallel::with_threads(threads, || add_round(&mut actual, &layout, &matrices)).unwrap();
                assert_eq!(actual.ptr, reference.ptr);
                assert_eq!(actual.idx, reference.idx);
                assert!(actual.val.iter().zip(&reference.val).all(|(a, b)| a.to_bits() == b.to_bits()));
            }
        }
    }

    #[test]
    fn row_task_size_preserves_product_bits_and_partial_output() {
        let n = 8193;
        let mut a = Csr { n, ptr: vec![0], idx: Vec::new(), val: Vec::new() };
        for r in 0..n {
            for c in r.saturating_sub(2)..(r + 3).min(n) {
                a.idx.push(c);
                a.val.push(if c % 2 == 0 { 1e12 } else { -1e-12 });
            }
            a.ptr.push(a.idx.len());
        }
        let x: Vec<f64> = (0..n).map(|i| (i % 23) as f64 - 11.0).collect();
        for len in [n - 1, n, n + 1] {
            let mut reference = vec![7.0; len];
            a.mul(&x, &mut reference);
            for threads in [1, 2, 7] {
                let mut result = vec![7.0; len];
                fieldcraft_parallel::with_threads(threads, || a.mul_with_rows(&x, &mut result, 512));
                assert!(reference.iter().zip(result).all(|(a, b)| a.to_bits() == b.to_bits()));
            }
        }
    }

    #[test]
    fn assemble_mirror_multiply() {
        // Two "elements" sharing node 1.
        let g = NodeGraph::new(3, [vec![0usize, 1], vec![1, 2]].into_iter());
        assert_eq!(g.adjacency, vec![vec![0, 1], vec![0, 1, 2], vec![1, 2]]);
        let layout = DofLayout::uniform(0..3, DofSet::TRANSLATIONS);
        let mut k = g.zero_matrix(&layout);
        let mut e = vec![0.0; 36];
        for i in 0..6 {
            e[i * 6 + i] = 2.0;
            e[i * 6 + (i + 3) % 6] = -1.0;
        }
        add_block(&mut k, &layout, &[0, 1], DofSet::TRANSLATIONS, &e).unwrap();
        add_block(&mut k, &layout, &[1, 2], DofSet::TRANSLATIONS, &e).unwrap();
        assert_eq!(k.diagonal(), vec![2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0]);
        let x = vec![1.0; 9];
        let mut y = vec![0.0; 9];
        k.mul(&x, &mut y);
        assert_eq!(y, vec![1.0; 3].into_iter().chain([2.0; 3]).chain([1.0; 3]).collect::<Vec<_>>());
        // Mirroring keeps a symmetric matrix as it is.
        let mut mirrored = k.clone();
        mirrored.mirror_lower();
        assert_eq!(mirrored, k);
        assert!(add_block(&mut k, &layout, &[0, 2], DofSet::TRANSLATIONS, &e).is_err(), "nodes 0 and 2 share no element");
        // Scalar (per component) and diagonal additions land on each component.
        let mut m = g.zero_matrix(&layout);
        add_scalar(&mut m, &layout, &[0, 1], DofSet::TRANSLATIONS, &[2.0, 1.0, 1.0, 2.0]).unwrap();
        add_diagonal(&mut m, &layout, &[2], DofSet::TRANSLATIONS, &[5.0]).unwrap();
        assert_eq!(m.diagonal(), vec![2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0]);
        assert_eq!(m.row(0).1, &[2.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        let s = k.add_scaled(-1.0, &k).unwrap();
        assert!(s.val.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn mixed_kinds_per_node_assemble_into_their_blocks() {
        // A solid on nodes 0 and 1 (translations) and a beam-like element on 1 and 2 (six each):
        // node 1 carries six, and the solid's block lands on its translations.
        let mut b = DofLayout::builder();
        b.node(10, DofSet::TRANSLATIONS).node(11, DofSet::STRUCTURAL).node(12, DofSet::STRUCTURAL);
        let layout = b.build();
        let g = NodeGraph::new(3, [vec![0usize, 1], vec![1, 2]].into_iter());
        let mut k = g.zero_matrix(&layout);
        assert_eq!((k.n, k.idx.len(), g.entries(&layout)), (15, 3 * 9 + 6 * 15 + 6 * 12, 189));
        let solid: Vec<f64> = (0..36).map(|i| i as f64).collect();
        add_block(&mut k, &layout, &[0, 1], DofSet::TRANSLATIONS, &solid).unwrap();
        let beam = vec![1.0; 144];
        add_block(&mut k, &layout, &[1, 2], DofSet::STRUCTURAL, &beam).unwrap();
        // Row Uy of node 1 (global 4): the solid's row 4 on node 0's translations and node 1's.
        let (cols, vals) = k.row(4);
        let at = |c: usize| vals[cols.binary_search(&c).unwrap()];
        assert_eq!([at(0), at(1), at(2)], [24.0, 25.0, 26.0]);
        assert_eq!([at(3), at(4), at(5)], [27.0 + 1.0, 28.0 + 1.0, 29.0 + 1.0]);
        assert_eq!(at(6), 1.0, "rotation Rx of node 1 from the beam only");
        // A translations element cannot claim rotations its node lacks.
        assert!(add_block(&mut k, &layout, &[0, 1], DofSet::STRUCTURAL, &beam).is_err());
    }
}
