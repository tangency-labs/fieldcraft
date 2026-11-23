//! Resonant modes of an empty rectangular cavity with perfectly conducting walls, from the
//! Nédélec elements of `fieldcraft_fem::nedelec`, against the closed form (D. M. Pozar,
//! *Microwave Engineering*, 4th ed., Wiley 2012, §6.4; docs/rfc/0001-electromagnetics.md §4):
//!
//! k²_mnl = (mπ/a)² + (nπ/b)² + (lπ/d)²,
//!
//! one mode for each (m, n, l) with at most one index zero (TE or TM to the axis of that index),
//! two when none is (a TE and a TM mode). The weak form curl curl E = k² E with E × n = 0 on the
//! walls gives K x = k² M x, K the curl–curl and M the mass matrix of the free edge unknowns.
//!
//! What is checked:
//! - **the null space** (RFC §8.4): every gradient of a function that vanishes on the walls is a
//!   solution with k = 0. Its dimension is the number of interior nodes at order 1, plus the
//!   number of free edges at order 2 (the gradients of the edges' quadratics). The inertia of
//!   K − σM for 0 < σ < k²₁ counts exactly those as its negative eigenvalues (Sylvester's law):
//!   no fewer (a mode lost) and no more (a spurious mode below the first resonance);
//! - **the resonances** above the shift by shift-invert Lanczos on the LDLᵀ factor, their
//!   multiplicities included, against the closed form;
//! - **convergence**: the eigenvalue error falls as h² at order 1 and h⁴ at order 2, the
//!   rate 2k of Nédélec elements of degree k (D. Boffi, *Finite element approximation of
//!   eigenvalue problems*, Acta Numer. 19:1–120, 2010, §19–20).
//!
//! The system is assembled as a solver assembles one: the edges and faces numbered by
//! `DofLayout`, the pattern of their blocks (`NodeGraph`), the element matrices added by
//! `add_block`, the unknowns ordered by nested dissection over the edges' midpoints and the faces'
//! centroids, and the walls' unknowns removed. The discrete gradient (`nedelec::gradient`) checks
//! the numbering: K G = 0.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use fieldcraft_fem::dof::{DofLayout, Entity};
use fieldcraft_fem::nedelec::{EdgeTet, Order, Potential, declare, gradient, slots};
use fieldcraft_solver::cholesky::UpperColumns;
use fieldcraft_solver::eigen::{Operator, Wanted, lanczos};
use fieldcraft_solver::ldlt::LdltFactor;
use fieldcraft_solver::ordering::fill_reducing_order;
use fieldcraft_solver::sparse::{Csr, NodeGraph, add_block};

/// The cavity's sides along x, y and z (metres).
const SIDES: [f64; 3] = [1.0, 0.5, 0.75];

/// The box split into `cells` cubes per axis, each into six tetrahedra around its main diagonal
/// (the same diagonal everywhere, so the faces match): the nodes, and the tetrahedra by node.
fn mesh(cells: [usize; 3]) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let [nx, ny, nz] = cells;
    let node = |i: usize, j: usize, k: usize| (i + (nx + 1) * (j + (ny + 1) * k)) as u32;
    let mut xyz = Vec::new();
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                xyz.push([SIDES[0] * i as f64 / nx as f64, SIDES[1] * j as f64 / ny as f64, SIDES[2] * k as f64 / nz as f64]);
            }
        }
    }
    let mut tets = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                // From corner (0, 0, 0) to (1, 1, 1) of the cube, one axis at a time.
                for axes in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
                    let mut at = [i, j, k];
                    let mut tet = [node(i, j, k); 4];
                    for (step, axis) in axes.into_iter().enumerate() {
                        at[axis] += 1;
                        tet[step + 1] = node(at[0], at[1], at[2]);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (xyz, tets)
}

/// The walls a point lies on, as (axis, side) bits.
fn walls(p: [f64; 3]) -> u8 {
    (0..3).fold(0, |w, a| w | if p[a] == 0.0 { 1 << (2 * a) } else { 0 } | if p[a] == SIDES[a] { 2 << (2 * a) } else { 0 })
}

/// A cavity's matrices on its free unknowns and the null space's dimension.
struct Cavity {
    k: Csr,
    m: Csr,
    null: usize,
}

/// The rows and columns `keep` of `a`, in that order.
fn restrict(a: &Csr, keep: &[usize]) -> Csr {
    let mut new = vec![usize::MAX; a.n];
    for (k, d) in keep.iter().enumerate() {
        new[*d] = k;
    }
    let (mut ptr, mut idx, mut val) = (vec![0], Vec::new(), Vec::new());
    for d in keep {
        let (cols, vals) = a.row(*d);
        let mut row: Vec<(usize, f64)> = cols.iter().zip(vals).filter(|(c, _)| new[**c] != usize::MAX).map(|(c, v)| (new[*c], *v)).collect();
        row.sort_by_key(|e| e.0);
        idx.extend(row.iter().map(|e| e.0));
        val.extend(row.iter().map(|e| e.1));
        ptr.push(idx.len());
    }
    Csr { n: keep.len(), ptr, idx, val }
}

fn cavity(order: Order, cells: [usize; 3]) -> Cavity {
    let (xyz, tets) = mesh(cells);
    let at = |l: u32| xyz[l as usize];
    let mut b = DofLayout::builder();
    for t in &tets {
        declare(&mut b, *t, order).unwrap();
    }
    let layout = b.build();
    let element_slots: Vec<Vec<usize>> = tets.iter().map(|t| slots(&layout, *t, order).unwrap()).collect();
    let graph = NodeGraph::new(layout.slots(), element_slots.iter());
    let (mut k, mut m) = (graph.zero_matrix(&layout), graph.zero_matrix(&layout));
    for (t, s) in tets.iter().zip(&element_slots) {
        let (e, _) = EdgeTet::of_labels(order, t.map(at), *t).unwrap();
        add_block(&mut k, &layout, s, order.kinds(), &e.curl_curl([1.0; 3]).unwrap().data).unwrap();
        add_block(&mut m, &layout, s, order.kinds(), &e.mass([1.0; 3]).unwrap().data).unwrap();
    }
    // Every column of the discrete gradient is in the kernel of K, before the walls are removed:
    // one per node, and at order 2 one per edge. A column has a few entries and K is symmetric,
    // so K g is the sum of K's rows at them.
    let mut columns: BTreeMap<Potential, Vec<(usize, f64)>> = BTreeMap::new();
    for (r, c, v) in gradient(&layout) {
        columns.entry(c).or_default().push((r, v));
    }
    let edges = columns.keys().filter(|c| matches!(c, Potential::Edge(_))).count();
    assert_eq!((columns.len() - edges, edges), (xyz.len(), if order == Order::Two { layout.edges() } else { 0 }));
    let scale = k.val.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let mut kg = vec![0.0; layout.dofs()];
    for (column, entries) in &columns {
        kg.fill(0.0);
        for (r, v) in entries {
            let (cols, vals) = k.row(*r);
            for (c, kv) in cols.iter().zip(vals) {
                kg[*c] += v * kv;
            }
        }
        assert!(kg.iter().all(|v| v.abs() < 1e-12 * scale), "{order:?}: K G ≠ 0 for {column:?}");
    }
    // The order of the slots by their points: nodes, edges' midpoints, faces' centroids.
    let centre = |ls: &[u32]| {
        let n = ls.len() as f64;
        ls.iter().fold([0.0; 3], |p, l| [0, 1, 2].map(|a| p[a] + at(*l)[a] / n))
    };
    let points: Vec<[f64; 3]> = (0..layout.slots())
        .map(|s| match layout.entity(s).unwrap() {
            Entity::Node(l) | Entity::Point(l) => at(l),
            Entity::Edge(e) => centre(&e),
            Entity::Face(f) => centre(&f),
        })
        .collect();
    let ordered = fill_reducing_order(&points, &graph.adjacency);
    // An edge or face is held when it lies in a wall: its nodes share one.
    let on_wall = |ls: &[u32]| ls.iter().fold(0x3f, |w, l| w & walls(at(*l))) != 0;
    let held = |s: usize| match layout.entity(s).unwrap() {
        Entity::Edge(e) => on_wall(&e),
        Entity::Face(f) => on_wall(&f),
        _ => true,
    };
    let free: Vec<usize> = ordered.iter().filter(|s| !held(**s)).flat_map(|s| layout.range(*s)).collect();
    // The null space: the interior nodes' potentials, and at order 2 the free edges' functions.
    let interior = xyz.iter().filter(|p| walls(**p) == 0).count();
    let free_edges = (0..layout.slots()).filter(|s| matches!(layout.entity(*s), Some(Entity::Edge(_))) && !held(*s)).count();
    Cavity { k: restrict(&k, &free), m: restrict(&m, &free), null: interior + if order == Order::Two { free_edges } else { 0 } }
}

/// The closed form's k², ascending, each as often as its multiplicity.
fn resonances(count: usize) -> Vec<f64> {
    let mut out = Vec::new();
    for m in 0..6usize {
        for n in 0..6usize {
            for l in 0..6usize {
                let zeros = [m, n, l].iter().filter(|i| **i == 0).count();
                if zeros > 1 {
                    continue;
                }
                let k2 = (m as f64 * PI / SIDES[0]).powi(2) + (n as f64 * PI / SIDES[1]).powi(2) + (l as f64 * PI / SIDES[2]).powi(2);
                out.extend(std::iter::repeat_n(k2, if zeros == 0 { 2 } else { 1 }));
            }
        }
    }
    out.sort_by(f64::total_cmp);
    out.truncate(count);
    out
}

/// Shift-invert: M products and solves with K − σM.
struct ShiftInvert<'a> {
    m: &'a Csr,
    factor: LdltFactor,
}

impl Operator for ShiftInvert<'_> {
    fn n(&self) -> usize {
        self.m.n
    }

    fn mass(&self, x: &[f64], y: &mut [f64]) {
        self.m.mul(x, y);
    }

    fn solve(&self, x: &mut [f64]) -> fieldcraft_solver::Result<()> {
        self.factor.solve(x)
    }
}

/// The negative eigenvalues of K − σM and the first `count` k² above σ.
fn modes(c: &Cavity, sigma: f64, count: usize) -> (usize, Vec<f64>) {
    let shifted = c.k.add_scaled(-sigma, &c.m).unwrap();
    let factor = LdltFactor::new(&UpperColumns::from_symmetric(&shifted), usize::MAX).unwrap();
    let negative = factor.inertia().negative;
    let op = ShiftInvert { m: &c.m, factor };
    let pairs = lanczos(&op, sigma, Wanted { count, min: sigma, max: f64::INFINITY }, 6, c.m.n, 1e-10).unwrap();
    (negative, pairs.values)
}

/// The largest relative error of the computed k² against the closed form.
fn error(computed: &[f64], exact: &[f64]) -> f64 {
    assert_eq!(computed.len(), exact.len(), "{computed:?}");
    computed.iter().zip(exact).map(|(c, e)| ((c - e) / e).abs()).fold(0.0, f64::max)
}

/// The six lowest k² of [`SIDES`] on its mesh of `[4, 2, 3]` cells, at orders 1 and 2: our
/// result, frozen once it matched the closed form (largest errors 5.5 % and 0.81 %, falling at
/// the rates of `the_eigenvalue_error_falls_as_h_to_twice_the_degree`). A change that moves one is
/// a change to the elements or the eigen solver, and must refreeze it on purpose.
const FROZEN: [[f64; 6]; 2] = [
    [27.028759694956854, 46.74075182797388, 53.864836466725954, 55.20300406406444, 65.66686752107707, 67.22310359797476],
    [27.45171229550246, 49.48044527307016, 57.283205382145034, 57.29289446123316, 67.36601280397858, 67.43677661134987],
];

/// How far a frozen k² may move: twice the eigen solver's bound. [`modes`] accepts a value θ of
/// the shift-inverted problem at a residual below 1e-10 θ, so θ is within 1e-10 θ of the discrete
/// eigenvalue (a symmetric problem's residual bounds its error), and k² = σ + 1/θ (σ > 0) within
/// as much relatively; the frozen value and a new one may each be that far, on opposite sides.
const FROZEN_TOLERANCE: f64 = 2e-10;

#[test]
fn a_rectangular_cavity_has_its_resonances_and_no_spurious_modes() {
    let exact = resonances(6);
    // (m, n, l) = (1, 0, 1) first, then (1, 1, 0), (0, 1, 1) and (2, 0, 1), and (1, 1, 1) twice:
    // k² ≈ 27.4, 49.3, 57.0 (twice), 66.9 (twice). Against the closed form they converge at the
    // rate of the next test; here they are frozen ([`FROZEN`]).
    let sigma = 0.5 * exact[0];
    for (order, frozen) in [(Order::One, FROZEN[0]), (Order::Two, FROZEN[1])] {
        let c = cavity(order, [4, 2, 3]);
        let (negative, k2) = modes(&c, sigma, exact.len());
        assert_eq!(negative, c.null, "{order:?}: the null space and nothing else lies below σ");
        let moved = k2.iter().zip(frozen).map(|(k, f)| ((k - f) / f).abs()).fold(f64::NAN, f64::max);
        assert!(k2.len() == frozen.len() && moved < FROZEN_TOLERANCE, "{order:?}: {k2:?}, frozen {frozen:?}, error {}", error(&k2, &exact));
    }
}

/// Heavy: the finer mesh of the study takes most of the file's time. The trunk gate keeps the
/// frozen values of the coarse mesh (`a_rectangular_cavity_has_its_resonances_and_no_spurious_modes`).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_eigenvalue_error_falls_as_h_to_twice_the_degree() {
    let exact = resonances(6);
    let sigma = 0.5 * exact[0];
    for (order, rate) in [(Order::One, 2.0), (Order::Two, 4.0)] {
        let errors: Vec<f64> = [[4, 2, 3], [8, 4, 6]]
            .into_iter()
            .map(|cells| {
                let c = cavity(order, cells);
                let (negative, k2) = modes(&c, sigma, exact.len());
                assert_eq!(negative, c.null);
                error(&k2, &exact)
            })
            .collect();
        // Halving h divides the error by 2^rate in the limit. Elements one degree off converge at
        // 2(k ± 1) = 2k ± 2, so within 1 of it, halfway to theirs, the observed rate is told
        // from theirs.
        let observed = (errors[0] / errors[1]).log2();
        assert!((observed - rate).abs() < 1.0, "{order:?}: errors {errors:?}, rate {observed}");
    }
}
