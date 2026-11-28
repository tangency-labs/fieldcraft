//! Cavity steps against theory (docs/rfc/0001-electromagnetics.md §4): models built in code, a
//! box of Kuhn-split tetrahedra, solved by `solve_cavity`.
//!
//! - **An empty rectangular cavity** with perfectly conducting walls (D. M. Pozar, *Microwave
//!   Engineering*, 4th ed., Wiley 2012, §6.4): k² = (mπ/a)² + (nπ/b)² + (lπ/d)², one mode for each
//!   (m, n, l) with at most one index zero, two when none is. The step's null space is the
//!   interior nodes' potentials (and the free edges' at order 2), no mode is missed, and the
//!   eigenvalue error falls as h² at order 1 and h⁴ at order 2, the rate 2k of Nédélec elements of
//!   degree k (D. Boffi, *Finite element approximation of eigenvalue problems*, Acta Numer.
//!   19:1–120, 2010).
//! - **The same cavity filled** with a material of εᵣ and μᵣ: k₀ = k / √(εᵣ μᵣ), so every
//!   frequency is divided by √(εᵣ μᵣ).
//! - **A cavity half filled with a dielectric** slab across its length: the lowest mode,
//!   E = ŷ sin(πx/a) f(z), with f = sin(β₁z) in the slab (0 < z < t) and A sin(β₂(d − z)) in the
//!   air; E_y and ∂E_y/∂z (H_x) continuous at z = t give the transverse resonance condition
//!   β₁ tan(β₂(d − t)) + β₂ tan(β₁t) = 0, β₁² = εᵣk₀² − (π/a)², β₂² = k₀² − (π/a)² (Pozar ch. 3
//!   and 6; R. F. Harrington, *Time-Harmonic Electromagnetic Fields*, McGraw-Hill 1961, ch. 4).
//! - **Two conductors**: perfect conductors on the walls x = 0 and x = a only, the rest magnetic
//!   walls (n × H = 0). φ = x/a is constant on each conductor, so ∇φ = x̂/a is a field with
//!   k = 0 that is no gradient of a potential vanishing on the conductors: the step counts it in
//!   the null space (one potential for the second conductor, RFC §8.4) and returns no mode at
//!   zero. The first mode is E = x̂ cos(πz/d), k² = (π/d)².
//! - **Two separate cavities**, the box with its slab 0.5 < x < 0.75 left out: each part's
//!   constant potential is left out of the null space on its own. With the box's walls
//!   conducting, the first mode is the part 0 < x < 0.5's, conducting at x = 0 and magnetic at
//!   x = 0.5: E = ŷ sin(πx) sin(πz/d), k² = π² + (π/d)².
//! - **Copper walls** (a conductor of σ on every wall): TE₁₀₁'s quality factor
//!   Q = (k a d)³ b η / (2π² Rₛ) / (2a³b + 2bd³ + a³d + ad³), Rₛ = √(ω μ₀ / (2σ)) (Pozar §6.7,
//!   eq. 6.46, and eq. 1.131); Q grows as √σ and falls as 1/√μ_c of the metal; the walls' losses
//!   add, 1/Q = Σ 1/Q_wall; filling the cavity with εᵣ and μᵣ scales Q by μᵣ^¾ εᵣ^-¼ (the same
//!   formula with the medium's k, η = η₀ √(μᵣ/εᵣ) and Rₛ at the lower frequency).
//! - **A cube** has degenerate modes (on its mesh of 2 × 2 × 2 cells at order 2, the second and
//!   third agree to round-off): asking for the second alone is no missed mode.

mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

use common::*;
use fieldcraft_fem::nedelec::{EdgeTet, Order, Potential, declare, gradient, slots};
use fieldcraft_model::*;
use fieldcraft_solver::analysis::TOLERANCE;
use fieldcraft_solver::cavity::{SPEED_OF_LIGHT, VACUUM_PERMEABILITY};
use fieldcraft_solver::sparse::{NodeGraph, add_block};
use fieldcraft_solver::{CavityResult, Limits, solve_cavity};

/// The cavity's sides along x, y and z (metres).
const SIDES: [f64; 3] = [1.0, 0.5, 0.75];

/// A box of `cells` cubes per axis, each six tetrahedra, with sides `sides`; the material of each
/// element by its centroid (`material(centroid)` is (εᵣ, μᵣ)); a cavity step of `order` asking
/// for `modes`; perfect conductors on the walls where `wall(axis, side)` holds.
fn cavity(
    cells: [usize; 3],
    sides: [f64; 3],
    order: u8,
    modes: u32,
    material: impl Fn([f64; 3]) -> (f64, f64),
    wall: impl Fn(usize, usize) -> bool,
) -> Grid {
    let mut g = grid(cells, ElementKind::Tet4, |p| [0, 1, 2].map(|k| sides[k] * p[k]));
    sections(&mut g, material);
    g.model.add_step(Step::new("Step", Procedure::Cavity { modes, min_hz: None, max_hz: None, order }), None).unwrap();
    // Each wall's faces: those whose nodes all lie in its plane (a face lies in one plane).
    let mut faces = Vec::new();
    for (a, side) in (0..3).flat_map(|a| [(a, 0), (a, 1)]) {
        if wall(a, side) {
            let at = side as f64;
            faces.extend(g.faces_where(|p| p[a] == at));
        }
    }
    g.model.add_set("walls", Set { faces, ..Set::default() }).unwrap();
    g.bc("walls", BoundaryCondition::PerfectConductor { target: Target::Set("walls".into()) });
    g
}

/// A section per material, the material of each element by its centroid ((εᵣ, μᵣ)).
fn sections(g: &mut Grid, material: impl Fn([f64; 3]) -> (f64, f64)) {
    let mut regions: BTreeMap<(u64, u64), Vec<u32>> = BTreeMap::new();
    for (label, e) in g.model.elements() {
        let c = e.nodes.iter().fold([0.0; 3], |a, n| {
            let p = g.model.nodes()[n].xyz;
            [0, 1, 2].map(|k| a[k] + p[k] / 4.0)
        });
        let (eps, mu) = material(c);
        regions.entry((eps.to_bits(), mu.to_bits())).or_default().push(*label);
    }
    for (i, ((eps, mu), elements)) in regions.into_iter().enumerate() {
        let (set, name) = (format!("region{i}"), format!("medium{i}"));
        let em = Electromagnetic { permittivity: f64::from_bits(eps), permeability: f64::from_bits(mu) };
        g.model.add_set(&set, Set { elements, ..Set::default() }).unwrap();
        g.model.add_material(&name, Material { electromagnetic: Some(em), ..Material::default() }).unwrap();
        g.model.add_section(&name, Section { set, material: name.clone(), kind: SectionKind::Solid {}, formulation: full() }).unwrap();
    }
}

fn every_wall(_: usize, _: usize) -> bool {
    true
}

fn vacuum(_: [f64; 3]) -> (f64, f64) {
    (1.0, 1.0)
}

fn solve(g: &Grid) -> CavityResult {
    solve_cavity(&g.model, None, &Limits::default()).unwrap()
}

/// k² of a frequency.
fn k2(hz: f64) -> f64 {
    let k = 2.0 * PI * hz / SPEED_OF_LIGHT;
    k * k
}

/// The closed form's k² of an empty box, ascending, each as often as its multiplicity.
fn resonances(sides: [f64; 3], count: usize) -> Vec<f64> {
    let mut out = Vec::new();
    for m in 0..6usize {
        for n in 0..6usize {
            for l in 0..6usize {
                let zeros = [m, n, l].iter().filter(|i| **i == 0).count();
                if zeros > 1 {
                    continue;
                }
                let k2: f64 = [m, n, l].iter().zip(sides).map(|(i, s)| (*i as f64 * PI / s) * (*i as f64 * PI / s)).sum();
                out.extend(std::iter::repeat_n(k2, if zeros == 0 { 2 } else { 1 }));
            }
        }
    }
    out.sort_by(f64::total_cmp);
    out.truncate(count);
    out
}

/// The largest relative error of the computed k² against `exact`.
fn error(r: &CavityResult, exact: &[f64]) -> f64 {
    let k: Vec<f64> = r.modes.iter().map(|m| k2(m.frequency)).collect();
    assert_eq!(k.len(), exact.len(), "{k:?}");
    let errors: Vec<f64> = k.iter().zip(exact).map(|(c, e)| ((c - e) / e).abs()).collect();
    // f64::max passes over a NaN: every error must be a number first.
    assert!(errors.iter().all(|e| e.is_finite()), "{k:?}");
    errors.into_iter().fold(0.0, f64::max)
}

/// The interior nodes, and the edges not in a wall, of a box of `cells` with every wall a
/// conductor.
fn interior(g: &Grid) -> (usize, usize) {
    let inside = |p: [f64; 3]| p.iter().all(|c| *c > 0.0 && *c < 1.0);
    let nodes = g.param.values().filter(|p| inside(**p)).count();
    let mut edges = BTreeSet::new();
    for e in g.model.elements().values() {
        for (i, a) in e.nodes.iter().enumerate() {
            for b in &e.nodes[i + 1..] {
                edges.insert((*a.min(b), *a.max(b)));
            }
        }
    }
    let in_wall = |a: u32, b: u32| (0..3).any(|k| [0.0, 1.0].iter().any(|at| g.param[&a][k] == *at && g.param[&b][k] == *at));
    (nodes, edges.iter().filter(|(a, b)| !in_wall(*a, *b)).count())
}

/// The shift `solve_cavity` takes for a box of `sides` with no lower bound, its materials' 1/μᵣ
/// at least `nu` and εᵣ at most `eps`: σ = −(π/D)² · ν / ε, D the box's diagonal.
fn shift(sides: [f64; 3], nu: f64, eps: f64) -> f64 {
    -PI * PI / sides.iter().map(|s| s * s).sum::<f64>() * nu / eps
}

/// How far the eigen solver may leave a k² from its discrete value, relative to it: it accepts θ
/// = 1/(k² − σ) at a residual below `TOLERANCE` · θ, which bounds θ's error (a symmetric
/// problem), so k² is within `TOLERANCE` · (k² − σ).
fn solver_bound(k2: f64, sigma: f64) -> f64 {
    TOLERANCE * (k2 - sigma) / k2
}

/// The largest relative error against the closed form of the frozen k² of [`SIDES`] on its
/// `[4, 2, 3]` mesh at `order` (5.5 % at order 1, 0.81 % at order 2). The error of a mode falls
/// with (k h)^2k, so on that mesh a mode no higher than the sixth, of a cavity whose cells are no
/// larger, is that close to its closed form: the bound of the tests below on that mesh.
fn mesh_error(order: u8) -> f64 {
    let exact = resonances(SIDES, 6);
    let frozen = FROZEN[usize::from(order - 1)];
    frozen.iter().zip(&exact).map(|(k, e)| ((k - e) / e).abs()).fold(0.0, f64::max)
}

/// The six lowest k² of [`SIDES`] on its mesh of `[4, 2, 3]` cells, at orders 1 and 2: our
/// result, frozen once it matched the closed form (largest errors 5.5 % and 0.81 %, falling at
/// the rates of `the_eigenvalue_error_falls_as_h_to_twice_the_degree`). A change that moves one is
/// a change to the elements or the eigen solver, and must refreeze it on purpose.
const FROZEN: [[f64; 6]; 2] = [
    [27.028759694956854, 46.74075182797379, 53.8648364667258, 55.203004064064295, 65.66686752107692, 67.22310359797497],
    [27.451712295502478, 49.48044527307019, 57.28320538214481, 57.292894461233175, 67.36601280397865, 67.43677661134991],
];

/// How far a frozen k² may move: twice the eigen solver's bound ([`solver_bound`]) at the step's
/// shift, as the frozen value and a new one may each be that far, on opposite sides. With
/// σ = −(π/D)² below zero the bound is above `TOLERANCE` (about 1.2 `TOLERANCE` for the first k²).
fn frozen_tolerance(k2: f64) -> f64 {
    2.0 * solver_bound(k2, shift(SIDES, 1.0, 1.0))
}

#[test]
fn a_rectangular_cavity_has_its_resonances_and_no_spurious_modes() {
    // (m, n, l) = (1, 0, 1) first, then (1, 1, 0), (0, 1, 1) and (2, 0, 1), and (1, 1, 1) twice:
    // k² ≈ 27.4, 49.3, 57.0 (twice), 66.9 (twice). Against the closed form they converge at the
    // rate of the next test; here they are frozen ([`FROZEN`]).
    let exact = resonances(SIDES, 6);
    for (order, frozen) in [(1, FROZEN[0]), (2, FROZEN[1])] {
        let g = cavity([4, 2, 3], SIDES, order, 6, vacuum, every_wall);
        let r = solve(&g);
        // No warning: in particular none that a mode below the highest one was missed.
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
        let k: Vec<f64> = r.modes.iter().map(|m| k2(m.frequency)).collect();
        let within = k.len() == frozen.len() && k.iter().zip(frozen).all(|(k, f)| ((k - f) / f).abs() < frozen_tolerance(f));
        assert!(within, "order {order}: {k:?}, frozen {frozen:?}, error {}", error(&r, &exact));
        let (nodes, edges) = interior(&g);
        assert_eq!(r.null_space, if order == 1 { nodes } else { nodes + edges }, "order {order}");
        // TE101 is E = ŷ sin(πx/a) sin(πz/d): the y component carries its energy, the field
        // follows the closed form's at the nodes, and it is scaled to a largest |E| of 1.
        let field = &r.modes[0].field;
        let energy = |axis: usize| field.iter().map(|f| f[axis] * f[axis]).sum::<f64>();
        let share = energy(1) / (energy(0) + energy(1) + energy(2));
        let largest = field.iter().map(|f| f.iter().map(|c| c * c).sum::<f64>()).fold(0.0, f64::max);
        // Against the closed form's field at the nodes: their cosine.
        let nodes: Vec<[f64; 3]> = r.nodes.iter().map(|n| g.model.nodes()[n].xyz).collect();
        let exact_field: Vec<f64> = nodes.iter().map(|p| (PI * p[0] / SIDES[0]).sin() * (PI * p[2] / SIDES[2]).sin()).collect();
        let dot: f64 = field.iter().zip(&exact_field).map(|(f, e)| f[1] * e).sum();
        let squared: f64 = field.iter().map(|f| f.iter().map(|c| c * c).sum::<f64>()).sum();
        let cosine = dot / (squared.sqrt() * exact_field.iter().map(|e| e * e).sum::<f64>().sqrt());
        // The modes are orthogonal, so cosine² is TE101's share of the field: above ½, the field
        // is mostly TE101's and no other mode's, and its y component carries most of its energy.
        // (0.92 and 0.93 at order 1, Whitney functions averaged at the corners; 0.997 at order 2.)
        let half = 0.5f64;
        assert!(share > half && cosine > half.sqrt() && (largest - 1.0).abs() < 1e-12, "order {order}: {share}, {cosine}, {largest}");
    }
}

/// Heavy: the finer mesh takes most of the file's time. The trunk gate keeps the frozen values
/// of the coarse one (`a_rectangular_cavity_has_its_resonances_and_no_spurious_modes`).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_eigenvalue_error_falls_as_h_to_twice_the_degree() {
    let exact = resonances(SIDES, 6);
    for (order, rate) in [(1, 2.0), (2, 4.0)] {
        let errors: Vec<f64> =
            [[4, 2, 3], [8, 4, 6]].into_iter().map(|cells| error(&solve(&cavity(cells, SIDES, order, 6, vacuum, every_wall)), &exact)).collect();
        // Halving h divides the error by 2^rate in the limit. Elements one degree off converge
        // at 2(k ± 1) = 2k ± 2, so within 1 of it, halfway to theirs, the observed rate is told
        // from theirs.
        let observed = (errors[0] / errors[1]).log2();
        assert!((observed - rate).abs() < 1.0, "order {order}: errors {errors:?}, rate {observed}");
    }
}

/// Every column of the discrete gradient is in the kernel of the curl–curl matrix assembled
/// over the whole mesh, before the walls are removed: one per node, and at order 2 one per edge.
/// A column has a few entries and K is symmetric, so K g is the sum of K's rows at them. Each
/// entry of K g sums a few dozen products that cancel exactly, so it is zero to a few dozen
/// rounding errors of K's largest entry (dozens of ε, about 1e-14): 1e-12 of it is well above.
#[test]
fn the_discrete_gradient_is_in_the_kernel_of_the_curl_curl_matrix() {
    let g = grid([2, 1, 2], ElementKind::Tet4, |p| [0, 1, 2].map(|k| SIDES[k] * p[k]));
    let tets: Vec<[u32; 4]> = g.model.elements().values().map(|e| [e.nodes[0], e.nodes[1], e.nodes[2], e.nodes[3]]).collect();
    let xyz = |l: u32| g.model.nodes()[&l].xyz;
    for order in [Order::One, Order::Two] {
        let mut b = fieldcraft_fem::dof::DofLayout::builder();
        for t in &tets {
            declare(&mut b, *t, order).unwrap();
        }
        let layout = b.build();
        let element_slots: Vec<Vec<usize>> = tets.iter().map(|t| slots(&layout, *t, order).unwrap()).collect();
        let graph = NodeGraph::new(layout.slots(), element_slots.iter());
        let mut k = graph.zero_matrix(&layout);
        for (t, s) in tets.iter().zip(&element_slots) {
            let (e, _) = EdgeTet::of_labels(order, t.map(xyz), *t).unwrap();
            add_block(&mut k, &layout, s, order.kinds(), &e.curl_curl([1.0; 3]).unwrap().data).unwrap();
        }
        let mut columns: BTreeMap<Potential, Vec<(usize, f64)>> = BTreeMap::new();
        for (r, c, v) in gradient(&layout) {
            columns.entry(c).or_default().push((r, v));
        }
        let edges = columns.keys().filter(|c| matches!(c, Potential::Edge(_))).count();
        assert_eq!((columns.len() - edges, edges), (g.model.nodes().len(), if order == Order::Two { layout.edges() } else { 0 }));
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
    }
}

#[test]
fn a_filled_cavity_resonates_lower_by_the_root_of_its_permittivity_and_permeability() {
    let empty = solve(&cavity([4, 2, 3], SIDES, 2, 4, vacuum, every_wall));
    let filled = solve(&cavity([4, 2, 3], SIDES, 2, 4, |_| (4.0, 2.25), every_wall));
    assert_eq!((empty.modes.len(), filled.modes.len()), (4, 4));
    // K is divided by 2.25 and M multiplied by 4, so the discrete k₀² by 9 exactly; each computed
    // k₀² is within the solver's bound of its own (the shift scales with them), and a frequency
    // within half that, so their ratio within the sum.
    let sigma = shift(SIDES, 1.0, 1.0);
    for (e, f) in empty.modes.iter().zip(&filled.modes) {
        let ratio = e.frequency / f.frequency;
        let bound = solver_bound(k2(e.frequency), sigma);
        assert!((ratio / 3.0 - 1.0).abs() < bound, "{ratio}, bound {bound}");
    }
}

#[test]
fn a_half_filled_cavity_follows_its_transverse_resonance() {
    // a × b × d = 1 × 0.5 × 1, a slab of εᵣ = 2.25 over 0 < z < 0.5.
    let (sides, t, eps) = ([1.0, 0.5, 1.0], 0.5, 2.25);
    let kc2 = (PI / sides[0]) * (PI / sides[0]);
    // The lowest root in k₀², by bisection on the condition written with sines and cosines, so
    // that it has no poles: β₁ sin(β₂ l₂) cos(β₁ t) + β₂ sin(β₁ t) cos(β₂ l₂).
    let condition = |k0sq: f64| {
        let (b1, b2) = ((eps * k0sq - kc2).sqrt(), (k0sq - kc2).sqrt());
        let l2 = sides[2] - t;
        b1 * (b2 * l2).sin() * (b1 * t).cos() + b2 * (b1 * t).sin() * (b2 * l2).cos()
    };
    let (mut lo, mut hi) = (kc2 * 1.000_001, kc2 + (PI / sides[2]) * (PI / sides[2]));
    assert!(condition(lo) * condition(hi) < 0.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if condition(lo) * condition(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let exact = 0.5 * (lo + hi);
    let r = solve(&cavity([4, 2, 4], sides, 2, 1, |c| if c[2] < t { (eps, 1.0) } else { (1.0, 1.0) }, every_wall));
    let e = (k2(r.modes[0].frequency) - exact) / exact;
    // Its cells (0.25 × 0.25 × 0.25) are no larger than those of the frozen mesh, but in the slab
    // the wavenumber is √εᵣ k₀, and the error grows with (k h)⁴ at order 2: εᵣ² times the frozen
    // mesh's error bounds it (the mode is below the frozen sixth, k₀² < 67).
    let bound = eps * eps * mesh_error(2);
    assert!(exact < FROZEN[1][5] && e.abs() < bound, "k₀² {} against {exact}: {e}, bound {bound}", k2(r.modes[0].frequency));
}

#[test]
fn a_second_conductor_adds_a_potential_and_no_mode_at_zero() {
    for order in [1, 2] {
        // Conductors on x = 0 and x = a; then on x = 0 alone.
        let two = solve(&cavity([4, 2, 3], SIDES, order, 2, vacuum, |axis, _| axis == 0));
        let one = solve(&cavity([4, 2, 3], SIDES, order, 2, vacuum, |axis, side| axis == 0 && side == 0));
        // The second plate holds its 3 × 4 nodes (and at order 2 its edges), which were free,
        // and adds one potential.
        assert!(one.null_space > two.null_space, "order {order}: {} and {}", one.null_space, two.null_space);
        if order == 1 {
            assert_eq!(one.null_space - two.null_space, 3 * 4 - 1);
        }
        // On the frozen mesh and below its sixth k²: within its error ([`mesh_error`]).
        let exact = (PI / SIDES[2]) * (PI / SIDES[2]);
        let first = k2(two.modes[0].frequency);
        assert!(((first - exact) / exact).abs() < mesh_error(order), "order {order}: {first} against {exact}");
        assert!(two.warnings.is_empty(), "{:?}", two.warnings);
    }
}

/// The box of [`cavity`] (`[4, 2, 3]` cells, order `order`) with the slab 0.5 < x < 0.75
/// removed, in two parts; the walls of the box conducting when `walls`.
fn separate(order: u8, walls: bool) -> Grid {
    let mut g = grid([4, 2, 3], ElementKind::Tet4, |p| [0, 1, 2].map(|k| SIDES[k] * p[k]));
    let centroid_x = |g: &Grid, e: &Element| e.nodes.iter().map(|n| g.model.nodes()[n].xyz[0]).sum::<f64>() / 4.0;
    let gap: Vec<u32> = g.model.elements().iter().filter(|(_, e)| (0.5..0.75).contains(&centroid_x(&g, e))).map(|(l, _)| *l).collect();
    for l in gap {
        g.model.remove_element(l).unwrap();
    }
    let kept: Vec<u32> = g.model.elements().keys().copied().collect();
    let em = Electromagnetic { permittivity: 1.0, permeability: 1.0 };
    g.model.add_set("parts", Set { elements: kept.clone(), ..Set::default() }).unwrap();
    g.model.add_material("vacuum", Material { electromagnetic: Some(em), ..Material::default() }).unwrap();
    g.model
        .add_section("vacuum", Section { set: "parts".into(), material: "vacuum".into(), kind: SectionKind::Solid {}, formulation: full() })
        .unwrap();
    g.model.add_step(Step::new("Step", Procedure::Cavity { modes: 2, min_hz: None, max_hz: None, order }), None).unwrap();
    if walls {
        let faces = (0..3)
            .flat_map(|a| [(a, 0.0), (a, 1.0)])
            .flat_map(|(a, at)| g.faces_where(|p| p[a] == at))
            .filter(|f| kept.contains(&f.element))
            .collect();
        g.model.add_set("walls", Set { faces, ..Set::default() }).unwrap();
        g.bc("walls", BoundaryCondition::PerfectConductor { target: Target::Set("walls".into()) });
    }
    g
}

#[test]
fn separate_cavities_each_leave_out_their_own_constant_potential() {
    for order in [1, 2] {
        // Conducting walls: the far part's walls are a second conductor, and its nodes a part
        // with no first conductor.
        // The frozen mesh's cells, below its sixth k²: within its error ([`mesh_error`]).
        let r = solve(&separate(order, true));
        let exact = PI * PI + (PI / SIDES[2]) * (PI / SIDES[2]);
        let first = k2(r.modes[0].frequency);
        assert!(((first - exact) / exact).abs() < mesh_error(order), "order {order}: {first} against {exact}");
        assert!(r.warnings.is_empty(), "order {order}: {:?}", r.warnings);
        // No conductor at all: one potential left out in each part, and no mode at zero. Every
        // wall magnetic, the parts' modes are those of conducting boxes (E and H exchanged): the
        // lowest of each part, 0.5 × 0.5 × 0.75 and 0.25 × 0.5 × 0.75, is (0, 1, 1), k² =
        // (π/0.5)² + (π/0.75)², and it holds in both parts and in the first along x too.
        let r = solve(&separate(order, false));
        assert_eq!(r.modes.len(), 2, "order {order}");
        let exact = (PI / 0.5) * (PI / 0.5) + (PI / SIDES[2]) * (PI / SIDES[2]);
        let k: Vec<f64> = r.modes.iter().map(|m| k2(m.frequency)).collect();
        assert!(k.iter().all(|k| ((k - exact) / exact).abs() < mesh_error(order)), "order {order}: {k:?} against {exact}");
    }
}

#[test]
fn asking_for_one_of_a_degenerate_cluster_is_no_missed_mode() {
    let r = solve(&cavity([2, 2, 2], [1.0; 3], 2, 2, vacuum, every_wall));
    assert_eq!(r.modes.len(), 2);
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    // The third, not asked for, is the second's twin.
    let three = solve(&cavity([2, 2, 2], [1.0; 3], 2, 3, vacuum, every_wall));
    // Twins agree to the solver's bound on each (a frequency to half its k²'s), so to its sum.
    let (second, third) = (three.modes[1].frequency, three.modes[2].frequency);
    let bound = solver_bound(k2(second), shift([1.0; 3], 1.0, 1.0));
    assert!((third - second).abs() < bound * second, "{second} and {third}");
}

#[test]
fn a_lower_bound_of_zero_is_no_bound() {
    let mut g = cavity([2, 1, 2], SIDES, 1, 2, vacuum, every_wall);
    let free = solve(&g);
    let mut step = g.model.step("Step").unwrap().clone();
    step.procedure = Procedure::Cavity { modes: 2, min_hz: Some(0.0), max_hz: None, order: 1 };
    g.model.replace_step(step).unwrap();
    let from_zero = solve(&g);
    assert_eq!(from_zero.modes, free.modes);
}

/// With a `min_hz`, the modes between it and the highest one found are counted too.
#[test]
fn a_lower_bound_is_checked_for_missed_modes() {
    let mut g = cavity([4, 2, 3], SIDES, 1, 4, vacuum, every_wall);
    let free = solve(&g);
    assert!(free.warnings.is_empty(), "{:?}", free.warnings);
    let f: Vec<f64> = free.modes.iter().map(|m| m.frequency).collect();
    let mut step = g.model.step("Step").unwrap().clone();
    step.procedure = Procedure::Cavity { modes: 2, min_hz: Some(0.5 * (f[0] + f[1])), max_hz: None, order: 1 };
    g.model.replace_step(step.clone()).unwrap();
    let above = solve(&g);
    assert_eq!(above.modes.len(), 2);
    // Each within the solver's bound at its own shift (a frequency within half its k²'s).
    let bound = |hz: f64, sigma: f64| 0.5 * solver_bound(k2(hz), sigma);
    for (a, b) in above.modes.iter().zip(&free.modes[1..]) {
        let both = bound(b.frequency, shift(SIDES, 1.0, 1.0)) + bound(a.frequency, k2(0.5 * (f[0] + f[1])));
        assert!((a.frequency - b.frequency).abs() < both * b.frequency, "{} and {}", a.frequency, b.frequency);
    }
    assert!(above.warnings.is_empty(), "{:?}", above.warnings);
    // A basis of one Lanczos block stops short of the band's modes, and the eigensolver returns
    // its Ritz values unconverged: the check says modes were missed.
    let short = Limits { max_lanczos_vectors: 6, ..Limits::default() };
    match solve_cavity(&g.model, None, &short) {
        Ok(r) => assert!(r.warnings.iter().any(|w| w.contains("modes lie between")), "{r:?}"),
        Err(e) => assert!(e.to_string().contains("converge"), "{e}"),
    }
}

/// Asking for more modes than the cavity has is reported.
#[test]
fn a_shortfall_is_against_the_modes_asked_for() {
    let mut g = cavity([1, 1, 1], [1.0; 3], 1, 1, vacuum, every_wall);
    let one = solve(&g);
    assert_eq!(one.modes.len(), 1);
    assert!(one.warnings.is_empty(), "{:?}", one.warnings);
    let mut step = g.model.step("Step").unwrap().clone();
    step.procedure = Procedure::Cavity { modes: 2, min_hz: None, max_hz: None, order: 1 };
    g.model.replace_step(step).unwrap();
    let two = solve(&g);
    assert_eq!(two.modes.len(), 1);
    assert_eq!(two.warnings, ["1 of the 2 modes asked for were found"]);
}

#[test]
fn what_a_cavity_step_refuses() {
    // A material without electromagnetic properties.
    let mut g = cavity([1, 1, 1], SIDES, 1, 1, vacuum, every_wall);
    let mut m = g.model.material("medium0").unwrap().clone();
    m.electromagnetic = None;
    m.density = Some(1.0);
    g.model.replace_material("medium0", m).unwrap();
    assert!(solve_cavity(&g.model, None, &Limits::default()).is_err());
    // Hexahedra.
    let mut h = grid([1, 1, 1], ElementKind::Hex8, |p| p);
    sections(&mut h, vacuum);
    h.model.add_step(Step::new("Step", Procedure::Cavity { modes: 1, min_hz: None, max_hz: None, order: 1 }), None).unwrap();
    assert!(solve_cavity(&h.model, None, &Limits::default()).is_err());
    // A range with no mode in it: an error, not an empty result.
    let mut g = cavity([1, 1, 1], SIDES, 1, 1, vacuum, every_wall);
    let mut step = g.model.step("Step").unwrap().clone();
    step.procedure = Procedure::Cavity { modes: 1, min_hz: Some(1e3), max_hz: Some(2e3), order: 1 };
    g.model.replace_step(step).unwrap();
    let err = solve_cavity(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(err.contains("no resonant mode was found between 1000 and 2000 Hz"), "{err}");
    // Steps of order 3 and of no modes are invalid.
    assert!(Procedure::Cavity { modes: 1, min_hz: None, max_hz: None, order: 3 }.validate().is_err());
    assert!(Procedure::Cavity { modes: 0, min_hz: None, max_hz: None, order: 1 }.validate().is_err());
}

/// Copper, S/m (Pozar 2012, appendix F).
const COPPER: f64 = 5.813e7;

/// The box's walls where `wall(axis, side)` holds made a conductor of `conductivity` and `permeability`.
fn conducting(mut g: Grid, conductivity: f64, permeability: f64) -> Grid {
    let bc = BoundaryCondition::Conductor { target: Target::Set("walls".into()), conductivity, permeability };
    g.model.replace_boundary_condition("walls", bc).unwrap();
    g
}

/// The lowest mode's quality factor.
fn quality(g: &Grid) -> f64 {
    let r = solve(g);
    r.modes[0].quality.unwrap_or_else(|| panic!("no quality factor: {:?}", r.warnings))
}

/// Pozar's eq. 6.46: TE₁₀₁ of an empty a × b × d box with walls of σ.
fn te101_quality(sides: [f64; 3], sigma: f64) -> f64 {
    let [a, b, d] = sides;
    let k = PI * (1.0 / (a * a) + 1.0 / (d * d)).sqrt();
    let omega = k * SPEED_OF_LIGHT;
    let rs = (omega * VACUUM_PERMEABILITY / (2.0 * sigma)).sqrt();
    let eta = VACUUM_PERMEABILITY * SPEED_OF_LIGHT;
    (k * a * d).powi(3) * b * eta / (2.0 * PI * PI * rs) / (2.0 * a.powi(3) * b + 2.0 * b * d.powi(3) + a.powi(3) * d + a * d.powi(3))
}

#[test]
fn copper_walls_give_te101_its_quality_factor() {
    // TE₁₀₁ is the lowest mode of 1 × 0.5 × 0.75 (k² = π² (1 + 1/0.75²), then TE₁₁₀ at 5π²).
    // On 8 × 4 × 6 cells at order 1 Q is 1.4 % above the closed form (see the convergence below).
    let exact = te101_quality(SIDES, COPPER);
    let q = quality(&conducting(cavity([8, 4, 6], SIDES, 1, 1, vacuum, every_wall), COPPER, 1.0));
    assert!((q / exact - 1.0).abs() < 0.02, "Q {q} against {exact}");
    // A perfect conductor absorbs nothing: no quality factor.
    assert!(solve(&cavity([2, 1, 2], SIDES, 1, 1, vacuum, every_wall)).modes[0].quality.is_none());
}

/// The quality factor's error converges as h² at both orders: the losses take the curl's trace on
/// the walls. Measured on 4 × 2 × 3 → 8 × 4 × 6: order 1 5.4 % → 1.4 %, order 2 −7.4 % → −2.2 %.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_quality_factor_converges_to_pozars() {
    let exact = te101_quality(SIDES, COPPER);
    for order in [1, 2] {
        let error =
            |n: usize| (quality(&conducting(cavity([2 * n, n, 3 * n / 2], SIDES, order, 1, vacuum, every_wall), COPPER, 1.0)) / exact - 1.0).abs();
        let (coarse, fine) = (error(2), error(4));
        assert!(fine < coarse / 3.0 && fine < 0.025, "order {order}: {coarse} → {fine}");
    }
}

#[test]
fn the_quality_factor_scales_with_the_metal() {
    // Q ∝ 1/Rₛ ∝ √(σ/μ_c), on the same mesh and mode: exact to rounding.
    let q = |sigma, mu| quality(&conducting(cavity([2, 1, 2], SIDES, 1, 1, vacuum, every_wall), sigma, mu));
    let base = q(COPPER, 1.0);
    assert!((q(4.0 * COPPER, 1.0) / base - 2.0).abs() < 1e-9);
    assert!((q(COPPER, 9.0) / base - 1.0 / 3.0).abs() < 1e-9);
}

#[test]
fn the_walls_losses_add() {
    // 1/Q of every wall conducting = the sum over the walls each conducting alone (the rest
    // perfect): the same field, its losses wall by wall.
    let lossy = |pick: &dyn Fn(usize, usize) -> bool| {
        let mut g = cavity([2, 1, 2], SIDES, 1, 1, vacuum, every_wall);
        let mut faces = Vec::new();
        for (a, side) in (0..3).flat_map(|a| [(a, 0), (a, 1)]) {
            if pick(a, side) {
                faces.extend(g.faces_where(|p| p[a] == side as f64));
            }
        }
        g.model.add_set("metal", Set { faces, ..Set::default() }).unwrap();
        g.bc("metal", BoundaryCondition::Conductor { target: Target::Set("metal".into()), conductivity: COPPER, permeability: 1.0 });
        quality(&g)
    };
    let all = lossy(&|_, _| true);
    let sum: f64 = (0..3).flat_map(|a| [(a, 0), (a, 1)]).map(|w| 1.0 / lossy(&|a, s| (a, s) == w)).sum();
    assert!((all * sum - 1.0).abs() < 1e-9, "1/Q {} against Σ {sum}", 1.0 / all);
}

#[test]
fn a_filled_cavitys_quality_factor_scales_with_its_medium() {
    let q = |eps, mu| quality(&conducting(cavity([2, 1, 2], SIDES, 2, 1, |_| (eps, mu), every_wall), COPPER, 1.0));
    let (empty, filled) = (q(1.0, 1.0), q(2.0, 3.0));
    let expected = 3.0f64.powf(0.75) * 2.0f64.powf(-0.25);
    assert!((filled / empty / expected - 1.0).abs() < 1e-9, "{} against {expected}", filled / empty);
}
