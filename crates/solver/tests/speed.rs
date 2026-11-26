//! Speed of a linear static solve on ordinary solid meshes: a block of quadratic tetrahedra and
//! one of quadratic hexahedra, held at the base and pushed at the top, solved on one thread, on
//! four and on the machine's.
//!
//! The times are printed, not asserted: they are the machine's (run with `--nocapture`, on the
//! commits to compare). What is asserted is what must hold on every machine: the reactions
//! balance the load, and the result is the same bit for bit on every number of threads.

mod common;

use common::*;
use fieldcraft_model::{Distribution, ElementKind, Load, Procedure, Target};
use fieldcraft_solver::{Limits, StaticResult, solve_static};

/// A unit block of `n` cells a side, held at z = 0, with a force along x on its top.
fn block(n: usize, kind: ElementKind) -> Grid {
    let mut g = grid([n; 3], kind, |p| p);
    g.define(iso(210e9, 0.3), None, full(), Procedure::Static {});
    g.fixed("base", |p| near(p[2], 0.0));
    g.node_set("top", |p| near(p[2], 1.0));
    g.load("push", Load::Force { target: Target::Set("top".into()), vector: [1000.0, 0.0, 0.0], distribution: Distribution::Total });
    g
}

fn bits(r: &StaticResult) -> Vec<u64> {
    r.displacement.iter().chain(&r.reaction).flatten().chain(&r.von_mises).map(|x| x.to_bits()).collect()
}

/// Solve on each number of threads, print the times, and check the balance and the bits.
fn timed(name: &str, n: usize, kind: ElementKind, threads: &[usize]) {
    let g = block(n, kind);
    let mut first: Option<Vec<u64>> = None;
    for &count in threads {
        let start = std::time::Instant::now();
        let r = fieldcraft_parallel::with_threads(count, || solve_static(&g.model, None, &Limits::default())).unwrap();
        let seconds = start.elapsed().as_secs_f64();
        eprintln!(
            "static solve, {name}, {n} cells a side: {} unknowns, {} factor entries, {count} threads, {seconds:.2} s",
            r.stats.equations, r.stats.factor_entries
        );
        for (applied, reaction) in r.applied_total.iter().zip(&r.reaction_total) {
            assert!((applied + reaction).abs() < 1e-6 * 1000.0, "{name}: applied {:?}, reactions {:?}", r.applied_total, r.reaction_total);
        }
        let found = bits(&r);
        match &first {
            Some(expected) => assert!(*expected == found, "{name}: {count} threads give another result"),
            None => first = Some(found),
        }
    }
}

/// A small block of each kind in the gate: the same checks, in a moment.
#[test]
fn static_solve_is_the_same_on_any_number_of_threads() {
    timed("tet10", 4, ElementKind::Tet10, &[1, 3]);
    timed("hex20", 4, ElementKind::Hex20, &[1, 3]);
}

/// The blocks whose times are compared between commits.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn static_solve_speed() {
    let machine = fieldcraft_parallel::default_threads();
    timed("tet10", 14, ElementKind::Tet10, &[1, 4, machine]);
    timed("hex20", 18, ElementKind::Hex20, &[1, 4, machine]);
}
