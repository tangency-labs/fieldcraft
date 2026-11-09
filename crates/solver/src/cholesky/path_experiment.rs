//! Test-only real-factorization comparison; all variants reside in the same v3 binary.
use super::*;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
#[path = "path_candidate.rs"]
mod candidate;
static MODE: AtomicUsize = AtomicUsize::new(0);
static CALLS: AtomicU64 = AtomicU64::new(0);
static FLOPS: AtomicU64 = AtomicU64::new(0);
pub(super) fn enabled(mu: usize, wu: usize, kw: usize) -> bool {
    let eligible = wu >= PANEL && mu.saturating_mul(wu).saturating_mul(kw) >= PARALLEL_WORK;
    if eligible {
        CALLS.fetch_add(1, Ordering::Relaxed);
        let flops: u64 = (0..mu).step_by(TILE).map(|r| (2 * wu.min(r + TILE.min(mu - r)).div_ceil(NR) * NR * TILE * kw) as u64).sum();
        FLOPS.fetch_add(flops, Ordering::Relaxed);
    }
    eligible && MODE.load(Ordering::Relaxed) != 0
}
pub(super) fn update(
    l: &Panel<'_>,
    cols: std::ops::Range<usize>,
    p: usize,
    mu: usize,
    wu: usize,
    m: usize,
    first: usize,
    relative: &[usize],
    values: &mut [f64],
) -> std::result::Result<(), FactorError> {
    let target = candidate::Target { m, first, relative, values };
    if MODE.load(Ordering::Relaxed) == 2 {
        candidate::update::<true>(l, cols, p, mu, wu, target)
    } else {
        candidate::update::<false>(l, cols, p, mu, wu, target)
    }
}
#[test]
#[ignore = "explicit real Cholesky path experiment; run alone, fixed v3 binary"]
fn real_factor_paths() {
    let child = 384;
    let separator = 768;
    let n = 2 * child + separator;
    let mut a = UpperColumns { n, ptr: vec![0], idx: vec![], val: vec![] };
    for c in 0..n {
        for r in 0..=c {
            if c >= 2 * child || r / child == c / child {
                a.idx.push(r);
                a.val.push(if r == c { n as f64 + 1.0 } else { ((r * 17 + c * 31) % 29) as f64 / 29.0 - 0.5 });
            }
        }
        a.ptr.push(a.idx.len());
    }
    MODE.store(0, Ordering::Relaxed);
    let reference = parallel::with_threads(1, || Factor::new(&a, usize::MAX)).unwrap();
    let mut reference_x = vec![1.0; n];
    reference.solve(&mut reference_x).unwrap();
    let threads = std::env::var("FC_PATH_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or(1);
    let mode = std::env::var("FC_PATH_MODE").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
    assert!(mode <= 2);
    CALLS.store(0, Ordering::Relaxed);
    FLOPS.store(0, Ordering::Relaxed);
    MODE.store(mode, Ordering::Relaxed);
    let clock = std::time::Instant::now();
    let factor = parallel::with_threads(threads, || Factor::new(&a, usize::MAX)).unwrap();
    let factor_seconds = clock.elapsed().as_secs_f64();
    let mut x = vec![1.0; n];
    factor.solve(&mut x).unwrap();
    let command_seconds = clock.elapsed().as_secs_f64();
    assert_eq!(factor.perm, reference.perm);
    assert!(factor.values.iter().flatten().zip(reference.values.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
    assert!(x.iter().zip(&reference_x).all(|(a, b)| a.to_bits() == b.to_bits()));
    let hash = x.iter().fold(0xcbf29ce484222325u64, |h, v| (h ^ v.to_bits()).wrapping_mul(0x100000001b3));
    eprintln!(
        "PATH_RESULT mode={mode} threads={threads} n={n} factor_seconds={factor_seconds:.9} factor_solve_seconds={command_seconds:.9} calls={} padded_flops={} bitwise=true relative_error=0 hash={hash:016x}",
        CALLS.load(Ordering::Relaxed),
        FLOPS.load(Ordering::Relaxed)
    );
    MODE.store(0, Ordering::Relaxed);
}
