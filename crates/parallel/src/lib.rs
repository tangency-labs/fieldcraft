//! Deterministic data parallelism for FieldCraft's CPU solvers.
//!
//! The solvers promise the same result bit for bit on every run, platform and thread count
//! (AGENTS.md, Determinism). Threads keep that promise when every floating-point value is
//! computed by exactly one task in an order that does not depend on how many threads run them:
//!
//! - **Fixed partitions.** Work is split into blocks whose sizes are constants of the algorithm
//!   (never derived from the number of threads), and each block is computed by one task.
//! - **Ordered results.** [`map`] returns the blocks' results in block order; whatever combines
//!   them (sums, minima, scatters into a global vector) does so on the calling thread in that
//!   order, or in the order the sequential algorithm used.
//! - **No shared accumulators.** Tasks write only into the disjoint slices they are handed
//!   ([`for_each_mut`], [`chunks_mut`]); there are no atomics on floating-point values and no
//!   reductions whose shape depends on scheduling (work stealing decides only *when* a block
//!   runs, never *what* it computes).
//!
//! Under those rules a run on one thread and a run on many execute the same operations on the
//! same operands, so the results are identical; the browser build (`wasm32`) runs the same
//! blocks one after another on the calling thread.
//!
//! **Thread count.** [`with_threads`] runs a closure on a pool of that many threads (pools are
//! cached when idle; concurrent commands lease separate pools, each with at most
//! [`MAX_THREADS`]). Nested work keeps its command's pool. Outside it, the helpers use the default pool:
//! [`ENV`] (`FIELDCRAFT_THREADS`) when it holds a whole number from 1 to [`MAX_THREADS`], else
//! the available cores.
//!
//! Rayon (MIT OR Apache-2.0) provides the pool (docs/adr/0001-libraries.md).
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

/// Most threads a solve may use.
pub const MAX_THREADS: usize = 256;

/// The environment variable that sets the default thread count.
pub const ENV: &str = "FIELDCRAFT_THREADS";

/// A thread count from text: a whole number from 1 to [`MAX_THREADS`].
pub fn parse_threads(text: &str) -> Option<usize> {
    text.trim().parse::<usize>().ok().filter(|n| (1..=MAX_THREADS).contains(n))
}

/// The cores this process may use (at most [`MAX_THREADS`]; 1 in the browser).
pub fn available() -> usize {
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::thread::available_parallelism().map_or(1, |n| n.get()).clamp(1, MAX_THREADS)
    }
    #[cfg(target_arch = "wasm32")]
    {
        1
    }
}

/// The default thread count: [`ENV`] when it is valid, else [`available`].
pub fn default_threads() -> usize {
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::env::var(ENV).ok().and_then(|v| parse_threads(&v)).unwrap_or_else(available)
    }
    #[cfg(target_arch = "wasm32")]
    {
        1
    }
}

/// Run `f` with `threads` worker threads (clamped to 1..=[`MAX_THREADS`]) serving the helpers of
/// this crate. The result does not depend on `threads`; only the time does.
pub fn with_threads<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
    pool::with_threads(threads.clamp(1, MAX_THREADS), f)
}

/// The number of threads the helpers would use here.
pub fn current_threads() -> usize {
    pool::current_threads()
}

/// `f(i)` for `i` in `0..count`, run in parallel, results in index order.
pub fn map<T: Send>(count: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    pool::map(count, f)
}

/// `f(i)` for `i` in `0..count`, results in index order, run in parallel in blocks of `grain`
/// consecutive indices (one task per block: for cheap `f`).
pub fn map_grained<T: Send>(count: usize, grain: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    let ranges: Vec<std::ops::Range<usize>> = blocks(count, grain).collect();
    let parts = map(ranges.len(), |k| ranges.get(k).cloned().unwrap_or(0..0).map(&f).collect::<Vec<T>>());
    parts.into_iter().flatten().collect()
}

/// `f(i, item)` for every item, run in parallel; each task gets its own item mutably (use it to
/// hand tasks disjoint parts of the output).
pub fn for_each_mut<T: Send>(items: &mut [T], f: impl Fn(usize, &mut T) + Sync + Send) {
    pool::for_each_mut(items, f);
}

/// `f(k, chunk)` for the consecutive chunks of `chunk` items of `data` (the last may be
/// shorter), run in parallel.
pub fn chunks_mut<T: Send>(data: &mut [T], chunk: usize, f: impl Fn(usize, &mut [T]) + Sync + Send) {
    let mut parts: Vec<&mut [T]> = data.chunks_mut(chunk.max(1)).collect();
    for_each_mut(&mut parts, |k, part| f(k, part));
}

/// The ranges of `len` items in blocks of `block` (the last may be shorter).
pub fn blocks(len: usize, block: usize) -> impl Iterator<Item = std::ops::Range<usize>> + Clone {
    let block = block.max(1);
    (0..len.div_ceil(block)).map(move |k| k * block..((k + 1) * block).min(len))
}

#[cfg(not(target_arch = "wasm32"))]
mod pool {
    use std::sync::{Arc, Mutex, OnceLock};

    use rayon::prelude::*;
    use rayon::{ThreadPool, ThreadPoolBuilder};

    /// Pools kept for reuse, most recently used last; older ones are dropped (their threads
    /// end) so a caller cycling through thread counts cannot pile up threads.
    const CACHED: usize = 4;

    /// Worker stack size: the element kernels keep small matrices on the stack.
    const STACK: usize = 8 << 20;

    static POOLS: Mutex<Vec<(usize, Arc<ThreadPool>)>> = Mutex::new(Vec::new());
    static DEFAULT: OnceLock<Option<Arc<ThreadPool>>> = OnceLock::new();

    fn build(threads: usize) -> Option<Arc<ThreadPool>> {
        ThreadPoolBuilder::new().num_threads(threads).stack_size(STACK).thread_name(|i| format!("fieldcraft-worker-{i}")).build().ok().map(Arc::new)
    }

    fn cached(threads: usize) -> Option<Arc<ThreadPool>> {
        let mut pools = POOLS.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        // A running command owns an Arc in addition to the cache's. Sharing that pool
        // lets a blocked command monopolize another command's workers (including the
        // worker that must observe cancellation). Reuse only idle pools.
        if let Some(at) = pools.iter().position(|(n, pool)| *n == threads && Arc::strong_count(pool) == 1) {
            let entry = pools.remove(at);
            let pool = entry.1.clone();
            pools.push(entry);
            return Some(pool);
        }
        let pool = build(threads)?;
        if pools.len() >= CACHED {
            pools.remove(0);
        }
        pools.push((threads, pool.clone()));
        Some(pool)
    }

    fn default_pool() -> Option<Arc<ThreadPool>> {
        DEFAULT.get_or_init(|| build(super::default_threads())).clone()
    }

    pub fn with_threads<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
        // Nested work of this command keeps its pool; independent callers lease idle pools.
        if rayon::current_thread_index().is_some() && rayon::current_num_threads() == threads {
            return f();
        }
        match cached(threads) {
            Some(pool) => pool.install(f),
            // No threads to be had: the helpers fall back to the default pool or to the calling
            // thread, with the same result.
            None => f(),
        }
    }

    /// Run `op` inside a pool: the current one, else the default one, else on this thread.
    fn inside<R: Send>(op: impl FnOnce() -> R + Send) -> R {
        if rayon::current_thread_index().is_some() {
            return op();
        }
        match default_pool() {
            Some(pool) => pool.install(op),
            None => op(),
        }
    }

    pub fn current_threads() -> usize {
        if rayon::current_thread_index().is_some() {
            return rayon::current_num_threads();
        }
        default_pool().map_or(1, |p| p.current_num_threads())
    }

    pub fn map<T: Send>(count: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
        if count <= 1 {
            return (0..count).map(f).collect();
        }
        inside(|| (0..count).into_par_iter().with_max_len(1).map(&f).collect())
    }

    pub fn for_each_mut<T: Send>(items: &mut [T], f: impl Fn(usize, &mut T) + Sync + Send) {
        if items.len() <= 1 {
            items.iter_mut().enumerate().for_each(|(i, x)| f(i, x));
            return;
        }
        inside(|| items.par_iter_mut().with_max_len(1).enumerate().for_each(|(i, x)| f(i, x)));
    }
}

#[cfg(target_arch = "wasm32")]
mod pool {
    pub fn with_threads<R: Send>(_threads: usize, f: impl FnOnce() -> R + Send) -> R {
        f()
    }

    pub fn current_threads() -> usize {
        1
    }

    pub fn map<T: Send>(count: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
        (0..count).map(f).collect()
    }

    pub fn for_each_mut<T: Send>(items: &mut [T], f: impl Fn(usize, &mut T) + Sync + Send) {
        items.iter_mut().enumerate().for_each(|(i, x)| f(i, x));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A sum whose value depends on the order of its terms in floating point: the helpers must
    /// give the sequential answer on any number of threads.
    fn work() -> (Vec<f64>, f64) {
        let terms: Vec<f64> = (0..10_000).map(|i| ((i as f64) * 0.618_033_988_7).sin() * 10f64.powi(i % 17 - 8)).collect();
        let partial = map(blocks(terms.len(), 97).count(), |k| {
            let r = blocks(terms.len(), 97).nth(k).unwrap_or(0..0);
            terms[r].iter().sum::<f64>()
        });
        let total = partial.iter().sum::<f64>();
        let mut scaled = terms.clone();
        chunks_mut(&mut scaled, 61, |k, c| c.iter_mut().for_each(|x| *x *= 1.0 + k as f64 * 1e-3));
        (scaled, total)
    }

    #[test]
    fn results_do_not_depend_on_the_thread_count() {
        let one = with_threads(1, work);
        for threads in [2, 3, 7, 16] {
            let other = with_threads(threads, work);
            assert!(one.0.iter().zip(&other.0).all(|(a, b)| a.to_bits() == b.to_bits()));
            assert_eq!(one.1.to_bits(), other.1.to_bits(), "{threads} threads");
        }
        assert_eq!(work().1.to_bits(), one.1.to_bits(), "the default pool");
    }

    #[test]
    fn thread_counts_are_capped_and_parsed() {
        assert_eq!(with_threads(0, current_threads), 1);
        assert_eq!(with_threads(3, current_threads), 3);
        assert_eq!(with_threads(usize::MAX, current_threads), MAX_THREADS);
        assert_eq!(parse_threads(" 8 "), Some(8));
        for bad in ["0", "-1", "257", "1e3", "", "lots", "18446744073709551616"] {
            assert_eq!(parse_threads(bad), None, "{bad}");
        }
        assert!((1..=MAX_THREADS).contains(&default_threads()));
        // Nested use runs on the pool it is in.
        assert_eq!(with_threads(2, || map(3, |_| current_threads())), vec![2, 2, 2]);
        assert_eq!(with_threads(3, || map_grained(10, 3, |i| i * i)), (0..10).map(|i| i * i).collect::<Vec<_>>());
        assert_eq!(blocks(10, 4).collect::<Vec<_>>(), vec![0..4, 4..8, 8..10]);
        assert_eq!(blocks(0, 4).count(), 0);
        assert_eq!(blocks(3, 0).count(), 3);
    }

    /// An unrelated command must not hold the only worker needed by a cancelled command.
    /// Release the blocker before asserting, so this regression cannot strand a test thread.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn concurrent_commands_do_not_wait_for_another_commands_pool() {
        use std::sync::mpsc;
        use std::time::Duration;
        let (started_tx, started_rx) = mpsc::channel();
        let (release_tx, release_rx) = mpsc::channel();
        let (done_tx, done_rx) = mpsc::channel();
        std::thread::scope(|scope| {
            scope.spawn(move || {
                with_threads(1, move || {
                    started_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                })
            });
            started_rx.recv().unwrap();
            scope.spawn(move || with_threads(1, || done_tx.send(()).unwrap()));
            let completed = done_rx.recv_timeout(Duration::from_secs(1));
            release_tx.send(()).unwrap();
            assert!(completed.is_ok(), "an unrelated command monopolized the cached pool");
        });
    }

    /// While join waits for its right branch, its worker can steal another command.
    /// Finishing our right branch must let our command return even while that other command
    /// waits for its own event. This is the return path a cancelled command also uses.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn completed_work_is_not_held_by_an_unrelated_command() {
        use std::sync::mpsc;
        use std::time::Duration;
        let (right_tx, right_rx) = mpsc::channel();
        let (left_tx, left_rx) = mpsc::channel();
        let (finish_tx, finish_rx) = mpsc::channel();
        let (other_tx, other_rx) = mpsc::channel();
        let (release_tx, release_rx) = mpsc::channel();
        let (done_tx, done_rx) = mpsc::channel();
        std::thread::scope(|scope| {
            scope.spawn(move || {
                with_threads(2, move || {
                    rayon::join(
                        move || {
                            right_rx.recv().unwrap();
                            left_tx.send(()).unwrap();
                        },
                        move || {
                            right_tx.send(()).unwrap();
                            finish_rx.recv().unwrap();
                        },
                    );
                });
                done_tx.send(()).unwrap();
            });
            left_rx.recv().unwrap();
            scope.spawn(move || {
                with_threads(2, move || {
                    other_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                })
            });
            other_rx.recv().unwrap();
            finish_tx.send(()).unwrap();
            let completed = done_rx.recv_timeout(Duration::from_secs(1));
            release_tx.send(()).unwrap();
            assert!(completed.is_ok(), "completed work waited for an unrelated command stolen by its join");
        });
    }

    #[test]
    fn many_thread_counts_do_not_pile_up_pools() {
        for threads in 1..=12 {
            assert_eq!(with_threads(threads, current_threads), threads);
        }
    }
}
