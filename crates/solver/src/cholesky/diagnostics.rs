//! Isolated, opt-in test experiment; no production solver selection or ISA changes.
//! Logical byte counts model kernel operands, not measured DRAM/cache traffic.
use super::*;
use std::sync::Mutex;
use std::time::Instant;
#[path = "diagnostic_baseline.rs"]
mod baseline;
#[path = "diagnostic_candidate.rs"]
mod candidate;

#[derive(Clone, Copy, Default)]
struct Metric {
    seconds: f64,
    calls: u64,
    bytes: u64,
}
thread_local! { static LOCAL: std::cell::RefCell<[Metric;4]> = const { std::cell::RefCell::new([Metric {seconds:0.0,calls:0,bytes:0};4]) }; }
struct Stats {
    enabled: bool,
    metrics: Mutex<[Metric; 4]>,
    dispatch: Mutex<Vec<(f64, f64, f64, usize)>>,
    pack_wall: Mutex<Vec<f64>>,
}
impl Stats {
    fn new(enabled: bool) -> Self {
        Self { enabled, metrics: Mutex::new([Metric::default(); 4]), dispatch: Mutex::new(vec![]), pack_wall: Mutex::new(vec![]) }
    }
    fn flush(&self) {
        if self.enabled {
            let values = LOCAL.with(|cell| std::mem::take(&mut *cell.borrow_mut()));
            let mut total = self.metrics.lock().unwrap();
            for (a, b) in total.iter_mut().zip(values) {
                a.seconds += b.seconds;
                a.calls += b.calls;
                a.bytes += b.bytes;
            }
        }
    }
    fn clock(&self) -> Option<Instant> {
        self.enabled.then(Instant::now)
    }
    fn record(&self, kind: usize, start: Instant, calls: u64, bytes: u64) {
        self.record_clock(kind, self.enabled.then_some(start), calls, bytes);
    }
    fn record_clock(&self, kind: usize, start: Option<Instant>, calls: u64, bytes: u64) {
        if let Some(start) = start {
            let elapsed = start.elapsed().as_secs_f64();
            LOCAL.with(|cell| {
                let mut m = cell.borrow_mut();
                m[kind].seconds += elapsed;
                m[kind].calls += calls;
                m[kind].bytes += bytes;
            });
        }
    }
}
fn dispatch_map<T: Send>(n: usize, single: bool, stats: &Stats, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    let setup = stats.clock();
    let start = Instant::now();
    let spans = Mutex::new(Vec::new());
    let op = |i| {
        let begin = Instant::now();
        let value = f(i);
        stats.flush();
        if stats.enabled {
            spans.lock().unwrap().push((begin.duration_since(start).as_secs_f64(), begin.elapsed().as_secs_f64()));
        }
        value
    };
    stats.record_clock(2, setup, 1, (n * std::mem::size_of::<usize>()) as u64);
    let result = if single { (0..n).map(op).collect() } else { parallel::map(n, op) };
    if stats.enabled {
        let s = spans.into_inner().unwrap();
        stats.dispatch.lock().unwrap().push((start.elapsed().as_secs_f64(), s.iter().map(|x| x.0).sum(), s.iter().map(|x| x.1).sum(), n));
    }
    result
}
fn dispatch_each<T: Send>(items: &mut [T], single: bool, stats: &Stats, f: impl Fn(usize, &mut T) + Sync + Send) {
    let setup = stats.clock();
    let start = Instant::now();
    let spans = Mutex::new(Vec::new());
    let n = items.len();
    let op = |i, x: &mut T| {
        let begin = Instant::now();
        f(i, x);
        stats.flush();
        if stats.enabled {
            spans.lock().unwrap().push((begin.duration_since(start).as_secs_f64(), begin.elapsed().as_secs_f64()));
        }
    };
    stats.record_clock(2, setup, 1, (n * std::mem::size_of::<usize>()) as u64);
    if single {
        for (i, x) in items.iter_mut().enumerate() {
            op(i, x)
        }
    } else {
        parallel::for_each_mut(items, op)
    }
    if stats.enabled {
        let s = spans.into_inner().unwrap();
        stats.dispatch.lock().unwrap().push((start.elapsed().as_secs_f64(), s.iter().map(|x| x.0).sum(), s.iter().map(|x| x.1).sum(), n));
    }
}
fn dispatch_chunks<T: Send>(items: &mut [T], chunk: usize, single: bool, stats: &Stats, f: impl Fn(usize, &mut [T]) + Sync + Send) {
    let mut parts: Vec<_> = items.chunks_mut(chunk).collect();
    dispatch_each(&mut parts, single, stats, |i, p| f(i, p));
}

#[test]
#[ignore = "explicit Dell v3 diagnostic; output is kernel-command time, not solve.static"]
fn four_phase_update_experiment() {
    let tiny = std::env::var("FC_DIAG_SIZE").as_deref() == Ok("tiny");
    let (mu, wu, kw) = if tiny { (129, 67, 65) } else { (2048, 1024, 512) };
    let lm = kw + mu;
    let data: Vec<f64> = (0..trapezoid(lm, kw)).map(|i| ((i * 37 % 101) as f64 - 50.0) / 101.0).collect();
    let relative: Vec<_> = (0..mu).map(|i| 3 + i).collect();
    let m = mu + 5;
    let initial = vec![0.75; trapezoid(m, relative[wu - 1] + 1)];
    let panel = Panel { m: lm, data: &data };
    let mut reference = initial.clone();
    let prod = products(&panel, 0..kw, kw, mu, wu, None);
    subtract(&mut reference, m, &prod, mu, wu, |i| relative[i]);
    drop(prod);
    let variant = std::env::var("FC_DIAG_VARIANT").unwrap_or_else(|_| "baseline".into());
    let threads = std::env::var("FC_DIAG_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or(1);
    let enabled = std::env::var("FC_DIAG_PROFILE").as_deref() != Ok("0");
    let count = if tiny { 2 } else { 4 };
    let stats = Stats::new(enabled);
    let mut output = vec![initial; count];
    let start = Instant::now();
    parallel::with_threads(threads, || {
        parallel::for_each_mut(&mut output, |_, out| {
            if variant == "baseline" {
                let prod = baseline::baseline_products(&panel, 0..kw, kw, mu, wu, None, &stats, false);
                let clock = stats.clock();
                subtract(out, m, &prod, mu, wu, |i| relative[i]);
                stats.record_clock(3, clock, 1, (wu * (2 * mu - wu + 1) / 2 * 24) as u64);
            } else {
                candidate::update(
                    &panel,
                    0..kw,
                    kw,
                    mu,
                    wu,
                    candidate::Target { m, first: 0, relative: &relative, values: out },
                    &stats,
                    variant == "single",
                    variant == "kernel",
                )
                .unwrap();
            }
            stats.flush();
        })
    });
    let elapsed = start.elapsed().as_secs_f64();
    stats.flush();
    for out in &output {
        assert!(out.iter().zip(&reference).all(|(a, b)| a.to_bits() == b.to_bits()), "{variant} threads {threads}");
    }
    let hash = reference.iter().fold(0xcbf29ce484222325u64, |h, v| (h ^ v.to_bits()).wrapping_mul(0x100000001b3));
    eprintln!(
        "DIAG_SUMMARY variant={variant} threads={threads} profile={enabled} mu={mu} wu={wu} kw={kw} independent_updates={count} command_seconds={elapsed:.9} bitwise=true relative_error=0 hash={hash:016x}"
    );
    for (i, name) in ["pack", "dense", "schedule_setup", "writeback"].iter().enumerate() {
        let v = stats.metrics.lock().unwrap()[i];
        eprintln!("DIAG_PHASE name={name} span_sum_seconds={:.9} calls={} logical_bytes={}", v.seconds, v.calls, v.bytes);
    }
    let d = stats.dispatch.lock().unwrap();
    eprintln!(
        "DIAG_DISPATCH regions={} tasks={} envelope_sum_seconds={:.9} task_span_sum_seconds={:.9} ready_to_start_sum_seconds={:.9}",
        d.len(),
        d.iter().map(|v| v.3).sum::<usize>(),
        d.iter().map(|v| v.0).sum::<f64>(),
        d.iter().map(|v| v.2).sum::<f64>(),
        d.iter().map(|v| v.1).sum::<f64>()
    );
}
