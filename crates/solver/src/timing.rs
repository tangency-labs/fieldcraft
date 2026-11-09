//! Opt-in diagnostic wall times. No solver operands or scheduling depend on these clocks.
//! Records are inclusive (in particular, parallel solve spans must not be summed as wall time).
use std::io::Write;
use std::sync::OnceLock;
use std::time::Instant;

fn enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var("FIELDCRAFT_SOLVER_TIMING").is_ok_and(|v| v == "1"))
}

/// Successful triangular traversals: RHS columns, executed lanes, equations, coefficient
/// visits, scalar multiply/subtract updates and permuted scalar values. These are logical
/// work counts, not measured cache misses, memory transactions or CPU time.
pub(crate) fn triangular_work(rhs: usize, lanes: usize, n: usize, entries: usize) {
    if enabled() {
        let visits = entries.saturating_mul(2);
        let updates = entries.saturating_sub(n).saturating_mul(2).saturating_mul(lanes);
        let copies = n.saturating_mul(rhs).saturating_mul(2);
        let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_TRIANGULAR_WORK\t{rhs}\t{lanes}\t{n}\t{visits}\t{updates}\t{copies}");
    }
}

pub(crate) struct Span {
    name: &'static str,
    start: Option<Instant>,
}

impl Span {
    pub(crate) fn new(name: &'static str) -> Self {
        Self { name, start: enabled().then(Instant::now) }
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            let seconds = start.elapsed().as_secs_f64();
            // Diagnostics must never turn an otherwise successful solve into an I/O error.
            let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_TIMING\t{}\t{seconds:.9}", self.name);
        }
    }
}

/// Per-supernode diagnostics. Times of concurrent nodes are not additive wall time.
/// Product throughput includes packing, allocation and worker scheduling; it is not a
/// hardware instruction counter. FLOPs count multiply and add separately, including padding.
pub(crate) struct Supernode {
    start: Option<Instant>,
    j: usize,
    m: usize,
    w: usize,
    pub(crate) assembly_seconds: f64,
    product_seconds: f64,
    pub(crate) assembly_product_seconds: f64,
    product_flops: f64,
    below_seconds: f64,
}

impl Supernode {
    pub(crate) fn new(j: usize, m: usize, w: usize) -> Self {
        static ENABLED: OnceLock<bool> = OnceLock::new();
        let enabled = *ENABLED.get_or_init(|| std::env::var("FIELDCRAFT_SUPERNODE_TIMING").is_ok_and(|v| v == "1"));
        Self {
            start: enabled.then(Instant::now),
            j,
            m,
            w,
            assembly_seconds: 0.0,
            product_seconds: 0.0,
            assembly_product_seconds: 0.0,
            product_flops: 0.0,
            below_seconds: 0.0,
        }
    }

    pub(crate) fn clock(&self) -> Option<Instant> {
        self.start.map(|_| Instant::now())
    }

    pub(crate) fn elapsed(start: Option<Instant>) -> f64 {
        start.map_or(0.0, |s| s.elapsed().as_secs_f64())
    }

    pub(crate) fn product(&mut self, start: Option<Instant>, mu: usize, wu: usize, kw: usize, tile: usize, nr: usize) {
        self.product_seconds += Self::elapsed(start);
        if self.start.is_some() {
            for i0 in (0..mu).step_by(tile) {
                let cols = wu.min(i0 + tile.min(mu - i0));
                self.product_flops += 2.0 * (cols.div_ceil(nr) * nr) as f64 * tile as f64 * kw as f64;
            }
        }
    }

    pub(crate) fn finish_assembly(&mut self) {
        self.assembly_product_seconds = self.product_seconds;
    }

    pub(crate) fn below(&mut self, start: Option<Instant>) {
        self.below_seconds += Self::elapsed(start);
    }
}

impl Drop for Supernode {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            let seconds = start.elapsed().as_secs_f64();
            let _ = writeln!(
                std::io::stderr().lock(),
                "FIELDCRAFT_SUPERNODE\t{}\t{}\t{}\t{seconds:.9}\t{:.9}\t{:.9}\t{:.0}\t{:.9}\t{:.9}",
                self.j,
                self.m,
                self.w,
                self.assembly_seconds,
                self.product_seconds,
                self.product_flops,
                self.below_seconds,
                self.assembly_product_seconds
            );
        }
    }
}
