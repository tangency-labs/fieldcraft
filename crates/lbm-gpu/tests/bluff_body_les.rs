//! Bluff-body LES reruns (#1106, #785 step 7): a circular cylinder at Re 3900 and a square cylinder at Re 22 000 on the
//! GPU, measured against Kravchenko & Moin (Phys. Fluids 12:403, 2000) and Lyn et al. (J. Fluid Mech. 304:285, 1995). The
//! test measures and prints; it asserts nothing about the accuracy, and a large error is reported, never tuned away.
//!
//! One case per run, from the environment variable `FIELDCRAFT_BLUFF_CASE` (`id:key=value,key=value`):
//! - `shape`: `cylinder` (default) or `square`;
//! - `re`: the Reynolds number on the body's diameter D (default 3900; 22000 for the square);
//! - `n`: cells per D (default 16); `span`, `height`, `length`, `up`: the spanwise width, the box height, the box length
//!   and the body's distance from the inlet, in D (defaults 3.125, 12 for the cylinder and 14 for the square, 24, 6);
//! - `sigma`: the hybrid recursive regularised collision's σ, a number or `dynamic` (default 0.99);
//! - `cs`: Vreman's Smagorinsky-equivalent constant (default 0.18, as Jacob, Malaspinas & Sagaut 2018);
//! - `start`, `end`: the averaging window in D/U (defaults 100 and 240);
//! - `ramp`: the inlet ramp in fundamental periods of the box, T = m T₀ with T₀ = 4 L / c_s (default 2.5, the zero of the
//!   fundamental mode's residual in `docs/manual/fluid/boundary-conditions.md`).
//!
//! Every run starts at rest on the body (a divergence-free field that stops at the body's surface), with interpolated
//! walls, the velocity faces raised from that field by the ramp. Printed per run: C_d with its 8-block 95 % interval,
//! C_l rms, the lift spectrum's peaks (St in cycles per D/U), the density range over the fluid cells, the recirculation
//! length behind the body, and the wall time. Skipped without a GPU adapter unless `FIELDCRAFT_REQUIRE_GPU=1`.
#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::Policy;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Outflow, Procedure, Profile, Turbulence, WallTreatment};

/// The lattice sound speed of D3Q19 in cells per step: c_s² = 1/3.
const SOUND_SPEED: f64 = 0.577_350_269_189_625_8;
/// Lattice speed of the bulk flow: the Mach number's order of 0.05 √3 ≈ 0.087 on the stream.
const LATTICE_SPEED: f64 = 0.05;
/// Steps between two force readings: the two-step mean is the odd and even steps' average (Ladd 1994).
const FORCE_STRIDE: u64 = 30;
/// Steps between field samples (density range, mean field).
const FIELD_STRIDE: u64 = 90;
/// Blocks of the averaging window for the statistical interval (one block is at least 10 D/U at these windows).
const BLOCKS: usize = 8;
/// Student's t at 95 % with BLOCKS − 1 degrees of freedom.
const T_95_DF7: f64 = 2.365;

#[derive(Clone, Debug, PartialEq)]
struct Case {
    id: String,
    square: bool,
    re: f64,
    n: usize,
    span: f64,
    height: f64,
    length: f64,
    up: f64,
    sigma: HrrSigma,
    cs: f64,
    start: f64,
    end: f64,
    ramp: f64,
}

impl Case {
    fn parse(spec: &str) -> Result<Case, String> {
        let (id, rest) = spec.split_once(':').unwrap_or((spec, ""));
        let mut c = Case {
            id: id.trim().to_string(),
            square: false,
            re: 3900.0,
            n: 16,
            span: 3.125,
            height: 12.0,
            length: 24.0,
            up: 6.0,
            sigma: HrrSigma::Value(0.99),
            cs: 0.18,
            start: 100.0,
            end: 240.0,
            ramp: 2.5,
        };
        for kv in rest.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            let (k, v) = kv.split_once('=').ok_or_else(|| format!("expected key=value, got {kv}"))?;
            let num = || v.trim().parse::<f64>().map_err(|_| format!("{k} needs a number, got {v}"));
            match k.trim() {
                "shape" => match v.trim() {
                    "cylinder" => c.square = false,
                    "square" => {
                        c.square = true;
                        c.re = 22_000.0;
                        c.height = 14.0;
                    }
                    other => return Err(format!("shape {other}")),
                },
                "re" => c.re = num()?,
                "n" => c.n = num()? as usize,
                "span" => c.span = num()?,
                "height" => c.height = num()?,
                "length" => c.length = num()?,
                "up" => c.up = num()?,
                "sigma" => {
                    c.sigma = match v.trim() {
                        "dynamic" => HrrSigma::Dynamic,
                        _ => HrrSigma::Value(num()?),
                    }
                }
                "cs" => c.cs = num()?,
                "start" => c.start = num()?,
                "end" => c.end = num()?,
                "ramp" => c.ramp = num()?,
                other => return Err(format!("unknown key {other}")),
            }
        }
        Ok(c)
    }
}

/// Every number of a run's geometry and timing that the run derives from its case.
struct Derived {
    dx: f64,
    span: f64,
    nz: usize,
    dt: f64,
    per_time: u64,
}

impl Derived {
    fn of(c: &Case) -> Derived {
        let n = c.n as f64;
        let dx = 1.0 / n;
        let nz = (c.span * n).round().max(1.0) as usize;
        let dt = LATTICE_SPEED * dx;
        Derived { dx, span: nz as f64 * dx, nz, dt, per_time: (n / LATTICE_SPEED).round() as u64 }
    }
}

/// Hann-windowed amplitude spectrum of `x` sampled every `dt`: the three largest peaks (frequency, power).
fn peaks(x: &[f64], dt: f64) -> Vec<(f64, f64)> {
    let n = x.len();
    let mean = x.iter().sum::<f64>() / n as f64;
    let w: Vec<f64> = (0..n).map(|i| (x[i] - mean) * 0.5 * (1.0 - (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos())).collect();
    let mut p = Vec::with_capacity(n / 2);
    for k in 1..n / 2 {
        let (mut re, mut im) = (0.0, 0.0);
        for (i, v) in w.iter().enumerate() {
            let a = 2.0 * std::f64::consts::PI * (k * i % n) as f64 / n as f64;
            re += v * a.cos();
            im -= v * a.sin();
        }
        p.push(re * re + im * im);
    }
    let mut local: Vec<(f64, f64)> =
        (1..p.len().saturating_sub(1)).filter(|&i| p[i] > p[i - 1] && p[i] >= p[i + 1]).map(|i| ((i + 1) as f64 / (n as f64 * dt), p[i])).collect();
    local.sort_by(|a, b| b.1.total_cmp(&a.1));
    local.truncate(3);
    local
}

/// The body's surface distance and the divergence-free start at rest on it: psi = U (y − y_c)(1 − e^(−d/0.5D)).
fn start_velocity(c: &Case, x: f64, y: f64) -> [f64; 2] {
    let (xc, yc) = (c.up, c.height / 2.0);
    let dist = |x: f64, y: f64| {
        if c.square {
            let (ax, ay) = (((x - xc).abs() - 0.5).max(0.0), ((y - yc).abs() - 0.5).max(0.0));
            (ax * ax + ay * ay).sqrt()
        } else {
            (((x - xc).powi(2) + (y - yc).powi(2)).sqrt() - 0.5).max(0.0)
        }
    };
    let psi = |x: f64, y: f64| (y - yc) * (1.0 - (-dist(x, y) / 0.5).exp());
    let h = 0.25 / c.n as f64;
    let u = (psi(x, y + h) - psi(x, y - h)) / (2.0 * h);
    let v = -(psi(x + h, y) - psi(x - h, y)) / (2.0 * h);
    [u, v]
}

/// A run the scheme cannot finish is a result: its line says where it stopped, and no mean is printed.
fn report_unstable(c: &Case, step: u64, per_time: u64, error: &str) {
    println!("RESULT {}: UNSTABLE at {:.1} D/U (step {step}): {error}", c.id, step as f64 / per_time as f64);
}

fn run(c: &Case) -> Result<(), String> {
    let clock = std::time::Instant::now();
    let n = c.n as f64;
    let d = Derived::of(c);
    let (xc, yc) = (c.up, c.height / 2.0);
    let mut dom = domain([0.0; 3], [c.length, c.height, d.span], d.dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective, sigma: None };
    dom.faces.y_min = FaceCondition::Slip {};
    dom.faces.y_max = FaceCondition::Slip {};
    dom.faces.z_min = FaceCondition::Periodic {};
    dom.faces.z_max = FaceCondition::Periodic {};
    dom.obstacles = vec!["Body".into()];
    let body = if c.square {
        obstacle("Body", cuboid([xc - 0.5, yc - 0.5, -d.dx], [xc + 0.5, yc + 0.5, d.span + d.dx]), &CUBOID_FACES)
    } else {
        obstacle("Body", cylinder([xc, yc], 0.5, -d.dx, d.span + d.dx, 1024), &["side", "bottom", "top"])
    };
    // The inlet rises over m fundamental periods of the box, T₀ = 4 L / c_s with L the box's length in cells.
    let nx = (c.length * n).round();
    let ramp_time = c.ramp * 4.0 * nx / SOUND_SPEED * d.dt;
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(c.end),
        steps: None,
        time_step: Some(d.dt),
        frames: 1,
        magic: None,
        sigma: Some(c.sigma),
        collision: Some(Collision::Hrr),
        turbulence: Some(Turbulence::Vreman { constant: c.cs, prandtl: 0.9, wall_model: None }),
        steady: None,
        walls: WallTreatment::Interpolated,
        inlet_ramp: Some(ramp_time),
        trace_bulk_one: false,
    };
    let m = model(1.0, 1.0 / c.re, dom, procedure);
    let limits = Limits { max_cells: 120_000_000, max_cell_updates: 20_000_000_000_000, max_steps: 4_000_000, ..Limits::default() };
    let obstacles = [body];
    let setup = Setup::new(&m, None, &obstacles, &limits).map_err(|e| format!("{}: refused at setup: {e}", c.id))?;
    for w in &setup.warnings {
        println!("  warning: {w}");
    }
    let mut flow = Flow::new(&setup, &obstacles).map_err(|e| format!("{}: {e}", c.id))?;
    let grid = *flow.grid();
    let [nxc, nyc, nzc] = grid.cells;
    flow.set_equilibrium(|cell| {
        let [x, y, _] = grid.centre(cell);
        let [u, v] = start_velocity(c, x, y);
        (1.0, [LATTICE_SPEED * u, LATTICE_SPEED * v, 0.0])
    });
    let (tau, _) = flow.omegas();
    println!(
        "CASE {}: shape {}, Re {}, {} cells per D, grid {nxc} x {nyc} x {nzc} = {} cells ({} fluid), tau {:.6}, steps {}, \
         ramp {:.1} periods = {:.1} D/U",
        c.id,
        if c.square { "square" } else { "cylinder" },
        c.re,
        c.n,
        grid.len(),
        flow.fluid_cells(),
        1.0 / tau,
        setup.steps,
        c.ramp,
        ramp_time
    );
    let fluid: Vec<bool> = (0..grid.len()).map(|i| flow.is_fluid(i)).collect();
    let policy = Policy::from_env().map_err(|e| format!("{e}"))?;
    let mut gpu =
        GpuFlow::with_options(&flow, policy, Options { streaming: Streaming::InPlace, profile: false }).map_err(|e| format!("{}: {e}", c.id))?;
    println!("  GPU: {}", gpu.adapter_name());
    let groups: Vec<usize> = gpu.groups().iter().enumerate().filter(|(_, g)| g.body.as_deref() == Some("Body")).map(|(i, _)| i).collect();
    let per_time = d.per_time;
    let start = (c.start * per_time as f64).round() as u64;
    let coefficient = 2.0 / (LATTICE_SPEED * LATTICE_SPEED * n * d.nz as f64);
    let (mut cd, mut cl) = (Vec::new(), Vec::new());
    let mut rho_mean = vec![0.0f64; nxc * nyc];
    let mut u_mean = vec![0.0f64; nxc * nyc];
    let mut samples = 0usize;
    let (mut rho_lo, mut rho_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut k = 0u64;
    while k < setup.steps {
        let next = (k + FORCE_STRIDE).min(setup.steps);
        if let Err(e) = gpu.advance(next - k - 1) {
            report_unstable(c, k, per_time, &e.to_string());
            return Ok(());
        }
        let f_odd = gpu.forces().map_err(|e| format!("{e}"))?;
        if let Err(e) = gpu.advance(1) {
            report_unstable(c, k, per_time, &e.to_string());
            return Ok(());
        }
        k = next;
        let f_even = gpu.forces().map_err(|e| format!("{e}"))?;
        let sum = |f: &[[f64; 3]]| groups.iter().fold((0.0, 0.0), |a, &g| (a.0 + f[g][0], a.1 + f[g][1]));
        let ((ox, oy), (ex, ey)) = (sum(&f_odd), sum(&f_even));
        if k >= start {
            cd.push(0.5 * (ox + ex) * coefficient);
            cl.push(0.5 * (oy + ey) * coefficient);
        }
        if k.is_multiple_of(FIELD_STRIDE) {
            let fields = gpu.fields().map_err(|e| format!("{e}"))?;
            for (cell, &r) in fields.density.iter().enumerate() {
                if fluid[cell] {
                    rho_lo = rho_lo.min(r);
                    rho_hi = rho_hi.max(r);
                }
            }
            if k >= start {
                for cell in 0..grid.len() {
                    let [i, j, _] = grid.coords(cell);
                    rho_mean[i + nxc * j] += fields.density[cell];
                    u_mean[i + nxc * j] += fields.velocity[cell][0];
                }
                samples += 1;
            }
        }
    }
    let wall = clock.elapsed().as_secs_f64();
    if cd.len() < BLOCKS * 2 {
        return Err(format!("{}: window too short ({} force samples)", c.id, cd.len()));
    }
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let rms = |x: &[f64]| (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt();
    let blk = cd.len() / BLOCKS;
    let blocks: Vec<f64> = (0..BLOCKS).map(|b| mean(&cd[b * blk..(b + 1) * blk])).collect();
    let bm = mean(&blocks);
    let bs = (blocks.iter().map(|v| (v - bm).powi(2)).sum::<f64>() / (BLOCKS as f64 - 1.0)).sqrt();
    let half_width = T_95_DF7 * bs / (BLOCKS as f64).sqrt();
    let sample_dt = FORCE_STRIDE as f64 / per_time as f64;
    let lift_peaks = peaks(&cl, sample_dt);
    let st = lift_peaks.first().map_or(f64::NAN, |p| p.0);
    // Recirculation: the mean streamwise velocity on the wake centreline; L_r from the rear face.
    let at = |a: &[f64], x: f64, y: f64| -> f64 {
        let s = samples as f64;
        let (fx, fy) = ((x / d.dx - 0.5).clamp(0.0, nxc as f64 - 1.001), (y / d.dx - 0.5).clamp(0.0, nyc as f64 - 1.001));
        let (i, j) = (fx.floor() as usize, fy.floor() as usize);
        let (tx, ty) = (fx - i as f64, fy - j as f64);
        let g = |i: usize, j: usize| a[i + nxc * j] / s;
        (1.0 - tx) * (1.0 - ty) * g(i, j) + tx * (1.0 - ty) * g(i + 1, j) + (1.0 - tx) * ty * g(i, j + 1) + tx * ty * g(i + 1, j + 1)
    };
    let rear = xc + 0.5;
    let mut lr = f64::NAN;
    let mut prev: Option<(f64, f64)> = None;
    let mut x = rear + 1.5 * d.dx;
    while x < (xc + 8.0).min(c.length - d.dx) {
        let u = at(&u_mean, x, yc) / LATTICE_SPEED;
        if let Some((px, pu)) = prev
            && lr.is_nan()
            && pu < 0.0
            && u >= 0.0
        {
            lr = px + (x - px) * (-pu) / (u - pu);
        }
        prev = Some((x, u));
        x += 0.25 * d.dx * n;
    }
    let lr_d = lr - rear;
    let cd_mean = mean(&cd);
    println!(
        "RESULT {}: C_d {:.4} +- {:.4} (blocks {:?}, lag-1 of blocks {:.2}); C_l rms {:.3}; St {:.4} (lift peaks {:?}); \
         L_r/D {:.3} (from the rear, {} field samples); rho {:.4}..{:.4} over the fluid cells; wall {:.0} s ({:.2} h)",
        c.id,
        cd_mean,
        half_width,
        blocks.iter().map(|b| (b * 1000.0).round() / 1000.0).collect::<Vec<_>>(),
        blocks.windows(2).map(|w| (w[0] - bm) * (w[1] - bm)).sum::<f64>() / blocks.iter().map(|v| (v - bm).powi(2)).sum::<f64>(),
        rms(&cl.iter().map(|v| v - mean(&cl)).collect::<Vec<_>>()),
        st,
        lift_peaks.iter().map(|p| ((p.0 * 1e4).round() / 1e4, p.1.round())).collect::<Vec<_>>(),
        lr_d,
        samples,
        rho_lo,
        rho_hi,
        wall,
        wall / 3600.0
    );
    Ok(())
}

#[test]
#[ignore = "heavy: GPU bluff-body LES rerun (#1106); set FIELDCRAFT_BLUFF_CASE, run on the GB10 via gpurun"]
fn bluff_body_les_rerun() {
    let Ok(spec) = std::env::var("FIELDCRAFT_BLUFF_CASE") else {
        println!("set FIELDCRAFT_BLUFF_CASE (id:key=value,...)");
        return;
    };
    let case = Case::parse(&spec).expect("a valid FIELDCRAFT_BLUFF_CASE");
    run(&case).expect("the bluff-body run completes");
}

#[test]
fn a_case_spec_sets_the_reference_defaults_and_the_square_re() {
    let cyl = Case::parse("cyl").expect("a bare id is a cylinder");
    assert!(!cyl.square);
    assert_eq!((cyl.re, cyl.n, cyl.height, cyl.up), (3900.0, 16, 12.0, 6.0));
    assert_eq!(cyl.sigma, HrrSigma::Value(0.99));
    let sq = Case::parse("sq:shape=square,sigma=dynamic").expect("a square case");
    assert!(sq.square);
    assert_eq!((sq.re, sq.height, sq.sigma), (22_000.0, 14.0, HrrSigma::Dynamic));
    assert!(Case::parse("x:shape=triangle").is_err());
    assert!(Case::parse("x:n").is_err());
}

#[test]
fn the_inlet_ramp_is_m_fundamental_periods_of_the_box() {
    // T₀ = 4 L / c_s steps: for a box of 384 cells, m · 4 · 384 · √3 steps at m = 2.5 (the default).
    let c = Case::parse("cyl").expect("a case");
    let nx = (c.length * c.n as f64).round();
    let steps = c.ramp * 4.0 * nx / SOUND_SPEED;
    assert!((steps - 2.5 * 4.0 * 384.0 * 3f64.sqrt()).abs() < 1e-9);
}
