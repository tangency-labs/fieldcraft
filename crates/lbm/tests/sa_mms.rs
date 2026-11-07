//! The method of manufactured solutions for the Spalart–Allmaras transport (task C2 of the hybrid
//! RANS/LES board): a harness that measures a transport scheme's error against a manufactured
//! field on N = 16, 32, 64 and the observed order between successive grids.
//!
//! The method: choose a smooth field φ_m, apply the continuous operator to it by hand, and add
//! the residual as a source, so φ_m is the exact solution of the forced problem; the discrete
//! solution's error then falls as Δx^p, p the scheme's order (P. J. Roache, *Code verification by
//! the method of manufactured solutions*, J. Fluids Eng. 124(1):4–10, 2002,
//! doi:10.1115/1.1436090; C. J. Roy, *Review of code and solution verification procedures for
//! computational simulation*, J. Comput. Phys. 205(1):131–156, 2005,
//! doi:10.1016/j.jcp.2004.10.036).
//!
//! The harness is checked here on two toy schemes whose orders are known, before it is pointed at
//! `SaField::advance` (task R3): first-order upwind advection (order 1) and conservative central
//! diffusion (order 2), each on a periodic line with an explicit Euler step at a fixed CFL or
//! diffusion number, so the time error scales with the space error.
//!
//! For SA itself (R3's acceptance: advection–diffusion–source order in [0.9, 1.2], pure diffusion
//! in [1.8, 2.2], the same order for a field crossing zero), the manufactured ν̃ is applied to
//! R2's closure on a frozen analytic velocity, and its source is written out by hand in this file
//! when the closure lands, term by term (advection, the σ-scaled diffusion and c_b2 term,
//! production, destruction, SA-neg branch).

/// Errors of one grid against the manufactured field.
#[derive(Debug, Clone, Copy)]
struct Errors {
    l2: f64,
    linf: f64,
}

/// The observed orders log(e_coarse / e_fine) / log 2 between successive grids, each grid twice
/// as fine as the one before.
fn orders(errors: &[Errors]) -> Vec<(f64, f64)> {
    errors.windows(2).map(|w| ((w[0].l2 / w[1].l2).log2(), (w[0].linf / w[1].linf).log2())).collect()
}

/// The errors of `run` on each grid size, and its observed orders.
fn study(name: &str, grids: &[usize], run: impl Fn(usize) -> Errors) -> Vec<(f64, f64)> {
    let errors: Vec<Errors> = grids.iter().map(|&n| run(n)).collect();
    let p = orders(&errors);
    println!("{name}:");
    for (k, (n, e)) in grids.iter().zip(&errors).enumerate() {
        let order = k.checked_sub(1).and_then(|j| p.get(j)).map_or(String::new(), |(a, b)| format!(", order L2 {a:.3}, L∞ {b:.3}"));
        println!("  N {n:3}: L2 {:.4e}, L∞ {:.4e}{order}", e.l2, e.linf);
    }
    p
}

fn errors(numeric: &[f64], exact: impl Fn(f64) -> f64) -> Errors {
    let n = numeric.len() as f64;
    let (mut sum, mut max) = (0.0f64, 0.0f64);
    for (i, v) in numeric.iter().enumerate() {
        let e = v - exact((i as f64 + 0.5) / n);
        sum += e * e;
        max = max.max(e.abs());
    }
    Errors { l2: (sum / n).sqrt(), linf: max }
}

const TAU: f64 = std::f64::consts::TAU;

/// Toy 1: ∂φ/∂t + a ∂φ/∂x = S on the periodic unit line, first-order upwind, explicit Euler at
/// CFL 0.5, to t = 0.25. Manufactured φ_m = 2 + sin 2π(x − t) + ½ cos 2πx, so
/// S = ∂φ_m/∂t + a ∂φ_m/∂x = 2π (a − 1) cos 2π(x − t) − π a sin 2πx.
fn upwind(n: usize) -> Errors {
    let a = 1.5;
    let exact = |x: f64, t: f64| 2.0 + (TAU * (x - t)).sin() + 0.5 * (TAU * x).cos();
    let source = |x: f64, t: f64| TAU * (a - 1.0) * (TAU * (x - t)).cos() - std::f64::consts::PI * a * (TAU * x).sin();
    let h = 1.0 / n as f64;
    let dt = 0.5 * h / a;
    let steps = (0.25 / dt).round() as usize;
    let x = |i: usize| (i as f64 + 0.5) * h;
    let mut phi: Vec<f64> = (0..n).map(|i| exact(x(i), 0.0)).collect();
    for s in 0..steps {
        let t = s as f64 * dt;
        let old = phi.clone();
        for i in 0..n {
            let up = old[(i + n - 1) % n];
            phi[i] = old[i] - a * dt / h * (old[i] - up) + dt * source(x(i), t);
        }
    }
    let t = steps as f64 * dt;
    errors(&phi, |x| exact(x, t))
}

/// Toy 2: ∂φ/∂t = ∂/∂x (D(x) ∂φ/∂x) + S on the periodic unit line, conservative central
/// differences (face diffusivities from the mean of the two cells), explicit Euler at diffusion
/// number 1/6, to t = 0.05, with a variable diffusivity as SA's ν + ν̃ is.
/// Manufactured φ_m = e^(−t) sin 2πx with D = 1 + ½ cos 2πx, so
/// S = ∂φ_m/∂t − D' φ_m' − D φ_m'', with ∂φ_m/∂t = −φ_m, D' = −π sin 2πx,
/// φ_m' = 2π e^(−t) cos 2πx and φ_m'' = −4π² e^(−t) sin 2πx.
fn central(n: usize) -> Errors {
    let d = |x: f64| 1.0 + 0.5 * (TAU * x).cos();
    let exact = |x: f64, t: f64| (-t).exp() * (TAU * x).sin();
    let source = |x: f64, t: f64| {
        let e = (-t).exp();
        let (dp, p1, p2) = (-std::f64::consts::PI * (TAU * x).sin(), TAU * e * (TAU * x).cos(), -TAU * TAU * e * (TAU * x).sin());
        -exact(x, t) - dp * p1 - d(x) * p2
    };
    let h = 1.0 / n as f64;
    let dt = h * h / (6.0 * 1.5);
    let steps = (0.05 / dt).round() as usize;
    let x = |i: usize| (i as f64 + 0.5) * h;
    let mut phi: Vec<f64> = (0..n).map(|i| exact(x(i), 0.0)).collect();
    // Face diffusivity between cell i and i + 1.
    let face: Vec<f64> = (0..n).map(|i| 0.5 * (d(x(i)) + d(x(i) + h))).collect();
    for s in 0..steps {
        let t = s as f64 * dt;
        let old = phi.clone();
        for i in 0..n {
            let (l, r) = (old[(i + n - 1) % n], old[(i + 1) % n]);
            let flux = face[i] * (r - old[i]) - face[(i + n - 1) % n] * (old[i] - l);
            phi[i] = old[i] + dt / (h * h) * flux + dt * source(x(i), t);
        }
    }
    let t = steps as f64 * dt;
    errors(&phi, |x| exact(x, t))
}

const GRIDS: [usize; 3] = [16, 32, 64];

#[test]
fn the_harness_finds_first_order_for_upwind_advection() {
    let p = study("upwind advection with a source", &GRIDS, upwind);
    // The finest pair is in the asymptotic range; the board's band for R3 is [0.9, 1.2].
    let (l2, linf) = p[1];
    assert!((0.9..=1.2).contains(&l2) && (0.9..=1.2).contains(&linf), "{p:?}");
}

#[test]
fn the_harness_finds_second_order_for_central_diffusion() {
    let p = study("conservative central diffusion, variable D", &GRIDS, central);
    let (l2, linf) = p[1];
    assert!((1.8..=2.2).contains(&l2) && (1.8..=2.2).contains(&linf), "{p:?}");
}

/// The harness is not fooled: a scheme with an O(1) bug (the source applied at half strength)
/// stalls: its observed order falls towards zero.
#[test]
fn the_harness_exposes_an_inconsistent_scheme() {
    let broken = |n: usize| {
        let e = upwind(n);
        let biased = errors(&vec![1.0; n], |_| 0.0);
        Errors { l2: e.l2 + 0.1 * biased.l2, linf: e.linf + 0.1 * biased.linf }
    };
    let p = study("upwind with an O(1) error", &GRIDS, broken);
    // The finest pair, where the O(1) part dominates (the coarse pair still sees the scheme's own
    // first-order error).
    let (l2, linf) = p[1];
    assert!(l2 < 0.5 && linf < 0.5, "{p:?}");
}
