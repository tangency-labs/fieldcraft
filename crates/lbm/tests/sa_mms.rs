//! Manufactured-solution order tests for the Spalart–Allmaras transport (Roache 2002, J. Fluids Eng. 124(1):4–10;
//! Roy 2005, J. Comput. Phys. 205(1):131–156). The harness is checked on two toy schemes of known order, then pointed
//! at `SaField::advance` through its wall and inflow faces with production and destruction off. Acceptance on the
//! finest pair of N = 16, 32, 64: advection–diffusion order in [0.9, 1.2], diffusion alone in [1.8, 2.2], in L2 and L∞.

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
    // The finest pair is in the asymptotic range; the accepted band is [0.9, 1.2].
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

// ---------------------------------------------------------------------------------------------
// SaField's transport through its wall and inflow faces.
// ---------------------------------------------------------------------------------------------

use fieldcraft_lbm::sa::{BoxFace, C_B2, SIGMA, SaField, SaInput};

/// A manufactured ν̃ on the unit square (X along x, Z along z), F = F₀ + A e^(−τ) G(X, Z), and the
/// first and second derivatives of G, written out by hand.
struct Manufactured {
    f0: f64,
    g: fn(f64, f64) -> [f64; 5],
}

/// The amplitude A, the lattice viscosity ν and the advection speed a (unit-square units).
const AMPLITUDE: f64 = 0.05;
const NU: f64 = 0.05;
const SPEED: f64 = 1.0;
const PI: f64 = std::f64::consts::PI;

/// Next to a wall: G = (1 + ½ sin 2πX) sin(πZ/2), periodic in X, zero at the wall Z = 0 and with
/// zero slope at the open face Z = 1. Returns [G, G_X, G_Z, G_XX, G_ZZ]:
/// G_X = π cos 2πX sin(πZ/2), G_Z = (π/2)(1 + ½ sin 2πX) cos(πZ/2),
/// G_XX = −2π² sin 2πX sin(πZ/2), G_ZZ = −(π²/4)(1 + ½ sin 2πX) sin(πZ/2).
fn wall_g(x: f64, z: f64) -> [f64; 5] {
    let (s, c) = ((TAU * x).sin(), (TAU * x).cos());
    let (sz, cz) = ((0.5 * PI * z).sin(), (0.5 * PI * z).cos());
    let w = 1.0 + 0.5 * s;
    [w * sz, PI * c * sz, 0.5 * PI * w * cz, -2.0 * PI * PI * s * sz, -0.25 * PI * PI * w * sz]
}

/// Behind an inflow face: G = sin(πX/2)(1 + ½ sin 2πZ), zero at the inflow X = 0 (so ν̃ is the
/// inflow value F₀ there), zero slope at the open outlet X = 1, periodic in Z. Returns
/// [G, G_X, G_Z, G_XX, G_ZZ]: G_X = (π/2) cos(πX/2)(1 + ½ sin 2πZ), G_Z = π sin(πX/2) cos 2πZ,
/// G_XX = −(π²/4) sin(πX/2)(1 + ½ sin 2πZ), G_ZZ = −2π² sin(πX/2) sin 2πZ.
fn inflow_g(x: f64, z: f64) -> [f64; 5] {
    let (sx, cx) = ((0.5 * PI * x).sin(), (0.5 * PI * x).cos());
    let (s, c) = ((TAU * z).sin(), (TAU * z).cos());
    let w = 1.0 + 0.5 * s;
    [sx * w, 0.5 * PI * cx * w, PI * sx * c, -0.25 * PI * PI * sx * w, -2.0 * PI * PI * sx * s]
}

impl Manufactured {
    fn value(&self, x: f64, z: f64, t: f64) -> f64 {
        self.f0 + AMPLITUDE * (-t).exp() * (self.g)(x, z)[0]
    }

    /// The source that makes F exact. SA's transport with the model off, on the unit square in
    /// the diffusive time τ:
    /// ∂F/∂τ + a ∂F/∂X = (1/σ)[∇·((ν + F)∇F) + c_b2 |∇F|²] + S,
    /// and ∇·((ν + F)∇F) = (ν + F)∇²F + |∇F|², so
    /// S = ∂F/∂τ + a F_X − (1/σ)[(ν + F)(F_XX + F_ZZ) + (1 + c_b2)(F_X² + F_Z²)],
    /// with ∂F/∂τ = −A e^(−τ) G and every other derivative A e^(−τ) times G's.
    fn source(&self, x: f64, z: f64, t: f64, speed: f64) -> f64 {
        let e = AMPLITUDE * (-t).exp();
        let [g, gx, gz, gxx, gzz] = (self.g)(x, z);
        let f = self.f0 + e * g;
        let (fx, fz) = (e * gx, e * gz);
        -e * g + speed * fx - ((NU + f) * e * (gxx + gzz) + (1.0 + C_B2) * (fx * fx + fz * fz)) / SIGMA
    }
}

/// Run `field` on an N × 1 × N box with the given faces to τ = 1/8 and measure it against `m`.
///
/// Lattice units, Δx = 1: the cell (i, k) is at X = (i + ½)/N, Z = (k + ½)/N; one lattice step is
/// Δτ = 1/N² (diffusive scaling, so ν̃ and ν are the same lattice numbers on every grid and the
/// time error, O(Δτ), falls as fast as the second-order space error); the lattice velocity is
/// a/N and the source per lattice step S/N².
fn run_box(m: &Manufactured, faces: [BoxFace; 6], speed: f64, n: usize) -> Errors {
    let cells = n * n;
    let at = |c: usize| (((c % n) as f64 + 0.5) / n as f64, ((c / n) as f64 + 0.5) / n as f64);
    let mut field = SaField::boxed([n, 1, n], &vec![true; cells], faces, m.f0, 0.0).unwrap();
    let initial: Vec<f64> = (0..cells)
        .map(|c| {
            let (x, z) = at(c);
            m.value(x, z, 0.0)
        })
        .collect();
    field.set_nu_tilde(&initial).unwrap();
    let u = vec![[speed / n as f64, 0.0, 0.0]; cells];
    let zeros = vec![0.0; cells];
    let steps = n * n / 8;
    let dt = 1.0 / (n * n) as f64;
    for s in 0..steps {
        let t = s as f64 * dt;
        let source: Vec<f64> = (0..cells)
            .map(|c| {
                let (x, z) = at(c);
                m.source(x, z, t, speed) * dt
            })
            .collect();
        let input = SaInput { u: &u, vorticity: &zeros, length: &zeros, nu: NU, model: false, source: Some(&source), substeps: None };
        field.advance(&input).unwrap();
    }
    let t = steps as f64 * dt;
    let (mut sum, mut max) = (0.0f64, 0.0f64);
    for (c, v) in field.nu_tilde().iter().enumerate() {
        let (x, z) = at(c);
        let e = v - m.value(x, z, t);
        sum += e * e;
        max = max.max(e.abs());
    }
    Errors { l2: (sum / cells as f64).sqrt(), linf: max }
}

const WALL: Manufactured = Manufactured { f0: 0.0, g: wall_g };
const INFLOW: Manufactured = Manufactured { f0: 0.02, g: inflow_g };
/// x periodic, a wall at z = 0 and an open face at z = 1 (y one periodic cell).
const WALL_FACES: [BoxFace; 6] = [BoxFace::Periodic, BoxFace::Periodic, BoxFace::Periodic, BoxFace::Periodic, BoxFace::Wall, BoxFace::Open];
/// An inflow face at x = 0, an open outlet at x = 1, z periodic.
const INFLOW_FACES: [BoxFace; 6] = [BoxFace::Inflow, BoxFace::Open, BoxFace::Periodic, BoxFace::Periodic, BoxFace::Periodic, BoxFace::Periodic];

/// The hand-written derivatives agree with central differences of G (h = 1e-4, so the check
/// holds to about 1e-6 of the derivative's scale).
#[test]
fn the_manufactured_derivatives_are_right() {
    let h = 1e-4;
    for g in [wall_g as fn(f64, f64) -> [f64; 5], inflow_g] {
        for &(x, z) in &[(0.1, 0.2), (0.37, 0.81), (0.9, 0.55), (0.5, 0.05)] {
            let [g0, gx, gz, gxx, gzz] = g(x, z);
            let numeric = [
                (g(x + h, z)[0] - g(x - h, z)[0]) / (2.0 * h),
                (g(x, z + h)[0] - g(x, z - h)[0]) / (2.0 * h),
                (g(x + h, z)[0] - 2.0 * g0 + g(x - h, z)[0]) / (h * h),
                (g(x, z + h)[0] - 2.0 * g0 + g(x, z - h)[0]) / (h * h),
            ];
            for (want, got) in numeric.iter().zip([gx, gz, gxx, gzz]) {
                assert!((want - got).abs() < 1e-5 * (1.0 + got.abs()), "({x}, {z}): {want} vs {got}");
            }
        }
    }
}

/// Advection, diffusion and the c_b2 term through a wall face: first order (upwind advection).
#[test]
fn sa_transport_next_to_a_wall_is_first_order() {
    let p = study("SaField, wall face, advection + diffusion + c_b2", &GRIDS, |n| run_box(&WALL, WALL_FACES, SPEED, n));
    let (l2, linf) = p[1];
    assert!((0.9..=1.2).contains(&l2) && (0.9..=1.2).contains(&linf), "{p:?}");
}

/// Diffusion and the c_b2 term alone through a wall face: second order.
#[test]
fn sa_diffusion_next_to_a_wall_is_second_order() {
    let p = study("SaField, wall face, diffusion + c_b2", &GRIDS, |n| run_box(&WALL, WALL_FACES, 0.0, n));
    let (l2, linf) = p[1];
    assert!((1.8..=2.2).contains(&l2) && (1.8..=2.2).contains(&linf), "{p:?}");
}

/// Advection, diffusion and the c_b2 term behind an inflow face (ν̃ = 0.02 there) to an open
/// outlet: first order.
#[test]
fn sa_transport_behind_an_inflow_is_first_order() {
    let p = study("SaField, inflow face, advection + diffusion + c_b2", &GRIDS, |n| run_box(&INFLOW, INFLOW_FACES, SPEED, n));
    let (l2, linf) = p[1];
    assert!((0.9..=1.2).contains(&l2) && (0.9..=1.2).contains(&linf), "{p:?}");
}

/// Diffusion and the c_b2 term alone behind an inflow face: second order.
#[test]
fn sa_diffusion_behind_an_inflow_is_second_order() {
    let p = study("SaField, inflow face, diffusion + c_b2", &GRIDS, |n| run_box(&INFLOW, INFLOW_FACES, 0.0, n));
    let (l2, linf) = p[1];
    assert!((1.8..=2.2).contains(&l2) && (1.8..=2.2).contains(&linf), "{p:?}");
}
