//! The model's own steady one-dimensional channel solution (see `sa_law_1d.rs`), shared by the tests that read it.

use fieldcraft_lbm::sa;

pub const NU: f64 = 1e-3;

/// The bulk C_f of the solution on the 1200-interval grid of `grid(1200, 0.05 * NU, 1.0)` (the trapezoid rule over the half channel).
pub const NEWTON_CF: f64 = 5.07708e-3;

/// The grid of `n` intervals, geometric from `y1` to `top`.
pub fn grid(n: usize, y1: f64, top: f64) -> Vec<f64> {
    let (mut lo, mut hi) = (1.0 + 1e-9, 1.5f64);
    for _ in 0..200 {
        let r = 0.5 * (lo + hi);
        if y1 * (r.powi(n as i32) - 1.0) / (r - 1.0) < top { lo = r } else { hi = r }
    }
    let r = 0.5 * (lo + hi);
    (0..=n).map(|i| if i == 0 { 0.0 } else { y1 * (r.powi(i as i32) - 1.0) / (r - 1.0) }).collect()
}

/// The steady equation's residual at node `i` (1 ≤ i ≤ n) of ν̃ on grid `y`: the finite-volume
/// diffusion, the c_b2 cross term (the gradient exact for quadratics on a stretched grid), P and D;
/// the total stress is `1 - y`; the top node has no flux through its upper face.
fn residual(i: usize, y: &[f64], nt: &[f64], nu: f64) -> f64 {
    let n = y.len() - 1;
    let (lo, mid) = (nt[i - 1], nt[i]);
    let hm = y[i] - y[i - 1];
    let dm = 0.5 * (sa::diffusivity(mid, nu) + sa::diffusivity(lo, nu)) / sa::SIGMA;
    let stress = (1.0 - y[i]).max(0.0);
    let shear = stress / (nu + sa::nu_t(mid, nu));
    let s = sa::sources(mid, nu, shear, y[i]);
    if i == n {
        let flux_in = dm * (mid - lo) / hm;
        return (-flux_in) / (0.5 * hm) + s.production - s.destruction;
    }
    let (hp, up) = (y[i + 1] - y[i], nt[i + 1]);
    let dp = 0.5 * (sa::diffusivity(mid, nu) + sa::diffusivity(up, nu)) / sa::SIGMA;
    let width = 0.5 * (hm + hp);
    let grad = (hm * hm * up - hp * hp * lo + (hp * hp - hm * hm) * mid) / (hm * hp * (hm + hp));
    (dp * (up - mid) / hp - dm * (mid - lo) / hm) / width + sa::C_B2 / sa::SIGMA * grad * grad + s.production - s.destruction
}

/// The steady problem as a boundary-value problem: Newton's method on the full residual with a
/// finite-difference tridiagonal Jacobian and a backtracking line search that keeps ν̃ ≥ 0, started from the algebraic layer
/// ν̃ = κ y. Returns (ν̃, y, the iterations, the largest residual relative to its scale κ²(1 + c_b2)/σ).
pub fn solve_newton(nu: f64, y: &[f64]) -> (Vec<f64>, Vec<f64>, usize, f64) {
    let n = y.len() - 1;
    let mut nt: Vec<f64> = (0..=n).map(|i| if i == 0 { 0.0 } else { sa::KAPPA * y[i] }).collect();
    let scale = sa::KAPPA * sa::KAPPA / sa::SIGMA * (1.0 + sa::C_B2);
    let res_all = |nt: &[f64]| -> Vec<f64> { (1..=n).map(|i| residual(i, y, nt, nu) / scale).collect() };
    let norm = |r: &[f64]| r.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let mut r = res_all(&nt);
    let mut iters = 0;
    for _ in 0..200 {
        iters += 1;
        if norm(&r) < 1e-12 {
            break;
        }
        // Tridiagonal Jacobian by central differences in each of the three arguments.
        let (mut a, mut b, mut c) = (vec![0.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1]);
        let mut p = nt.clone();
        for i in 1..=n {
            let mut local = [0.0; 3];
            for (q, j) in [i - 1, i, i + 1].into_iter().enumerate() {
                if j > n || (j == 0) {
                    continue;
                }
                let h = (1e-6 * nt[j].abs()).max(1e-9 * nu);
                p[j] = nt[j] + h;
                let rp = residual(i, y, &p, nu);
                p[j] = nt[j] - h;
                let rm = residual(i, y, &p, nu);
                p[j] = nt[j];
                local[q] = (rp - rm) / (2.0 * h) / scale;
            }
            a[i] = local[0];
            b[i] = local[1];
            c[i] = local[2];
        }
        // Solve J d = -r by the Thomas algorithm over i = 1..=n.
        let mut d: Vec<f64> = std::iter::once(0.0).chain(r.iter().map(|v| -v)).collect();
        let mut bb = b.clone();
        for i in 2..=n {
            let w = a[i] / bb[i - 1];
            bb[i] -= w * c[i - 1];
            d[i] -= w * d[i - 1];
        }
        let mut step = vec![0.0; n + 1];
        step[n] = d[n] / bb[n];
        for i in (1..n).rev() {
            step[i] = (d[i] - c[i] * step[i + 1]) / bb[i];
        }
        let mut lambda = 1.0;
        let before = norm(&r);
        loop {
            let trial: Vec<f64> = (0..=n).map(|i| (nt[i] + lambda * step[i]).max(0.0)).collect();
            let rt = res_all(&trial);
            if norm(&rt) < before || lambda < 1e-6 {
                nt = trial;
                r = rt;
                break;
            }
            lambda *= 0.5;
        }
    }
    (nt, y.to_vec(), iters, norm(&r))
}

/// U from ν̃ on the grid: the trapezoid rule on the shear.
pub fn velocity(y: &[f64], nt: &[f64], nu: f64) -> Vec<f64> {
    let shear: Vec<f64> = (0..y.len()).map(|i| (1.0 - y[i]).max(0.0) / (nu + sa::nu_t(nt[i], nu))).collect();
    let mut u = vec![0.0; y.len()];
    for i in 1..y.len() {
        u[i] = u[i - 1] + 0.5 * (shear[i] + shear[i - 1]) * (y[i] - y[i - 1]);
    }
    u
}
