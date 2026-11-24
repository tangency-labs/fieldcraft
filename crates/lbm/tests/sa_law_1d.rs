//! The wall law beside the model's own one-dimensional channel solution: Spalart–Allmaras on a stretched grid to the wall (y⁺ 0.05 at the
//! first node, 1200 intervals), steady, at Re_τ 1000 (ν = 1e-3, u_τ = δ = 1), the total stress (ν + ν_t) U' = 1 − y, solved by Newton's
//! method on the full discrete equation (`support/sa_1d.rs`). The constant-stress law ([`fieldcraft_lbm::sa_law::u_plus`]) is within 0.1 of
//! the solution in u⁺ up to y⁺ 60 (the stress and ν̃ taper of the channel drift it away above that), and the solution's bulk C_f is the
//! reference the wall-layer tests of `sa_wall_law.rs` measure against.

#[path = "support/sa_1d.rs"]
#[allow(dead_code)]
mod sa_1d;

use fieldcraft_lbm::sa;
use fieldcraft_lbm::sa_law;
use sa_1d::NU;

#[test]
fn the_wall_law_and_the_bulk_friction_beside_the_one_dimensional_sa_channel() {
    let y = sa_1d::grid(1200, 0.05 * NU, 1.0);
    let (nt, y, iterations, residual) = sa_1d::solve_newton(NU, &y);
    assert!(residual < 1e-9, "Newton did not converge: {residual:e} after {iterations} iterations");
    let u = sa_1d::velocity(&y, &nt, NU);
    let at = |yp: f64, f: &[f64]| {
        let target = yp * NU;
        let k = y.partition_point(|v| *v < target).clamp(1, y.len() - 1);
        f[k - 1] + (f[k] - f[k - 1]) * (target - y[k - 1]) / (y[k] - y[k - 1])
    };
    for yp in [15.0, 30.0, 60.0] {
        let (law, one) = (sa_law::u_plus(yp), at(yp, &u));
        assert!((law - one).abs() < 0.1, "the law is {:+.3} from the solution at y+ {yp}", law - one);
        // ν̃ is κ u_τ y in the layer, to a few percent, which is what the law assumes.
        let ratio = at(yp, &nt) / (sa::KAPPA * yp * NU);
        assert!((ratio - 1.0).abs() < 0.1, "ν̃ / (κ y) is {ratio:.3} at y+ {yp}");
    }
    let bulk: f64 = (1..y.len()).map(|i| 0.5 * (u[i] + u[i - 1]) * (y[i] - y[i - 1])).sum();
    let cf = 2.0 / (bulk * bulk);
    assert!((cf / sa_1d::NEWTON_CF - 1.0).abs() < 1e-4, "C_f {cf:.5e}");
}
