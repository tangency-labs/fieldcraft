# The characteristic outlet's reflection of a plane pulse

The characteristic pressure outlet (`outflow: "characteristic"`, `Flow::characteristic_outlet`)
relaxes the face's pressure towards its reference at the rate K = σ c_s / L, after the linearised
characteristic (LODI) treatment of Poinsot and Lele (1992), with σ ≈ 0.25 after Rudy and Strikwerda
(1980). Selle, Nicoud and Poinsot (2004) derive its reflection coefficient for a normal wave,
R(ω) = −1 / (1 − 2iω/K), so |R| = K / √(K² + 4ω²). In time, that is a one-pole low-pass of rate
K/2 applied to the wave a fully reflecting end returns.

**The duct** (`solve.outlet_reflection`, `crates/lbm/src/reflection.rs`, lattice units): 2048 × 1 × 1
cells, periodic sides, a velocity inlet at rest at x_min and a pressure outlet at x_max, TRT
τ = 0.51. A right-going Gaussian pulse ρ′ = 10⁻⁴ exp(−(x − 1000)² / (2 · 40²)), u′ = c_s ρ′ / ρ₀,
passes a probe at cell 1300, reaches the outlet and comes back. The duct runs 7400 steps twice,
with the uniform outlet and with the characteristic outlet.

**The check, in time.** The uniform outlet's reflection, filtered by y[t] = a y[t−1] + (1 − a) u[t]
with a = e^(−K/2), predicts the characteristic outlet's. The case gates at σ 1.0 and at σ 0.25
(the default) on:

- the largest departure of the measured reflection from that prediction, over the incident peak,
  up to step 7268, before the inlet's reflection of the reflection returns;
- the decay rate fitted (least squares of ln |r| on t) to the reflection's tail over K/2;
- the control: the uniform outlet's reflected over incident peak (a fully reflecting end).

**What the departure is.** The largest departure is the same at the three measured σ values
(2.28 × 10⁻³ of the incident peak at σ 0.25, 0.5 and 1.0), so it does not come from the
relaxation. It is the
closure's own direct reflection: the outlet's floor, 0.23 % of the incident peak.

**Why not |R| in the relaxation band.** This pulse experiment cannot resolve |R| at ω ≈ K/4 to 4K
before the inlet reflection returns: with σ ≤ 1, the relaxation's time constant 2/K is at least
the round trip, so the inlet's reflection returns before the tail decays. The time-domain check above tests the same
one-pole response instead. The record holds 0.54 e-folds of the tail at σ 1.0 and 0.13 at σ 0.25.

**The windowed spectral ratio** at ω 0.01, 0.02 and 0.04 radians per step (σ 0.25) is the ratio
of the discrete Fourier amplitudes of the reflected and the incident pulse over rectangular windows
of ±5w/c_s = ±346 steps. The window cuts off the relaxation's slow tail (time constant
2/K = 28,000 steps at σ 0.25), so these are not the outlet's |R| at those frequencies: they
include the truncated relaxation tail, and at ω 0.04 the value is almost all truncation (it
scales with σ: 0.0426 at σ 0.25, 0.0817 at 0.5, 0.1579 at 1.0). They are compared with the same
window applied to the one-pole prediction and frozen as this runner's fixture, so that a change
to the outlet shows.

| Quantity | σ 0.25 | σ 0.5 (not gated) | σ 1.0 |
|---|---|---|---|
| K (per step) | 7.048 × 10⁻⁵ | 1.410 × 10⁻⁴ | 2.819 × 10⁻⁴ |
| Tail level over incident peak | 0.00605 | 0.01195 | 0.02333 |
| Largest departure over incident peak | 2.28 × 10⁻³ | 2.28 × 10⁻³ | 2.28 × 10⁻³ |
| Fitted tail rate over K/2 | 1.0054 | 1.0015 | 1.0005 |
| Uniform outlet, reflected over incident peak | 0.9928 | 0.9928 | 0.9928 |

**Cost:** each run takes about 2.4 s in release; the case makes four. It is marked slow.

## References
- L. Selle, F. Nicoud, T. Poinsot. *Actual impedance of nonreflecting boundary conditions:
  implications for computation of resonators.* AIAA J. 42(5):958–964, 2004.
- T. J. Poinsot, S. K. Lele. *Boundary conditions for direct simulations of compressible viscous
  flows.* J. Comput. Phys. 101:104–129, 1992.
- D. H. Rudy, J. C. Strikwerda. *A nonreflecting outflow boundary condition for subsonic
  Navier–Stokes calculations.* J. Comput. Phys. 36:55–70, 1980.
