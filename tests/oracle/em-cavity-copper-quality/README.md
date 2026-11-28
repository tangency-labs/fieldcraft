# Quality factor of a rectangular copper cavity

The 20 × 10 × 15 mm box of air of `em-cavity-rectangular`, its six walls copper (σ = 5.813e7 S/m,
Pozar appendix F) as a `conductor` boundary condition, meshed with 16 × 8 × 12 blocks of tetrahedra
and solved for its lowest mode, TE101, with first-order Nédélec edge elements (`solve.cavity`). The
walls absorb power through their surface resistance Rₛ = √(ω μ₀ / (2σ)) (Pozar eq. 1.131); the
closed form of the quality factor by the perturbation method (Pozar §6.7, eq. 6.46) is

Q = (k a d)³ b η / (2π² Rₛ) · 1 / (2a³b + 2bd³ + a³d + ad³), k = π √(1/a² + 1/d²), η = μ₀ c₀.

Refinement: the solver's tests (`crates/solver/tests/cavity.rs`) solve a box of the same proportions;
the error of Q falls as h² at both orders (order 1: 5.4 % on 4 × 2 × 3 cells, 1.4 % on 8 × 4 × 6,
0.65 % on 12 × 6 × 9), so on 16 × 8 × 12 it is about 0.35 %, inside the 1 % tolerance.
