# Natural convection in a square cavity with a centred conducting body, Ra = 10⁵

The de Vahl Davis cavity (side H, the left wall held hot, the right wall cold, the top and bottom
adiabatic, Pr = 0.71, Ra = g β ΔT H³ / (ν α) = 10⁵) with a square body of side H/2 at its centre.
The body conducts heat at k_s / k_f = 0.2, 1 or 5 and generates none: conjugate heat transfer, the
body's temperature coupled to the air's across its four faces. The quantity is the mean Nusselt
number of the hot wall, on the cavity's side and the fluid's conductivity, for each ratio. The
steady state does not depend on the body's heat capacity, which is the fluid's here.

**Model:** as in `de-vahl-davis-ra1e5`, dimensionless in SI (H = 1 m, α = 1 m²/s, ν = 0.71 m²/s,
β = 1/K, ΔT = 1 K), 2D as a one-cell-thick periodic slab. The body is `design.json`: a box
through the slab's thickness, a conducting solid of the fluid domain with its own material per
ratio (density 1, specific heat 1, conductivity 0.2, 1 or 5).

**Measure:** the heat the lattice gives the hot wall over k ΔT times the wall's area
(`/heat/0/mean`, scale −1/Δx = −64), as for the de Vahl Davis cases.

**Run:** 64 × 64 cells, the time step of `de-vahl-davis-ra1e5`, until the same steady stop rule.

**Grid refinement** (`crates/lbm/tests/conjugate.rs`, the same cavity and the same scheme, run for
a fixed time with the heat averaged over its last tenth):

| k_s / k_f | House et al. | 48 × 48 | 64 × 64 | extrapolated (48, 64; second order) |
|---|---|---|---|---|
| 0.2 | 4.6257 | 4.6511 | 4.6356 | 4.6157 (−0.22 %) |
| 1 | 4.5061 | 4.5299 | 4.5151 | 4.4961 (−0.22 %) |
| 5 | 4.3221 | 4.3433 | 4.3297 | 4.3122 (−0.23 %) |

The grid converges from above at second order. The tolerance, 0.6 %, covers the 64-cell grid's
distance to its own extrapolated limit (about 0.42 %) and a margin; House et al. do not state
their finest grid's accuracy. With no flow (the conduction limit) the published values are
0.7069, 1 and 1.4138; `crates/lbm/tests/conjugate.rs` checks those every run.

Reference: J. M. House, C. Beckermann, T. F. Smith, *Effect of a centered conducting body on
natural convection heat transfer in an enclosure*, Numerical Heat Transfer Part A 18 (1990)
213–225, doi:10.1080/10407789008944791, Table 3 (the finest grid).
