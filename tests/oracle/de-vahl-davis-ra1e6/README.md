# Natural convection in a square cavity (de Vahl Davis), Ra = 10⁶

The differentially heated square cavity: side H, the left wall (x = 0) held hot, the right wall
(x = H) held cold, the top and bottom adiabatic, no-slip everywhere, air (Pr = 0.71) under
Boussinesq buoyancy at a Rayleigh number g β ΔT H³ / (ν α) of 10⁶. The model is dimensionless
in SI: H = 1 m, α = 1 m²/s (k = 1 W/(m K), ρ = 1 kg/m³, c = 1 J/(kg K)), ν = 0.71 m²/s, β = 1/K,
ΔT = 1 K (1 °C hot, 0 °C cold, 0.5 °C at the start, which is the buoyancy's reference), and g =
Ra ν α / (β ΔT H³). Velocities are then in units of α/H and the local Nusselt number on the hot wall
is minus its heat flux (the heat into the wall). The reference's z (vertical) is the model's y.

**2D in a 3D solver:** the domain is one cell thick (H × H × Δx) and periodic across its thickness,
the existing convention of the lattice tests: every column across the thickness is the same, so
the flow is exactly two-dimensional.

**Measures:** the mean Nusselt number is the heat the lattice gives the hot wall (its anti-bounce-back
links, conservative) over k ΔT times the wall's area. The profiles (`solve.probe` with a `line`)
are read at the cell centres along the hot wall (`wall_heat_flux`: the one-sided second-order
difference through the held wall temperature and the two cell layers beside it), along the vertical
mid-plane (u) and the horizontal mid-plane (w); the extremes are located between the samples by
the parabola through the extreme sample and its neighbours, as de Vahl Davis located his by
polynomial interpolation. The smallest Nusselt number is at the top corner (z = H) in the
reference; ours is read half a cell below it, at the top cell's centre.

**Run:** 128 × 128 cells, time step for a lattice speed of at most 0.05 at de Vahl Davis's largest
velocity and a thermal relaxation time τ_T ≤ 0.8, until the steady stop rule (the velocity and
temperature residuals and the heat into both walls changing by less than 10⁻⁵ over two windows of
0.02 H²/α) ends it: 123552 time steps, which take several minutes in `cargo xtask oracle` on an
Apple M4 Max shared with other builds (103–134 s, 431 s and 927 s in runs at load averages of
about 40, 30 and 90: the time follows what else the machine runs). The case therefore runs nightly (`"nightly": true`,
`cargo xtask oracle --nightly`); other runs check it by its fixture.

**Grid refinement** (the tolerance of each quantity is de Vahl Davis's stated accuracy at this Ra,
plus the change between the two grids, plus, for Nu_min, the profile's change over the half cell
to the corner; positions take, instead of his accuracy, half of the finest mesh he solved at this
Ra (h = 0.025 for Ra ≤ 10⁴, 0.0125 above), since his positions come from polynomial interpolation
and his bench mark itself is extrapolated on a coarser 21 × 21 or 41 × 41 subset):

| Quantity | de Vahl Davis | 128 × 128 (this case) | 256 × 256 |
|---|---|---|---|
| Nu (mean, hot wall) | 8.8 | 8.8476 | 8.8258 |
| Nu_max | 17.925 | 17.44 | 17.546 |
| z of Nu_max | 0.0378 | 0.040902 | 0.039275 |
| Nu_min | 0.989 | 1.0423 | 0.94515 (smallest sample) |
| z of Nu_min | 1 | 0.98953 | 0.99805 (its cell) |
| u_max on x = H/2 | 64.63 | 65.46 | 65.184 |
| z of u_max | 0.85 | 0.85189 | 0.85166 |
| w_max on z = H/2 | 219.36 | 221.36 | 220.84 |
| x of w_max | 0.0379 | 0.038066 | 0.037863 |

The largest local Nusselt number misses de Vahl Davis's value: 2.7 % low on 128 cells and 2.1 %
on 256, against a tolerance of 1.6 %. Here the benchmark is the outlier. De Vahl Davis's 17.925
is the one value he could not extrapolate like the others: the usual pair of grids gave a
maximum above the finer grid's, so he extrapolated it from his 61 × 61 and 81 × 81 solutions
instead and located it by numerical differentiation (p. 255, Fig. 2). Later reference solutions give
17.536 with a mean Nusselt number of 8.82513 (Hortmann, Perić & Scheuerer, *Finite volume
multigrid prediction of laminar natural convection: bench-mark solutions*, Int. J. Numer. Methods
Fluids 11 (1990) 189–207, doi:10.1002/fld.1650110206; the values as quoted in Table 1 of the
ERCOFTAC knowledge base, UFR 4-10, a secondary source, as the paper was not read). Le Quéré
(*Accurate solutions to the square thermally driven cavity at high Rayleigh number*, Computers &
Fluids 20 (1991) 29–41, doi:10.1016/0045-7930(91)90025-d) is the spectral solution they agree with
(its values were not read here). Our 256 × 256 run gives 17.546 at z = 0.0393, 0.06 % from
Hortmann's, and a mean Nusselt number of 8.8258 (0.007 % from it). De Vahl Davis stays the
reference and the tolerance stays; the quantity is a known failure against it.

The thermal lattice's walls do not depend on the time step (#629): it has two relaxation times with
a fixed Λ_T = 1/12 from τ_T = 0.79 (docs/manual/fluid/collision.md). With the single relaxation time it had
before, at Ra = 10³ on 32 × 32 cells the largest local Nusselt number was 1.635 (8.6 % high) at τ_T
= 1.5 and 1.512 (0.5 %) at τ_T = 0.8; now it is 1.513 (+0.56 %) and 1.511 (+0.41 %), and the
smallest 0.3 % low at both (`local_wall_heat_flux_with_two_relaxation_times` in
`crates/lbm/tests/thermal.rs`). This case runs at τ_T = 0.62, below 0.79, where the lattice keeps
the single relaxation time (Λ_T = (τ_T − ½)², smaller than 1/12 there and as stable as before), so
its results did not change.

Reference: G. de Vahl Davis, *Natural convection of air in a square cavity: a bench mark numerical
solution*, International Journal for Numerical Methods in Fluids 3 (1983) 249–264,
doi:10.1002/fld.1650030305 (confirmed on Crossref), Table V, "The bench mark solution", p. 256;
its stated accuracy is better than 1 % at Ra = 10⁶ and down to 0.1 % at Ra = 10³ (taken as
0.1, 0.22, 0.46 and 1 % for Ra = 10³ … 10⁶, geometrically between his two stated ends).
