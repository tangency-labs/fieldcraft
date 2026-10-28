# Natural convection in a square cavity (de Vahl Davis), Ra = 10⁵

The differentially heated square cavity: side H, the left wall (x = 0) held hot, the right wall
(x = H) held cold, the top and bottom adiabatic, no-slip everywhere, air (Pr = 0.71) under
Boussinesq buoyancy at a Rayleigh number g β ΔT H³ / (ν α) of 10⁵. The model is dimensionless
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

**Run:** 64 × 64 cells, time step for a lattice speed of at most 0.05 at de Vahl Davis's largest
velocity and a thermal relaxation time τ_T ≤ 0.8, until the steady stop rule (the velocity and
temperature residuals and the heat into both walls changing by less than 10⁻⁵ over two windows of
0.02 H²/α) ends it: 26280 time steps, 10 s in `cargo xtask oracle` on an Apple M4 Max shared with other
builds (load average about 40).

**Grid refinement** (the tolerance of each quantity is de Vahl Davis's stated accuracy at this Ra,
plus the change between the two grids, plus, for Nu_min, the profile's change over the half cell
to the corner; positions take, instead of his accuracy, half of the finest mesh he solved at this
Ra (h = 0.025 for Ra ≤ 10⁴, 0.0125 above), since his positions come from polynomial interpolation
and his bench mark itself is extrapolated on a coarser 21 × 21 or 41 × 41 subset):

| Quantity | de Vahl Davis | 64 × 64 (this case) | 128 × 128 |
|---|---|---|---|
| Nu (mean, hot wall) | 4.519 | 4.5346 | 4.5228 |
| Nu_max | 7.717 | 7.7371 | 7.7534 |
| z of Nu_max | 0.081 | 0.082661 | 0.081089 |
| Nu_min | 0.729 | 0.7498 | 0.70909 |
| z of Nu_min | 1 | 0.99219 | 0.99609 |
| u_max on x = H/2 | 34.73 | 34.978 | 34.835 |
| z of u_max | 0.855 | 0.85501 | 0.85496 |
| w_max on z = H/2 | 68.59 | 68.687 | 68.665 |
| x of w_max | 0.066 | 0.066778 | 0.066115 |

The thermal lattice's walls do not depend on the time step (#629): it has two relaxation times with
a fixed Λ_T = 1/12 from τ_T = 0.79 (docs/manual/fluid/collision.md). With the single relaxation time it had
before, at Ra = 10³ on 32 × 32 cells the largest local Nusselt number was 1.635 (8.6 % high) at τ_T
= 1.5 and 1.512 (0.5 %) at τ_T = 0.8; now it is 1.513 (+0.56 %) and 1.511 (+0.41 %), and the
smallest 0.3 % low at both (`local_wall_heat_flux_with_two_relaxation_times` in
`crates/lbm/tests/thermal.rs`). This case runs at τ_T = 0.69, below 0.79, where the lattice keeps
the single relaxation time (Λ_T = (τ_T − ½)², smaller than 1/12 there and as stable as before), so
its results did not change.

Reference: G. de Vahl Davis, *Natural convection of air in a square cavity: a bench mark numerical
solution*, International Journal for Numerical Methods in Fluids 3 (1983) 249–264,
doi:10.1002/fld.1650030305 (confirmed on Crossref), Table V, "The bench mark solution", p. 256;
its stated accuracy is better than 1 % at Ra = 10⁶ and down to 0.1 % at Ra = 10³ (taken as
0.1, 0.22, 0.46 and 1 % for Ra = 10³ … 10⁶, geometrically between his two stated ends).
