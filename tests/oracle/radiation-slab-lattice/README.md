# A radiating slab in the lattice's conjugate heat transfer

A column of fluid 1 m tall between a hot plate at the bottom (y = 0, 373.15 K) and a cold plate
at the top (y = 1, 293.15 K) holds a conducting slab on the hot plate, 0.25 m thick. The slab's
faces that touch the fluid radiate to surroundings at 293.15 K (ε = 0.9) through the fluid, which
is transparent; nothing flows. In steady state the heat the hot plate gives the slab leaves its
top face two ways: conducted through the fluid to the cold plate, and radiated.

**Model:** one cell wide, periodic across x and z, 16 cells of Δx = 0.0625 m along y, in SI with
temperatures in kelvin. The fluid has ρ = 1, c = 1, k = 5 (α = 5 m²/s), ν = 0.1 m²/s; the slab
(`design.json`, a box reaching below the column so only its top face meets the fluid) has
ρ = 1, c = 1, k_s = 20. Δt = Δx² / (8 α), run to t = 8/α = 1.6 s (16 384 steps), long past the
column's diffusion time.

**Closed form:** the lattice holds the plates at the faces, half a cell from the first cell
centres, and couples the slab to the fluid across the face between their cells by the series
resistance of the two half cells, so its steady state is a thermal circuit with a radiating node,
the slab's top cell T_c:

(T₀ − T_c)/R₁ = (T_c − T₁)/R₂ + ε σ (T_c⁴ − T_sur⁴),
R₁ = (0.25 − Δx/2)/k_s, R₂ = Δx/(2 k_s) + 0.75/k.

Its root, by bisection to round-off, is T_c = 362.609625 K; the radiated power over the column's
cross-section is ε σ (T_c⁴ − T_sur⁴) Δx² = 1.97422 W, 52 % of the 3.76442 W the hot plate gives.

**Measures:** the slab's top-cell temperature (`solve.probe` at its centre) and the slab's
radiated power (the solve's `solids[0].radiated.final`). The discrete solution is the circuit
exactly, so 8 and 16 cells (`crates/lbm/tests/conjugate.rs`) both match to 10⁻⁹ K; the
tolerances, 10⁻⁸ relative for the temperature and 10⁻⁶ for the power, leave room for the
remaining transient and round-off.

Reference: F. P. Incropera, D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*, §3.1 (the
thermal circuit) and §1.3.1 (the surface energy balance with radiation to large surroundings).
