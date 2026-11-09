# Heat through a thin wall between two held plates

Two plates a distance L = 1 m apart are held at 1 °C (y = 0) and 0 °C (y = 1 m), with still fluid
between them (conductivity k_f = 0.05 W/(m·K), ρ c = 1 J/(m³·K)). A thin wall at y = ½ divides the
fluid into two layers: a zero-thickness sheet with conductance G = 0.05 W/(m²·K) per area (its
conductivity over its thickness). In the steady state the heat crosses the two fluid layers and the
sheet in series, q = ΔT / (L/k_f + 1/G) = 1 / (20 + 20) = 0.025 W/m², and the sheet's two sides sit
at 1 − q (L/2)/k_f = 0.75 °C and q (L/2)/k_f = 0.25 °C: the jump across it is q/G.

**Model:** a column of 8 cells (Δx = 0.125 m) along y, one cell across in x and z with periodic
faces, the sheet on the face between the fourth and fifth cells (`baffles`), the plates as held
temperatures on y_min and y_max.

**Measures:** the solve's `baffles` report for the sheet: the heat flux into it from the hot side
and the mean temperatures of its two sides.

**Run:** Δt = Δx²/(8 α), 100 s (2560 steps), five diffusion times L²/α. The discrete steady
profile is linear in each layer and the sheet's exchange is Fourier's law across it, so the
lattice reproduces the series solution to round-off; `crates/lbm/tests/baffles.rs` checks the
same on all three axes and for an insulating and a near-perfect sheet.
