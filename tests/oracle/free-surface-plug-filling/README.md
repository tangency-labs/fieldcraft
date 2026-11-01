# A plug of liquid filling a frictionless channel from a velocity inlet

A channel 1 m long, one cell (Δx = 1/64 m) across and periodic across, starts full of gas with a
velocity inlet at x = 0 injecting liquid at U = 0.05 m/s (ν = 10⁻³ m²/s, ρ = 1). With no wall
friction the liquid is a plug at U: the volume behind the front is the injected volume plus the
inlet cell's, so the front lies at x_f = U t + Δx and the liquid's volume is (U t + Δx) Δx².

**Measures:** the liquid's volume at the end (t = 12 s: 1.50293·10⁻⁴ m³) from the free surface's
history (`solve.probe` with `front`), the time the liquid's extent (the face of the last cell at
least half full) reaches x = 0.5 m, which it does when the front x_f passes half a cell before it,
(0.5 − 1.5 Δx) / U = 9.53125 s, and the solve's relative mass error (zero to round-off).

**Run:** Δt = 1/64 s (lattice speed 0.05, τ = 0.69), 768 steps.
