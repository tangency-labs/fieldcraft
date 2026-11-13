# Hagen–Poiseuille flow in a round pipe with interpolated walls

A pipe of radius R = 10 mm along z, the fluid driven by a body force a = 9.6 m/s² (ρ = 1 kg/m³,
μ = 10⁻³ Pa s, so ν = 10⁻³ m²/s), periodic along the pipe and one cell deep. Fully developed flow is
u(r) = a (R² − r²) / (4 ν), u_max = a R² / (4 ν) = 0.24 m/s on the axis.

**Geometry:** `design.json` cuts a cylinder of radius 10 mm out of a 30 × 30 × 10 mm block; the
block is the obstacle, its cylindrical face the pipe wall. The face's exact cylinder gives the link
fractions q (M. Bouzidi, M. Firdaouss, P. Lallemand, *Momentum transfer of a Boltzmann-lattice
fluid with boundaries*, Phys. Fluids 13(11):3452–3459, 2001, doi:10.1063/1.1399290).

**Run:** `walls: "interpolated"`, two-relaxation-time collision (Λ = 3/16), Δt = Δx² / (6 ν) so
τ = 1 on both grids, for 0.5 s: five viscous times R²/ν. The boxes are R + 2.5 Δx on each side of
the axis, so a cell centre lies on it; the measure is the largest velocity on a line of cell
centres across the pipe (`/line/max/velocity`).

**Grid refinement** (macOS aarch64):

| Walls | 8 cells per radius | 16 cells per radius |
|---|---|---|
| interpolated (this case) | 0.23963 (−0.15 %) | 0.23990 (−0.04 %) |
| halfway (the same recipe with `walls: "halfway"`) | 0.22678 (−5.5 %) | 0.23503 (−2.1 %) |

Interpolated walls converge at about second order (the error falls 3.7 times), halfway walls at
about first (2.7 times). The crate's test `pipe_poiseuille_converges_at_second_order_with_interpolated_walls`
(`crates/lbm/tests/interpolated.rs`) measures the L2 error over 4 to 32 cells per radius.
