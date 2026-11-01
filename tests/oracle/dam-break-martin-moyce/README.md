# Collapse of a liquid column (dam break) against Martin and Moyce

Martin and Moyce (1952) released columns of water from behind a gate onto a dry horizontal
plane and filmed the surge front and the falling column. Their Table 2 (rectangular section,
n² = 2, base a = 2¼ in = 57.15 mm, height n² a = 2a) gives the mean of nine records of the time
T = n t √(g/a) (n = √2) at which the front passed each station Z = z/a, z measured from the plane
of symmetry (Z = 1 is the column's original face); Table 6 gives the column's top H = h/(n² a)
against τ = t √(g/a). The runs were "the emptying of a rectangular corner, the motions being
similar because of plane symmetry", so the back wall stands for the plane of symmetry.

**Time origins.** The release could not be filmed, so Table 2's times are normalised to pass
Z = 1.44 at T = 1.19 in every record; the case therefore compares the front's travel *between*
stations: T(3.00) − T(1.44) = 1.32 and T(4.11) − T(1.44) = 2.14. Table 6 counts τ from the
release trigger and the motion began at or before τ = 0.16; the case compares H at τ = 2 (0.50)
and the fall between τ = 1 and τ = 2 (0.806 − 0.50 = 0.306, the first interpolated between the
rows 0.83 at τ = 0.93 and 0.78 at τ = 1.08).

**The model.** A tank 5a long and 3a high, one cell deep and periodic across (a plane problem);
the column a × 2a at the back wall; gravity 9.81 m/s²; ρ = 1000 kg/m³. The back wall is a
slip face: it is the plane of symmetry, and a no-slip wall would hold a film of the viscous
liquid (below). The floor and the far wall are no-slip. 50 cells per a (Δx = 1.143 mm,
250 × 1 × 150 cells), Δt = 3.116e-05 s (lattice speed 0.05 at the free-fall speed over the tank's
height), run to τ = 4 (T = 5.66) with frames at τ = 1, 2, 3, 4.

**Viscosity.** The two-relaxation-time lattice holds τ ≥ 0.505; this case runs at τ = 0.515,
a kinematic viscosity of 2.10e-04 m²/s (μ = 0.21 Pa·s, 200 times water's). The Reynolds number
on the front speed (1.3 m/s) and a is about 350 against 10⁵ in the experiment. Gravity and
inertia drive the collapse, so the column's fall is insensitive to it; the front, which runs
along the no-slip floor, is slowed by it.

**Measures.** The front from `solve.probe` `front` (the face of the last cell at least half
full along x, at every time step; arrival times at the stations by interpolation); the column's
height from `solve.probe` `depth` in the first column of cells (the fill levels summed), at
frames 1 and 2 (τ = 1 and 2); the relative mass error of the solve.

**Refinement** (τ = 0.515 at every resolution, so the physical viscosity halves with each
refinement):

| Quantity | Martin & Moyce | 25 cells per a | 50 (this case) | 100 |
|---|---|---|---|---|
| T(3.00) − T(1.44) | 1.32 | 1.569 (+18.9 %) | 1.489 (+12.8 %) | 1.450 (+9.9 %) |
| T(4.11) − T(1.44) | 2.14 | 2.487 (+16.2 %) | 2.358 (+10.2 %) | 2.288 (+6.9 %) |
| H(τ = 1) | 0.806 | 0.805 | 0.804 | 0.805 |
| H(τ = 2) | 0.50 | 0.502 | 0.502 | 0.501 |
| H(τ = 1) − H(τ = 2) | 0.306 | 0.303 | 0.302 | 0.303 |
| H(τ = 3) | 0.29 (between rows) | 0.314 | — | 0.304 |
| relative mass error | 0 | −5·10⁻¹³ | −1·10⁻¹² | −2·10⁻¹² |

The column's fall matches the measurements within 0.5 % at every resolution; the front's travel
times converge toward them from above as the viscosity falls. With a no-slip back wall at 50
cells per a the front is unchanged (1.502, 2.395) but the first column of cells keeps a film of
the viscous liquid: H(τ = 2) = 0.888 there, while half a column in (x = a/2) the depth is 0.467.
The front speed between Z = 1.44 and 3 is 1.05 a per unit T here against the measured 1.18
(their largest measured U = dZ/dT was 1.71, the shallow-water value 2).
