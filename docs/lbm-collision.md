# Lattice Boltzmann collision: two relaxation times or regularised

A fluid step picks its collision operator with `collision` (`model.step.create`, procedure
`fluid`; the step form's *Collision* field):

| `collision` | What it is | When |
|---|---|---|
| `"trt"` (default) | Two relaxation times (Ginzburg, Verhaeghe & d'Humières 2008) with Λ = 3/16 | Most flows. Bounce-back walls sit exactly halfway between cells, so it is the more accurate of the two on a given grid. |
| `"regularized"` | Recursive regularised (Latt & Chopard 2006; Coreixas et al. 2017) | Flows the two-relaxation-time scheme cannot hold on the grid you can afford: τ⁺ below about 0.52, fast inlets and fans, sharp corners. It costs 15–35 % more per step, its staircase walls are less accurate and its pressure outlets cost 3–8 % near the outlet (below). |

The same choice runs on the CPU and on the GPU (`backend: "gpu"`), with fans, porous zones,
buoyancy and heat. `solve.fluid` reports it as `collision`. With the two-relaxation-time scheme
below τ⁺ = 0.51 the solve warns and suggests the regularised collision; with the regularised
collision it warns below τ⁺ = 0.508, its measured limit on the flow-past sample, where the
reflecting inlet and outlet still fail. Setting `magic` with the regularised collision is pointed
out, since Λ means nothing to it.

## The scheme

`crates/lbm/src/regularized.rs` has the derivation; in short:

- Every population is rebuilt from the cell's density, its velocity `u = (Σ f c + F/2)/ρ`
  (Guo's half force) and its non-equilibrium momentum flux before relaxing. The ghost modes the
  two-relaxation-time scheme relaxes by `1 − ω ≈ −1` as τ⁺ → ½ are set to equilibrium every
  step instead.
- Second order: the equilibrium and the flux `P = Σ c c f − ρ (c_s² I + u u) + ½ (u F + F u)`.
  Third order, recursively from `P`: `u_a P_bc + u_b P_ac + u_c P_ab`, plus the equilibrium
  `ρ u_a u_b u_c`, for the six third moments D3Q19 holds (`xxy, yzz, xyy, xzz, xxz, yyz`). On
  D3Q19 each pair (`H_xxy`, `H_yzz`) enters by its sum and difference, which are orthogonal
  under the weights (norms 2/27 and 6/27), so the populations carry exactly those moments.
- Guo forcing, re-derived for the rebuilt populations:
  `f* = f^eq + (1 − ω) f^(1)[P] + ½ F_q`, `F_q = w_q [3 c·F + 9 (c·u)(c·F) − 3 u·F]`. Mass is
  kept, the momentum gains exactly `F`, and fans, porous zones (Guo & Zhao's implicit velocity)
  and buoyancy enter through `F` and `u` as before.
- Pressure outlets hold the equilibrium of the outlet density and the normal velocity Zou–He
  finds (`hold_equilibrium` in `flow.rs`, and the same in the GPU shader). With the bounce-back
  non-equilibrium part a regularised outlet layer grows an oscillation that flips sign every
  step along the outlet (a 16-cell channel at τ⁺ = 0.52 and lattice speed 0.1 went unstable at
  step 301; `a_regularised_channel_with_a_pressure_outlet_stays_stable_at_tau_0_52` runs it 2000
  steps), which the two-relaxation-time scheme's slowly relaxing antisymmetric part damps. The
  equilibrium outlet (Krüger et al. 2017, ch. 5) does not impose the viscous stress on its own
  layer, and that costs accuracy near it, more as τ⁺ approaches ½: a fully developed channel
  profile is 3.2 % off in the cell next to the outlet at τ⁺ = 0.6 and 7–7.5 % at 0.52 (9.6 % in the
  outlet cell itself), against 0.4 % for the bounce-back outlet where that one stays stable
  (τ⁺ = 0.52, lattice speed 0.01). Mid-channel is unchanged. A wake leaving through the outlet is
  disturbed further in too: the pressure-outlet guard of #598
  (`a_wake_at_a_pressure_outlet_distorts_only_the_last_layer`, run with two relaxation times)
  measures 1.2 × 10⁻² interior departure under the regularised collision, beyond its 9 × 10⁻³
  limit. Keep outlets well downstream of what you measure; once the convective outflow
  (`outflow: convective`, #626) lands, prefer it with this collision. Velocity inlets keep Zou–He.
- Collisions run by fixed blocks of cells as before: 1 and N threads give the same result bit
  for bit (`results_are_bit_for_bit_the_same_on_one_and_many_threads`).

**Why this operator.** The cumulant collision (Geier et al. 2015) is defined on D3Q27's 27
moments; D3Q19 lacks the corner directions its higher cumulants need. The projected
regularisation (second order only) filters the ghost modes too, but drops the third moments; in
our runs it failed where the recursive form ran (a 64² lid-driven cavity at Re 2000, τ⁺ = 0.510:
projected unstable at step 786, recursive 20 000 steps; a periodic shear layer at lattice speed
0.2, τ⁺ = 0.505: projected unstable at step 603, recursive 5 000 steps), as Coreixas et al. (2017,
2019) found.

## Verification

Analytic and published references, `cargo test -p fieldcraft-lbm` (heavy ones with
`--include-ignored`); the two-relaxation-time results are unchanged.

| Case | Reference | Two relaxation times | Regularised |
|---|---|---|---|
| Plane Poiseuille, body force, 8/16/32 cells, τ = 1 | parabola | exact (round-off) | 5.2e-3, 1.3e-3, 3.3e-4 of the peak: order 2.00 |
| the same, τ = 0.6 | parabola | exact | 1.1e-2, 2.9e-3, 7.2e-4: order 2.00 |
| Couette, 16 cells | linear, wall shear μU/H | 1e-10, shear to 1e-6 | 2.3e-13, shear to 1e-6 |
| Taylor–Green decay, 16/32/64 cells | exp(−2νt) | 1.6e-2, 4.0e-3, 9.9e-4: order 2.0 | 1.7e-2, 4.3e-3, 1.1e-3: order 2.0; mass to 1e-13 |
| Sphere array, bounce-back, 1500 steps | mass, force balance | — | mass drift 6.7e-13; wall force = driving force to 1e-3 |
| Parabolic inlet, pressure outlet, 16 cells | Poiseuille; inflow = outflow | profile 3.0e-4, gradient +0.07 %, mass steady to 1e-12 | profile 2.3e-3, gradient +0.49 %, mass steady to 1e-12; 3.2e-2 next to the outlet |
| Porous plate (Darcy–Forchheimer, loss coefficient) | the law | within 0.5 % | within 0.5 % |
| Fan against a plate | operating point | flow 0.5 %, rise 1 %, mass 1e-12 | the same |
| de Vahl Davis cavity, Ra 10³, 16 / 32 cells | Nu 1.118 | 1.1161 / 1.1172 | 1.1173 / 1.1178 |
| de Vahl Davis cavity, Ra 10⁴, 48 cells | Nu 2.243 | 2.2449 | 2.2426 |
| Lid-driven cavity, Re 100, 64 cells | Ghia et al.: u/U −0.2109 at y 0.4531 | −0.2136 at 0.4585 | −0.2139 at 0.4586 |
| Schäfer–Turek 2D-1, Re 20, D/Δx 10 | c_D 5.57–5.59, c_L 0.0104–0.0110 | 5.989, 0.0120 | 6.303, 0.0180 |
| the same, D/Δx 20 | | 5.773 (+3.5 %), 0.0107 | 5.890 (+5.5 %), 0.0132 |
| Schäfer–Turek 2D-2, Re 100, D/Δx 20 | c_D,max 3.22–3.24, c_L,max 0.99–1.01, St 0.295–0.305 | 3.447, 1.102, 0.2987 | 3.496, 1.081, 0.2975 |

The regularised collision's walls are the price: without Λ = 3/16 its bounce-back wall is off
the halfway point by O(Δx²) (second order, but a larger constant), which shows on staircase
cylinders as a few percent more drag on coarse grids. Use it where the two-relaxation-time
scheme does not run, not in place of it.

**CPU and GPU** (`cargo test -p fieldcraft-lbm-gpu --test verification`, every case with both
collisions and both streaming layouts, Apple M4 Max, Metal; historical verification before
native GPU support became Vulkan-only): regularised velocity fields within
5e-7–7.4e-6 of the largest speed, densities within 2.6e-7, wall shear within 1.9e-5 of itself,
heats within 2.1e-6, mass drift 0: inside the backend's stated tolerances (2e-5, 1e-6, 5e-5, 1e-4).

## The thermal lattice (#629)

With heat, the D3Q7 temperature lattice (`flow.rs`, `lbm.wgsl`) collides with two relaxation
times too, whatever `collision` the flow uses: each axis pair's antisymmetric part, which
carries the heat flux, relaxes with `τ_T = 4 α + ½` (that sets the diffusivity), the symmetric
part and the rest population with `τ_T⁺ = ½ + Λ_T / (τ_T − ½)`. With `Λ_T` fixed a steady
temperature field does not depend on `τ_T` (d'Humières & Ginzburg 2009), and neither do the
walls: the anti-bounce-back of a held temperature, the bounce-back of an adiabatic wall or a flux,
and the interface links of conducting solids keep their rules, and the conjugate links still give
the solid exactly the heat the fluid loses (`conjugate.rs`). Before, the lattice had a single
relaxation time, `Λ_T = (τ_T − ½)²`, so the walls moved with the time step: on the automatic time
step's τ_T ≈ 1.44 the de Vahl Davis cavity's largest local Nusselt number was 7.5 % high and its
smallest 11.8 % low on 32 cells.

**Λ_T = 1/12** (`setup::THERMAL_MAGIC`). For a conducting field along a flat held wall (a harmonic
temperature, `T = Re e^{i k y − k x}`, solved to steady state on the lattice in a standalone
calculation of the scheme), the heat each wall link takes is fourth-order accurate at 1/12 —
8.1e-4, 5.0e-5 and 3.2e-6 of the flux on 8, 16 and 32 cells per wavelength — and second order at
every other value (6.4e-3 on 16 cells at 1/24 and 1/8, 1.3e-2 at 1/6, 2.5e-2 at 1/4); with a flow
along the wall the optimum moves to smaller values (1/24 at a cell Péclet number of 0.3). In
the cavity the wall-flux error grows with Λ_T, and 1/12 puts it within 0.2 % of the reference on
64 cells:

| Ra | cells | Λ_T | τ_T | Nu | Nu_max | Nu_min | Nu_max (links) | Nu_min (links) |
|---|---|---|---|---|---|---|---|---|
| 10³ | 32² | (τ_T − ½)², before | 0.6 | +0.02 % | −0.19 % | +0.71 % | +0.45 % | +0.28 % |
| 10³ | 32² | (τ_T − ½)², before | 0.8 | +0.01 % | +0.46 % | −0.40 % | +0.24 % | −0.39 % |
| 10³ | 32² | (τ_T − ½)², before | 1.0 | −0.00 % | +1.78 % | −2.66 % | +1.12 % | −2.26 % |
| 10³ | 32² | (τ_T − ½)², before | 1.5 | −0.07 % | +8.63 % | −13.51 % | +4.61 % | −9.66 % |
| 10³ | 32² | (τ_T − ½)², before | 1.439 (automatic) | −0.06 % | +7.54 % | −11.79 % | +4.15 % | −8.61 % |
| 10³ | 32² | 1/24 | 0.6 / 1.5 | +0.02 / −0.07 % | +0.06 / +0.23 % | +0.28 / +0.26 % | +0.35 / +0.52 % | +0.12 / +0.08 % |
| 10³ | 32² | 1/12 | 0.6 / 1.5 | +0.02 / −0.07 % | +0.40 / +0.56 % | −0.31 / −0.32 % | +0.24 / +0.41 % | −0.31 / −0.32 % |
| 10³ | 32² | 1/8 | 0.6 / 1.5 | +0.02 / −0.07 % | +0.73 / +0.89 % | −0.89 / −0.90 % | +0.34 / +0.55 % | −0.81 / −0.82 % |
| 10³ | 32² | 1/6 | 0.6 / 1.5 | +0.02 / −0.07 % | +1.06 / +1.22 % | −1.48 / −1.49 % | +0.59 / +0.80 % | −1.30 / −1.31 % |
| 10³ | 32² | 3/16 | 0.6 / 1.5 | +0.02 / −0.07 % | +1.24 / +1.39 % | −1.77 / −1.78 % | | |
| 10³ | 32² | 1/4 | 0.6 / 1.5 | +0.02 / −0.07 % | +1.74 / +1.88 % | −2.66 / −2.67 % | +1.07 / +1.28 % | −2.26 / −2.27 % |
| 10³ | 32² | shipped | 0.8 | +0.01 % | +0.41 % | −0.31 % | +0.26 % | −0.31 % |
| 10³ | 32² | shipped | 1.0 | −0.00 % | +0.44 % | −0.31 % | +0.28 % | −0.31 % |
| 10³ | 32² | shipped | 1.5 | −0.07 % | +0.56 % | −0.32 % | +0.41 % | −0.32 % |
| 10³ | 32² | shipped | 1.439 (automatic) | −0.06 % | +0.54 % | −0.31 % | +0.39 % | −0.32 % |
| 10³ | 64² | (τ_T − ½)², before | 0.8 | −0.01 % | +0.17 % | −0.17 % | +0.13 % | −0.18 % |
| 10³ | 64² | (τ_T − ½)², before | 1.5 | −0.03 % | +2.09 % | −3.42 % | +1.62 % | −3.28 % |
| 10³ | 64² | 1/24 | 1.5 | −0.03 % | +0.11 % | −0.01 % | +0.20 % | −0.01 % |
| 10³ | 64² | 1/12 | 1.5 | −0.03 % | +0.19 % | −0.15 % | +0.17 % | −0.16 % |
| 10³ | 64² | 1/6 | 1.5 | −0.03 % | +0.37 % | −0.44 % | +0.12 % | −0.45 % |
| 10³ | 64² | 1/4 | 1.5 | −0.03 % | +0.54 % | −0.73 % | +0.27 % | −0.75 % |
| 10³ | 64² | shipped | 0.8 | −0.01 % | +0.16 % | −0.15 % | +0.13 % | −0.16 % |
| 10⁵ | 64² | (τ_T − ½)², before | 0.8 | +0.16 % | +2.26 % | −11.06 % | +0.87 % | −10.96 % |
| 10⁵ | 64² | 1/24 | 0.8 | +0.16 % | +0.45 % | −4.03 % | +1.26 % | −6.66 % |
| 10⁵ | 64² | 1/12 (shipped) | 0.8 | +0.16 % | +2.01 % | −10.09 % | +0.92 % | −10.12 % |
| 10⁵ | 64² | 1/6 | 0.8 | +0.16 % | +5.22 % | −22.19 % | +0.27 % | −20.53 % |
| 10⁵ | 64² | 1/4 | 0.8 | +0.16 % | +8.48 % | −34.22 % | −0.35 % | −30.76 % |

Errors against de Vahl Davis (1983): Nu 1.118, Nu_max 1.505, Nu_min 0.692 at Ra 10³; 4.519,
7.717, 0.729 at 10⁵. Nu_max and Nu_min as `solve.probe` reads them (the one-sided difference
through the wall, the extreme on the parabola through three samples; Nu_min at the top cell's
centre, half a cell below the corner where the reference is) and from each wall link's own heat.
Run to 1.5 H²/α with `cargo test --release -p fieldcraft-lbm --test thermal -- --ignored
local_wall_heat_flux_with_two_relaxation_times --nocapture` (the 1/8 and 3/16 rows from the same
harness with those values). At a fixed Λ_T the rows move with τ_T by at most 0.16 points on 32
cells and 0.04 on 64 (the single relaxation time: 14 and 3.4): that remainder is the linear
equilibrium's numerical diffusion along the flow, 4 u² of the diffusivity in lattice units, which
grows with the time step's lattice speed (0.03 at τ_T 1.5 on 32 cells). At Ra 10⁵ on 64 cells the
top corner is under-resolved and smaller Λ_T reads closer there; 1/12 is the value that is exact
for conduction and keeps the stability of the single relaxation time (below), and the cavity
cases of the oracle, run at τ_T 0.62–0.8, stay within their tolerances.

**Below τ_T = 0.79** (`setup::thermal_magic`) the lattice keeps the single relaxation time
(`Λ_T = (τ_T − ½)²`, which is smaller there). A fixed Λ_T makes `τ_T⁺` grow without bound as
`τ_T → ½` (a fluid of high Prandtl number, such as water, on a fine grid), and a von Neumann
analysis of the D3Q7 lattice with its linear equilibrium finds that unstable from lattice speeds
of 0.02–0.03 at τ_T 0.505–0.51 for any Λ_T below 1/4, where the single relaxation time holds 0.25.
With `τ_T⁺ ≤ τ_T` the largest stable lattice speed is at least 0.25 at every τ_T, as with the
single relaxation time (0.25 at 0.505, 0.36 at 0.7, 0.44 at 0.8, 0.46–0.47 from 1 to 3, against
0.25, 0.36, 0.45 and 0.50–0.54), well above the characteristic lattice speeds a run is allowed
(Mach 0.3, 0.17). The single relaxation time's wall error grows
with its Λ_T, so below 0.79 it is smaller than at 1/12.

**What it costs.** A transient Fourier mode at rest decays 0.105 (k Δx)² faster than α k² (4.0e-3
on 32 cells per wavelength; the single relaxation time had −8e-4 at τ_T = 1 and −1.5e-2 at 1.5),
and the conjugate array of insulating cylinders on 16, 32 and 64 cells (`conjugate.rs`), a
first-order staircase interface, reads −5.7, −3.0 and −1.2 % against −3.9, −2.2 and −0.85 % at
the single relaxation time's τ_T = 1 (Λ_T = 1/4) and −12.4 and −4.7 % at τ_T = 0.6 (Λ_T = 0.01):
its error now does not depend on the time step.

**CPU and GPU.** The shader collides the same way (`collide_thermal`, the symmetric rate in
`therm.w`).

## Stability

**The issue's measurement** (#557): the flow-past sample (the bracket in a 240 × 120 × 100 mm box,
air, the inlet raised, initial velocity 25 mm/s, 16 s, automatic time step), on the GPU (Apple
M4 Max); the CPU gives the same outcome where run (75 and 90 mm/s at 2.5 mm: the two-relaxation-
time scheme unstable at steps 951 and 786, the regularised running 16 s). Re on the 50 mm
upright, ν = 15 mm²/s.

| inlet (mm/s) | Re | spacing (mm) | τ⁺ | two relaxation times | regularised |
|---|---|---|---|---|---|
| 25 | 83 | 2.5 | 0.536 | runs | runs |
| 40 | 133 | 2.5 | 0.5225 | runs | runs |
| 50 | 167 | 2.5 | 0.518 | unstable, step 2520 | runs |
| 60 | 200 | 2.5 | 0.515 | unstable, step 1332 | runs |
| 75 | 250 | 2.5 | 0.512 | unstable, step 911 | runs |
| 90 | 300 | 2.5 | 0.510 | unstable, step 742 | runs |
| 100 | 333 | 2.5 | 0.509 | unstable, step 684 | runs |
| 110 | 367 | 2.5 | 0.508 | | runs |
| 120 | 400 | 2.5 | 0.5075 | unstable, step 612 | unstable, step 2032 |
| 150 | 500 | 2.5 | 0.506 | unstable, step 438 | unstable, step 684 |
| 90 | 300 | 2.0 | 0.5125 | unstable, step 968 | runs |
| 90 | 300 | 1.25 | 0.520 | unstable, step 4672 | runs |
| 180 | 600 | 1.25 | 0.510 | unstable, step 750 | runs, sheds |

On the gallery grid (2.5 mm) the highest Reynolds number that runs goes from 133 to 367 and the
lowest τ⁺ from 0.5225 to 0.508. At 1.25 mm the wake sheds at Re 600: the cross-stream velocity
50 mm behind the bracket swings by ±150 mm/s at 1.07 Hz, a Strouhal number of 0.30 on the
upright's height. At 2.5 mm, 100–110 mm/s, it wavers by ±9–13 mm/s without settling into one
frequency. Beyond τ⁺ ≈ 0.508 the inlet and the outlet still reflect the start-up's pressure
waves; absorbing layers there are the next step.

**The heat sink of the electronics-cooling examples** (#587: the six-fin sink held at 80 °C in a
duct, 20 °C air, 1 mm cells, 4 s, automatic time step; GPU, the regularised 150 mm/s run also on
the CPU with the same heat):

| inlet (mm/s) | τ⁺ | two relaxation times | regularised |
|---|---|---|---|
| 100 | 0.5225 | 4.56 W, 13.2 K/W | 4.45 W, 13.5 K/W |
| 150 | 0.515 | unstable, step 704 | 6.08 W, 9.9 K/W |
| 200 | 0.511 | | 6.89 W, 8.7 K/W |
| 250 | 0.509 | | 6.89 W, 8.7 K/W |
| 300 | 0.5075 | unstable, step 303 | unstable, step 1023 |

**Small cases** in the test suite: a 48² lid-driven cavity at Re 1500 (τ⁺ 0.5096): two relaxation
times unstable at step 767, regularised runs (`a_lid_driven_cavity_at_tau_0_51_…`); a cylinder of
ten cells in the Schäfer–Turek channel at Re 200 (τ⁺ 0.510): unstable at step 1324, regularised
runs 6000 steps (heavy).

## Cost

| | two relaxation times | regularised |
|---|---|---|
| CPU, one thread, flow-past sample (million updates per CPU second) | 13.8 | 9.2 (−33 %) |
| GPU M4 Max, in place, periodic 128³ (million updates per second) | 1160–1250 | 920–1030 (−20 %) |
| GPU M4 Max, in place, sphere 128³ | 945–955 | 810–820 (−14 %) |

The CPU figure is CPU time (the build machine was shared, so wall-clock throughput was not
repeatable); the GPU figures are the range of three runs of
`cargo test --release -p fieldcraft-lbm-gpu --test throughput -- --ignored --nocapture`.

## References

- J. Latt, B. Chopard, *Lattice Boltzmann method with regularized pre-collision distribution
  functions*, Math. Comput. Simul. 72 (2006) 165–168, doi:10.1016/j.matcom.2006.05.017.
- C. Coreixas, G. Wissocq, G. Puigt, J.-F. Boussuge, P. Sagaut, *Recursive regularization step for
  high-order lattice Boltzmann methods*, Phys. Rev. E 96 (2017) 033306,
  doi:10.1103/PhysRevE.96.033306.
- C. Coreixas, B. Chopard, J. Latt, *Comprehensive comparison of collision models in the lattice
  Boltzmann framework: Theoretical investigations*, Phys. Rev. E 100 (2019) 033305,
  doi:10.1103/PhysRevE.100.033305.
- Z. Guo, C. Zheng, B. Shi, *Discrete lattice effects on the forcing term in the lattice Boltzmann
  method*, Phys. Rev. E 65 (2002) 046308, doi:10.1103/PhysRevE.65.046308.
- M. Geier, M. Schönherr, A. Pasquali, M. Krafczyk, *The cumulant lattice Boltzmann equation in
  three dimensions: Theory and validation*, Comput. Math. Appl. 70 (2015) 507–547,
  doi:10.1016/j.camwa.2015.05.001.
- T. Krüger et al., *The Lattice Boltzmann Method: Principles and Practice*, Springer 2017,
  doi:10.1007/978-3-319-44649-3.
- Q. Zou, X. He, *On pressure and velocity boundary conditions for the lattice Boltzmann BGK
  model*, Phys. Fluids 9 (1997) 1591–1598, doi:10.1063/1.869307.
- M. Schäfer, S. Turek, *Benchmark computations of laminar flow around a cylinder*, NNFM 52
  (1996) 547–566, doi:10.1007/978-3-322-89849-4_39.
- U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow using the
  Navier–Stokes equations and a multigrid method*, J. Comput. Phys. 48 (1982) 387–411,
  doi:10.1016/0021-9991(82)90058-4.
- D. d'Humières, I. Ginzburg, *Viscosity independent numerical errors for Lattice Boltzmann
  models: From recurrence equations to "magic" collision numbers*, Comput. Math. Appl. 58 (2009)
  823–840, doi:10.1016/j.camwa.2009.02.008.
- I. Ginzburg, *Equilibrium-type and link-type lattice Boltzmann models for generic advection and
  anisotropic-dispersion equation*, Adv. Water Resour. 28 (2005) 1171–1195,
  doi:10.1016/j.advwatres.2005.03.004.
- G. de Vahl Davis, *Natural convection of air in a square cavity: a bench mark numerical
  solution*, Int. J. Numer. Methods Fluids 3 (1983) 249–264, doi:10.1002/fld.1650030305.

All confirmed on Crossref. No code of any lattice Boltzmann package was read.
