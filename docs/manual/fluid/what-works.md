# What works: verified fluid capabilities

The fluid capabilities whose results have been checked against theory, a published benchmark or an
experiment, with the numbers, and the ones that have not. Each number is quoted from the page or
test that records it; nothing here was rerun for this page. Use it to pick examples and to word
posts: a capability is described as verified only for the case and resolution it was verified on.

Large-eddy simulation of bluff bodies and wall-modelled LES are **not** on this list. Their state,
the measured errors and what blocks them are in [RFC 0003](../../rfc/0003-lbm-les-roadmap.md).

## Laminar flow (CPU, f64)

| Capability | Case | Reference | Measured | Tolerance | Source |
|---|---|---|---|---|---|
| D3Q19 TRT collision | plane Poiseuille, body force, 4/8/16 cells, τ = 1 | parabola (theory) | exact to round-off | 10⁻¹⁰ of the peak | [Lattice Boltzmann](lattice-boltzmann.md), Verification |
| | Taylor–Green vortex, 16/32/64 cells | e^(−2νt) | 1.6e-2, 4.0e-3, 9.9e-4: order 2.0; mass to 1e-13 | finest < 2e-3, order > 1.8 | same |
| | Couette, 16 cells | linear profile, wall shear μU/H | profile 1e-10, shear 1e-6 | 10⁻¹⁰; shear 10⁻⁶ relative | same |
| Regularised and hybrid (HRR) collisions | Poiseuille 8/16/32 cells | parabola | order 2.00 (both) | |order − 2| < 0.2 (HRR) | [Collision operators](collision.md) |
| | shear wave, HRR, σ 0.5–0.985 | exact rate of the scheme's linear map | within 2.5e-7 | 1 % | same |
| MRT collision | Taylor–Green 16/32/64 | e^(−2νt) | orders 2.015, 2.004 | finest < 2e-3, order > 1.8 | same |
| Lid-driven cavity | Re 100, 64² cells | Ghia, Ghia & Shin 1982: u/U −0.2109 at y 0.4531 | −0.2136 at 0.4585 (TRT), −0.2139 (regularised) | — | same |
| Confined cylinder, steady | Schäfer & Turek 2D-1, Re 20, D/Δx 20, staircase walls | c_D 5.57–5.59 | 5.773 (+3.5 %) | — | same |
| Confined cylinder, shedding | Schäfer & Turek 2D-2, Re 100, D/Δx 20, staircase walls | c_D,max 3.22–3.24, c_L,max 0.99–1.01, St 0.295–0.305 | St 0.2987 (inside the band); c_D,max 3.447 (+7 %), c_L,max 1.102 (+9 %) | c_D 8 %, c_L 12 %, St 2 % (the outlet-distance test) | same, and [Boundary conditions](boundary-conditions.md) (outlet at 7.5 D: 3.356, 1.072, 0.2989) |
| Interpolated (Bouzidi) walls | round pipe, 8 and 16 cells per radius (oracle `pipe-poiseuille-interpolated-walls`) | a R²/(4ν) | −0.15 %, −0.04 % (halfway walls: −5.5 %, −2.1 %) | 0.5 %, 0.15 % | [Oracle results](../../oracle.md) |
| | annular Couette, sliding rod, 6–16 cells per radius | ln profile; rod drag 2πμUL/ln(R₂/R₁) | order 1.98; rod drag within 0.11 % (halfway: +1.5 to +14 %) | order > 1.6; drag 0.5 % | [Boundary conditions](boundary-conditions.md) |
| | periodic Stokes sphere array, R = 3, 6, 9 cells | Hasimoto 1959, Sangani & Acrivos 1982 | −1.64 %, −0.85 %, −0.32 % (halfway: +11.7 %, +1.7 %, +2.3 %) | 2 %, closer at every refinement | same |
| Inlets and outlets | parabolic inlet, pressure outlet, 16 cells | Poiseuille profile, gradient −12μŪ/H² | profile 3.0e-4, gradient +0.07 % | 5e-3, 1 % | [Lattice Boltzmann](lattice-boltzmann.md) |
| | characteristic outlet, plane pulse in a 2048-cell duct (oracle `characteristic-outlet-reflection`) | Selle, Nicoud & Poinsot 2004 one-pole response | departure 2.28e-3 of the peak; tail rate 1.0005 of K/2 | departure 4e-3; rate 0.5 % | [Boundary conditions](boundary-conditions.md) |
| Inlet ramp (`inlet_ramp`) | empty 192-cell tube | one-dimensional solution by characteristics | within 0.13 % (m = 1), 0.08 % (m = 2) | — | [Boundary conditions](boundary-conditions.md), "Raising the inlet" |
| | 2-D cylinder Re 3900, 16 cells per D, far-field C_p rms | ideal tube × the run's measured decay | no ramp 20.19; ramps of 1, 2, 2.5 T₀: 7.020, 1.350, 0.049 (ideal tube × decay 7.008, 1.357, 0.058) | printed, not asserted | same |
| Shear-dependent viscosity | power-law Poiseuille, n = 0.5 and 1.5 (oracle) | Bird, Armstrong & Hassager | within 2.5e-5 to 3.8e-4 m/s of 0.04–0.30 m/s | 3.3e-4, 6.0e-4 m/s | [Oracle results](../../oracle.md) |
| Free surface | hydrostatic column; plug filling (oracle) | hydrostatics; volume conservation | 0.3 %; 0.11 %; mass to 1e-10 | 1.5 %; 0.5 %; 1e-10 | same |
| | dam break (oracle `dam-break-martin-moyce`, experiment) | Martin & Moyce 1952 | front times +10 to +13 %; column height +0.5 % | 14–16 %; 5 % | same |
| Determinism (threads) | a channel with an inlet, outlet, slip top, body force and a block, plain TRT with halfway, linear and quadratic interpolated walls; heat with buoyancy | one thread against 2, 3 and 7 | bit for bit | exact | `crates/lbm/tests/threads.rs` |
| Determinism (repeat) | one ball in a channel, solved twice | the first solve | bit for bit | exact | `results_are_bit_for_bit_reproducible` (`crates/lbm/tests/verification.rs`) |

## Heat (CPU, f64)

| Capability | Case | Reference | Measured | Tolerance | Source |
|---|---|---|---|---|---|
| Natural convection | de Vahl Davis cavity, Ra 10³ to 10⁶, Pr 0.71 (oracle) | de Vahl Davis 1983 | mean Nu −0.012 %, +0.15 %, +0.35 %, +0.54 %; Nu_max at Ra 10⁶ −2.7 % (known failure) | mean Nu 0.13 %, 0.27 %, 0.73 %, 1.3 %; every other quantity has its own in docs/oracle.md | [Heat transfer](heat-transfer.md), [Oracle results](../../oracle.md) |
| Conjugate heat transfer | cavity with a centred conducting body, k_s/k_f 0.2, 1, 5 (oracle) | House, Beckermann & Smith 1990 | +0.21 %, +0.20 %, +0.18 % | 0.6 % | [Heat transfer](heat-transfer.md) |
| Conduction through solids and contacts | composite walls, contact resistance, thin walls | series conduction | exact discrete solution | 1e-9 | same |
| Radiation from a conducting slab | 16-cell slab (oracle `radiation-slab-lattice`) | conduction–radiation balance | < 1e-6 relative | 1e-8, 1e-6 relative | same |
| Melting and freezing | two-phase Stefan problem (oracle) | Neumann's solution | front +0.27 % on 128 cells | 0.4 % | same |
| | channel frozen on its cold quarter (oracle) | plane Poiseuille in the liquid | +1.47 % | 2.2 % | same |

## Turbulence: what is verified

| Capability | Case | Reference | Measured | Tolerance | Source |
|---|---|---|---|---|---|
| Smagorinsky's model | Couette with Smagorinsky's model, HRR σ 0.9 | ν_t = (CΔx)²U/h | wall shear 3e-6, ν_t/ν 8e-4 relative | 2e-3, 2e-2 | [Collision operators](collision.md) |
| Spalart–Allmaras (RANS) | wall-resolved Re_τ 180 channel | the model's own 1-D solution | second-order convergence; the model itself +1.1 % in U_b against Moser, Kim & Mansour | — | [Turbulence](turbulence.md), verification.md |
| SA constant-stress wall law | Re_τ 1000 channel, first node at y⁺ 31–100 | wall-resolved 1-D SA | C_f +0.8 to +7.3 % | — | verification.md, "Spalart–Allmaras: its wall law in the channel" |

**Smoke-tested, no accuracy claim.** The matched on-node wall (`spalding_bounce_matched`) at
Re_τ 1000, 8 cells per δ, 2 δ/u_τ discarded and 2 averaged (the heavy job's bounded run): U_b⁺
19.792 against Lee & Moser's 19.990 (#1121; the run asserts stability and U_b⁺ within 30 %). Its
full study on the GPU (16, 24, 32 cells per δ to the statistical protocol) is **inconclusive (time
cap): not passed**, with C_f −2.9 % at 16 cells, −6.7 % at 24 and −8.7 % at 32 (#1151). Do not post
either as an accuracy result.

## GPU backend (wgpu, f32)

| What | Measured | Source |
|---|---|---|
| Agreement with the CPU | stated tolerances: velocity 2e-5 of the peak, density 1e-6, wall forces 5e-5. Measured on the M4 Max: velocity 5e-7 to 7.4e-6, wall shear 1.9e-5; on the GB10, Schäfer–Turek C_d 7e-5 % and St 4e-4 % off the CPU (#785) | [GPU backend](gpu.md) |
| Throughput, TRT, in place, periodic 128³ (million cell updates per second) | RTX 4090 3897, M4 Max 1364, GB10 1224 | same, "Performance" |
| Throughput, TRT, sphere 128³ | RTX 4090 3572, M4 Max 1256, GB10 1146 | same |
| Throughput, HRR σ 0.99 / + Vreman / dynamic σ, periodic 128³, GB10 | 777 / 773 / 772 (regularised 1207) | same |
| Against the CPU on a whole model | `flow-past` at 1.25 mm, 1.47 M cells, 800 steps: CPU 18.1 s (65 M/s), GPU 3.6 s (328 M/s), M4 Max; face forces within 1e-7 of the CPU's | same, "Worked example" |
| A large model | AIO radiator air side, 5.07 M cells, 20 513 steps: about four minutes on the GB10 | [Examples](../../examples.md), "AIO cooler radiator" |

The CPU figure for one thread is 13.8 million updates per CPU second (TRT, `flow-past`;
[Collision operators](collision.md), "Cost"). The CPU stops scaling at about four threads (#643);
quote GPU speed-ups against a stated CPU configuration, never as one ratio.

## OpenFOAM export and read-back

`deck.openfoam` writes a fluid step as an ESI OpenFOAM v2406 case, incompressible (pimpleFoam) or
compressible (rhoCentralFoam, a laminar ideal gas); `results.openfoam` reads the results back for
viewing ([OpenFOAM cases](openfoam.md)). What is tested: every written dictionary against our own
reader of the documented grammar, the round trip of what we write, and read-backs of files we build
in the documented ASCII format (reconstructed and decomposed cases, refusals of hostile counts, of
uncovered cells, of non-positive T and p). What is not: no case written by FieldCraft has been run
in OpenFOAM by the project. One compressible case still needs a run on v2406 by someone who has it
(#1094). An OpenFOAM run is an interoperability check; it never verifies FieldCraft's results.

## Proposed showcase examples

Run times marked *measured* are from the page cited; the others are estimates from cells × steps
over the throughput above and must be measured before they are posted.

| # | Model | What it shows | Run time | Where |
|---|---|---|---|---|
| 1 | Schäfer & Turek 2D-2: cylinder in a channel at Re 100, D/Δx 20, staircase (halfway) walls, convective outlet; animate the vorticity | periodic shedding with St 0.2987 inside the benchmark's 0.295–0.305 (c_D,max +7 % on this staircase grid); the CPU and the GPU agreeing on C_d to 7e-5 % (#785) | a few minutes on the CPU (estimate) | CPU, then GB10 for the side-by-side |
| 2 | de Vahl Davis cavity at Ra 10⁵, then the House et al. cavity with a conducting body | mean Nusselt number within 0.35 % (Ra 10⁵) and 0.18–0.21 % (House et al.) of the benchmarks; temperature and streamlines | about a minute each (estimate; oracle cases) | CPU |
| 3 | AIO cooler radiator, air side (`cooling-aio-half`) | 5 M cells of real geometry from SolidCraft in about four minutes; 61 W at effectiveness 0.93, grid change 0.4 % | about 4 min (measured) | GB10 |
| 4 | Inlet ramp, two runs of the 2-D cylinder in a 24 D × 12 D box: started at full speed and ramped over 2.5 T₀ | the start-up pressure wave and its removal: far-field C_p rms 20.19 without the ramp, 0.049 with it, against 0.058 from the characteristics solution times the run's decay | about 10–15 min per run on the CPU (estimate) | CPU |
| 5 | Staircase against interpolated walls: the round pipe at 8 cells per radius, or the sphere array at R = 3 | −5.5 % against −0.15 % on the pipe; +11.7 % against −1.6 % on the spheres | seconds to a minute (estimate) | CPU |
| 6 | Lid-driven cavity at Re 100, 64² (the pending oracle case's recipe) | the primary vortex and the centreline profile on Ghia et al.'s points | under a minute (estimate) | CPU |
| 7 | `cooling-heat-sink` | mixed convection with buoyancy, 1.93 W at 13.0 K/W | about 4 min (measured) | CPU |
| 8 | GPU throughput on the GB10, M4 Max and RTX 4090 | 1.2–3.9 billion cell updates per second, CPU/GPU agreement on the same model | the throughput test, minutes (`--test throughput`) | GB10, M4 Max |
| 9 | Dam break against Martin & Moyce | a free surface against measured front positions, mass exact | minutes (estimate) | CPU |
| 10 | A model exported with `deck.openfoam`, run in OpenFOAM by its user, read back with `results.openfoam` and shown with FieldCraft's slices | interoperability only, never a comparison of answers | export seconds; the OpenFOAM run is the user's | CPU; needs OpenFOAM |

Not for posts as accuracy: bluff-body LES (cylinder Re 3900 C_d 1.07–1.83 against 0.99 depending on
the setup and the start), wall-modelled LES channels, the Re_τ 5200 channel (not runnable below
about 42 cells per δ). See [RFC 0003](../../rfc/0003-lbm-les-roadmap.md).
