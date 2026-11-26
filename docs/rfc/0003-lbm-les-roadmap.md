# RFC 0003: Lattice Boltzmann LES, where it stands and what comes next

- **Status:** proposed; the LBM-LES work is parked as of 2026-11-26. Nothing in this RFC is decided.
- **Discussion:** [#785](https://github.com/tangency-labs/fieldcraft/issues/785) (tracking), with its sub-issues
- **Scope:** large-eddy simulation on the lattice Boltzmann solver (`fieldcraft-lbm`, `fieldcraft-lbm-gpu`):
  bluff bodies, wall-modelled channels, the collision near τ = ½, grid refinement and absorbing layers.
  Laminar flow, heat, the GPU backend and the OpenFOAM export are outside it; what of those is verified is
  in [What works](../manual/fluid/what-works.md).

## 1. Summary

FieldCraft's LES has every ingredient the published LBM-LES successes use except the resolution: the
hybrid recursive regularised collision with Vreman's model and dynamic σ (Jacob, Malaspinas & Sagaut 2018),
interpolated walls, a two-step force, an inlet ramp, a matched wall model. On the uniform grids that run in
minutes to an hour on the GB10 (12–32 cells per diameter, 16–32 per channel half-height) no bluff body and no wall-modelled channel
meets a published reference within 5 %. The four blockers, in the order they cost accuracy: resolution
(no refinement yet), the start-up and its instability, the collision's linear instability near τ = ½, and
wall models whose matching height the grids cannot place where Kawai & Larsson put it. Refinement is half
designed and partly built (#1107). This RFC records the numbers, the options with their costs, a
recommended order and the decisions the owner must make before the work resumes.

## 2. Where it stands (numbers)

All bluff-body numbers are single realisations on unconverged grids, block intervals of one run; none is
an acceptance run. References: circular cylinder at Re 3900, C_d 0.99 ± 0.05 (experiment), 1.04 (LES),
St 0.215 (Kravchenko & Moin 2000); square cylinder at Re 22 000, C_d ≈ 2.1, St 0.132 (Lyn et al. 1995),
2.18 (Trias et al. 2015 DNS); channel, Lee & Moser 2015.

### 2.1 Bluff bodies

| Run | Setup | C_d | Against the reference | Record |
|---|---|---|---|---|
| #787 original | cylinder, 16 cells per D, free stream everywhere, staircase, one-step force | 1.354 ± 0.040 | +37 % | #787 |
| expert probe | same, at rest on the body, halfway | 1.260 ± 0.012 | +27 % | verification.md, "Bluff-body LES" |
| expert probe | same, interpolated walls | **1.113 ± 0.012** | +12 % | same |
| expert probe | same, 24 cells per D | 1.067 ± 0.019 | +8 % | same |
| #1106 rerun (PR #1147) | HRR dynamic σ + Vreman, interpolated, at rest + ramp 2.5 T₀, 16 cells | **1.827 ± 0.073** | +85 %; St 0.057 against 0.215 | PR #1147 |
| #1106 rerun | HRR σ 0.99 + Vreman, same | unstable at 3.6 D/U (lattice speed 0.547) | — | PR #1147 |
| #1109 rerun (PR #1149) | HRR σ 0.99 + Smagorinsky, 12 cells, uniform start, two-step force | 1.799 ± 0.018 | +82 % | PR #1149 |
| #1109 rerun | same, dynamic σ | 1.620 ± 0.043 | +64 % | PR #1149 |
| #1109 rerun | same, Musker wall model | −0.026 ± 0.064, not stationary | — | PR #1149 |
| square, expert probe | 16 cells, at rest, two-step force | 1.81 ± 0.07 | −14 % | verification.md |
| square, #1106 rerun | HRR dynamic σ + Vreman, ramp | 2.176 ± 0.088 | +4 % (Lyn); St 0.057 against 0.132 | PR #1147 |
| square, #1106 rerun | σ 0.99 + Vreman, ramp | unstable at 3.2 D/U | — | PR #1147 |
| #1110 noise study (PR #1148) | cylinder, HRR + Smagorinsky, halfway, ramp; 16 and 32 cells, U_lat 0.05 and 0.10 | all four diverge, 12–73 D/U | no noise statistic | PR #1148 |

Reproduction on develop (2026-11-26, GB10, #1106 comment; Smagorinsky C 0.1, σ 0.99, no ramp, 16 cells per
D, span 3.125 D, 100–200 D/U):

| start | C_d ± 95 % | St | L_r/D |
|---|---|---|---|
| at rest on the body, no cross-flow wave | 1.437 ± 0.009 | 0.250 | 0.60 |
| at rest on the body, a 1 % cross-flow wave (the probe's start) | **1.113 ± 0.012** | 0.220 | 2.1 |

| case | ramp | result |
|---|---|---|
| cylinder, Vreman 0.18, σ 0.99 | 2.5 T₀ | unstable at 3.6 D/U (lattice speed 0.547) |
| cylinder, Smagorinsky 0.1, σ 0.99 | 2.5 T₀ | unstable at 137 D/U (lattice speed 0.685) |
| square, Vreman, σ 0.99 | 2.5 T₀ | unstable at 3.2 D/U (lattice speed 0.538) |
| cylinder, Vreman, σ 0.99 | none | stable to 240 D/U, C_d 1.422, St 0.243 |

What this settles and what it leaves:
- **The recorded 1.113 is start-sensitive, not a regression.** It reproduces to the digit with the probe's
  1 % cross-flow wave; without the wave the same setup gives 1.437 (+29 %). INFERRED, not checked: the seed decides
  whether the wake goes three-dimensional on a 16-cell grid. The #1109 baseline (1.80, uniform start) equals the old
  one-cell-span value, consistent with a wake that stays two-dimensional (INFERRED, not checked).
- **The dynamic-σ rows are not results yet.** Their 0.057 lift peak has odd harmonics at 0.171 and 0.293 and
  a lift range of ±5; it is not a shedding frequency, and the cylinder's 1.827 is not explained.
- **The ramp destabilises fixed σ 0.99 while the inlet is still near zero.** Without the ramp the same run is
  stable to 240 D/U. Cause not found; reproducer on the GB10 (about 13 s to the failure):
  `FIELDCRAFT_BLUFF_CASE='r:shape=cylinder,sgs=vreman,cs=0.18,sigma=0.99,ramp=2.5,start=5,end=40'` with the
  ignored test `bluff_body_les_rerun` (PR #1147).
- **The #1110 divergences are not explained.** Amplifying a 10⁻² departure by 100 takes ln 100 / rate steps:
  about 1.5 × 10⁵ at the slowest linear rate where these cases run (3e-5 per step), about 2.3 × 10⁴ at the
  fastest (2e-4). The #1110 runs diverged in 3 800 to 25 700 steps, one of them during its ramp, with lattice
  speeds above 0.5 near the body. So the linear mode cannot be ruled out at the fastest rate; a nonlinear,
  start-up-driven blow-up is the other candidate (INFERRED, not tested).

### 2.2 Wall-modelled channels

| Wall treatment | Grid | C_f against Lee & Moser | Record |
|---|---|---|---|
| none (halfway bounce-back), Re_τ 1000 | 16 cells per δ | +173 to +261 % | verification.md |
| Malaspinas & Sagaut (`musker`) | 16, 32 | −51 to −74 %, still drifting | same |
| Han, Ooka & Kikumoto (`spalding_bounce`) | 16, 32 | −46, −49 % | same |
| Han's wall, law sampled at rows 3–6 (probe) | 16, 24, 32 | +5.8, +5.4, +5.1 % | same |
| `spalding_bounce_matched`, bounded run, 8 cells, 2 + 2 δ/u_τ | 8 | U_b⁺ 19.792 against 19.990 (−1.0 %) | heavy job, #1121 |
| `spalding_bounce_matched`, GPU, full protocol, Re_τ 1000 | 16, 24, 32 | −2.85, −6.73, −8.71 %; time cap, inconclusive | PR #1151 |
| same, Re_τ 2000 | 16, 24, 32 | −1.83, −5.88, −8.47 %; time cap | PR #1151 |
| same, Re_τ 5200 | 16, 24, 32 | refused: τ 0.50004–0.50008 below the bound 0.5001 | PR #1151 |
| Han's wall at A20, hybrid collision | — | −52 to −55 % (MRT row void) | verification.md, #1057 |

The matched wall's C_f error grows with resolution at both Reynolds numbers, and the worst row above the
sample height has U⁺ +5.6 to +7.1 % above Lee & Moser (+5.58 % at Re_τ 2000, 16 cells, to +7.11 % at
Re_τ 1000, 32 cells; PR #1151). Every grid has h_wm/δ above Kawai & Larsson's 0.1 (4 cells are 0.25δ, 0.167δ, 0.125δ). The cause is
not isolated.

### 2.3 What is verified and holds

The force budget closes per step to 2–3 × 10⁻¹⁴; CPU and GPU agree to 7 × 10⁻⁵ % on the laminar control;
the Smagorinsky term has no factor-of-two error (#785 "Checked and found consistent"). The HRR, Vreman and
dynamic σ implementations were read against Jacob et al. 2018 equation by equation by the #785 expert, with
the differences already stated in `crates/lbm/src/hrr.rs`'s docs (unverified here: that reading is in the
maintainer's notes, not in the repository). The inlet ramp
matches the one-dimensional characteristics solution within 0.13 % (#1134). These are not in question.

## 3. What blocks published-accuracy bluff bodies

### 3.1 Resolution and refinement
The published LBM-LES of the cylinder at Re 3900 that the #785 expert read use D/Δx ≈ 80 at the wall through
nested levels (Spinelli et al. 2023; Jacob et al. 2018, D/80 with 10 M points; as read by the expert, not
re-checked here). FieldCraft runs 12–24 on one uniform grid. At 16 cells the separating shear layers sit on
one or two cells.

A uniform grid of 80 cells per D in the 24 D × 12 D × 4 D box is 24 · 12 · 4 · 80³ = 5.9 × 10⁸ cells. On the
GB10 (computed, not run):
- **Memory:** 76 bytes of populations per cell in place plus 12 of velocity for the hybrid collision
  ([GPU backend](../manual/fluid/gpu.md)) is about 52 GB, inside the GB10's 128 GB of unified memory. But
  the population state lives in one storage buffer (gpu.md, Limitations), and its 4.5 × 10¹⁰ bytes are ten
  times the 4 294 967 292-byte limit per buffer the backend reports (gpu.md, worked example). It cannot run
  until #699 splits the buffer.
- **Time:** at lattice speed 0.05 one D/U is 80 / 0.05 = 1 600 steps; 150 D/U (ramp, transient, averaging)
  is 2.4 × 10⁵ steps, 1.4 × 10¹⁴ cell updates. At the bluff-body throughput measured on the GB10 (0.56 G
  updates/s, HRR + Smagorinsky, #1107 design §4) that is 2.5 × 10⁵ s, about 70 hours; at the 773 M/s of the
  periodic HRR + Vreman benchmark (gpu.md) about 50 hours.
- **Refined** (finest box 4 D × 3 D × 4 D at D/80, about 25 M cells, the design's estimate): about 5 hours.

So uniform D/80 becomes possible once #699 lands, at 50–70 GB10-hours a run; refinement brings it to about 5,
which is what makes a resolution study (D/40, D/80) affordable.

State of #1107 (its hand-off comment of 2026-11-26 carries the design and the PR plan):
- Decided by the owner (2026-11-25/26): scheme C (cell-centred, reconstructed; Qi, Klimach & Roller 2019),
  regularised ghosts for RR/HRR and full rescaled ghosts for TRT, one nested box per level, coupling stencil B
  (eight sources), absorbing layers of Xu & Sagaut 2013 type II with a required reflection R.
- Built: PR #1150, the absorbing layer's model and command (commit 1 of 2); the collision term is the open
  second commit. Levels in the model and setup: built, uncommitted, branch `lbm-1107-refine` (not pushed).
- Not built: two-level coupling (briefed, ungated), many levels with HRR and LES across
  interfaces, refined runs in `solve.fluid`, the GPU, the D/Δx 40 and 80 cylinder.
- Measured in the design's prototypes (numpy, not in the repository; #1107): the Couette check confirms the derived rescaling τ_f/(2τ_c); Qi's printed form gives
  a 15 % stress jump. Regularised ghosts put TRT Poiseuille 5.6e-3 off, so TRT keeps full ghosts.

### 3.2 The start
A start with the free stream inside the body sends an acoustic transient that never leaves the box (force
rms 3.7 against 0.41 from rest, +4 % on the mean; the square went unstable at 161 D/U). The ramp (#1118)
removes the inlet's pulse (C_p rms 20.19 → 0.049 at 2.5 T₀). With the ramp, the stream-function start at
rest on the body has lattice speeds above the inlet's near the body, and every fixed-σ run and all four
#1110 runs diverge.

### 3.3 The collision near τ = ½
The regularised and hybrid collisions are linearly unstable in a stream below τ⁺ ≈ 0.52 (von Neumann maps
over 14 directions, k on a π/32 grid; numbers in #1154). At the default lattice speed 0.05 each collision has
its own bound, and with τ⁺ − ½ = 3 U N / Re at Re 3900:
- **regularised (σ 1):** τ⁺ ≥ 0.502, which needs N ≈ 52 cells per D;
- **hybrid (σ 0.99):** τ⁺ ≥ 0.5005, met from N ≈ 13 cells per D.

The hybrid cylinder on 16 cells runs at τ⁺ 0.5006, just inside its bound; a regularised run there is not.
`MIN_TAU_TURBULENT` 0.5001 admits both. A refinement hierarchy makes it worse on the coarse levels: at D/80
and speed 0.05 the finest τ⁺ − ½ is 3 · 0.05 · 80 / 3900 = 3.1 × 10⁻³, and with five levels the coarsest is
2⁻⁴ of it, 1.9 × 10⁻⁴, below the hybrid's bound at that speed (INFERRED from the maps; #1107 design §4). Density
fields stay at ρ ±7–17 % in the runs that complete.

### 3.4 Wall models
Bluff bodies at Re 3900 need no wall model at D/80 (y⁺ ≤ 1 at the wall, published). Wall models matter for
higher Re and for channels and plates. The first-node forms of the papers carry half the drag or less (C_f −46 to −74 %, §2.2); sampling three or
more rows up gives +5 %; the matched wall at 4 cells is inconclusive at the time cap (not passed), with an
error that grows with the grid (§2.2). Wall models on interpolated walls (#927) and on obstacles (Han stack #1046–#1048) wait on
this.

## 4. Options and costs

Effort is agent-hours (builder and expert work, reviews included) plus GB10 hours; the GB10 is borrowed.

| # | Option | What it buys | Effort | Cost if it is the wrong call |
|---|---|---|---|---|
| A | The start and the ramp: find why the 2.5 T₀ ramp destabilises fixed σ 0.99 with the inlet near zero (13 s reproducer); fix a seeded start the product documents (the 1 % cross-flow wave moves C_d from 1.437 to 1.113); explain the dynamic-σ rows' 0.057 peak | runs that complete with the ramp; a C_d that does not hang on an undocumented seed; the baseline every later step needs | M: 8–12 agent-h, 3–5 GB10-h | none: every later bluff-body number depends on it |
| B | #1110's 2 × 2 noise split rerun on the start A settles | the regularised instability's share of the ρ ±7–17 % noise against under-resolution | S: 3–4 agent-h, 2 GB10-h | none beyond its GB10 time |
| C | #785 step 6: warn (recommended by the expert) or refuse below the stream-stability bound | users are told when the run is inside the unstable region | S: 2–3 agent-h, no GPU | warn: noisy runs keep completing, flagged; refuse: at Re 3900 the regularised collision needs about 52 cells per D (bound 0.502), the hybrid σ 0.99 about 13 (bound 0.5005) |
| D | Finish #1107 on the CPU: layer term (PR #1150 commit 2), levels setup, two-level coupling, many levels with HRR and LES, refined `solve.fluid` | the route to 80 cells per D | L: 40–60 agent-h, no GPU until the end | a scheme combination (Qi + Schukmann + Xu–Sagaut + Bouzidi) that is ours in composition fails a check; each PR's bars catch it before the next |
| E | #1107 on the GPU (two buffers per level; AA is not compatible, Schornbaum & Rüde §4.1) | the D/Δx 80 run in hours rather than days | L: 20–30 agent-h, 4–6 GB10-h of parity runs | a second implementation to keep in step |
| F | The D/Δx 40 and 80 cylinder (#1107 PR 8) | the acceptance run against experiment (C_d 0.99 ± 0.05, St 0.215) | S agent work; GB10 ≈ 20 min at D/40, ≈ 5 h at D/80 refined (#1107 design §4, INFERRED from 0.56 G updates/s); uniform D/80 needs #699 and 50–70 h | the bar must be signed before the run, from experiment: C_d 0.99 ± 0.05, St 0.215 (Kravchenko & Moin 2000 and the experiments they cite) |
| G | Matched wall: isolate the C_f trend (CPU run of one GPU grid; h_wm fixed in outer units at ≥ 40 cells per δ; Re_τ 5200 at ≥ 42 cells per δ) | knows whether the trend is the model, the height or f32 | M: 6–8 agent-h, 6–10 GB10-h (32 cells took 3 464 s; 42 cells ≈ 2.9 h per Re_τ, estimated as N⁴ from that) | keeps a wall model whose error grows with the grid |
| H | A collision stable near τ = ½ (cumulant on D3Q27, Geier et al. 2015) | removes 3.3 at its root; the published alternative to HRR | L: 60+ agent-h, new lattice | the largest change; D3Q27 costs 42 % more memory and traffic |
| I | Stop LES at what works and document the limits | the team moves on; the manual states the measured errors (done by #1104/#1116) | none | users who need bluff-body drag at high Re go elsewhere |

## 5. Recommended order

1. **A** (the start and the ramp) first: no later bluff-body number means anything while C_d moves 29 % with
   an undocumented seed and the ramp kills fixed-σ runs. It uses the drivers already written (PR #1147, #1149).
2. **C**, the warning, at the same time: two hours, and it stops silent noisy runs.
3. **B**, #1110's 2 × 2 rerun on that start.
4. **D**, refinement on the CPU, in the design's PR order: layer term, levels, two-level coupling,
   many levels with HRR, refined `solve.fluid`. Each PR's bars are fixed in the design before it runs.
5. **F** at D/Δx 40 on the CPU or a GPU port of the uniform-grid parts, then **E** and **F** at D/Δx 80.
6. **G** in parallel with D (different crate areas, GB10 time only between D's milestones).
7. **H** only if F at D/Δx 40 shows the coarse levels going unstable (design §4).

Merging the measurement PRs (#1147, #1148, #1149, #1151) before any of this keeps the numbers on develop.
They change no solver code.

## 6. Decisions the owner must make

1. **#785 step 6:** warn (expert's recommendation) or refuse below the stream-stability bound, and at which
   speed (the characteristic one, known at setup).
2. **The start for bluff-body LES:** whether the product seeds a cross-flow wave (and how it is stated: the
   1 % seed is the probe's choice, not a cited value), once A has found why the ramp destabilises fixed σ.
3. **The D/Δx 80 run's bar and GB10 time:** the tolerance, set from experiment (C_d 0.99 ± 0.05, St 0.215:
   Kravchenko & Moin 2000 and the experiments they cite; never other solvers' spread), and the ≈ 5 h of a
   borrowed box.
4. **The matched wall's fate if G confirms the trend:** keep it with a warning, change the sample height
   rule, or withdraw it.
5. **#1117's laminar slip** (a(4Δx)²/(2ν) at the on-node wall): accept or veto.
6. **Contributor SA wall-layer stack (#1130–#1132):** the method choices listed in #1016's hand-off comment
   (row-5 cutoff derived from #1130's reference; the harmonic rule derived from Chapman–Enskog or dropped;
   the thermal path documented).
7. **Whether to resume LES at all (option I)** before the work now ahead of it, such as M24 (polymer melt and
   creeping flows, ROADMAP.md), which was queued to start after LES.
8. **M4's four forks ([#1158](https://github.com/tangency-labs/fieldcraft/issues/1158)),** outside CFD but
   parked with it: F1 what "M4 done" covers; F4 the formulations of #1068, #1071 and #1073; F5 the input form
   of a temperature that varies in space (#1074); F6C a drilling membrane for the twisted beam (#1077).

## 7. Hand-off issues

Every unfinished item has an issue with its branch, numbers, next step and effort
(labels `area:cfd`, `paused`; `help wanted` where no GPU box or owner decision is needed):

| Issue | Item | Effort |
|---|---|---|
| #785 | tracking, with the step table | — |
| #1153 | the start seed (1.113 / 1.437) and the ramp with fixed σ | M, 3–5 GB10-h |
| #1154 | step 6, warn or refuse (owner) | S |
| #1106 | bluff-body reruns, closes with PR #1147 | S |
| #1109 | #788 rows after the two-step force, PR #1149 | S, 2 GB10-h |
| #1110 | noise split, PR #1148 | S, 2 GB10-h |
| #1108 | matched wall on the GPU, Re_τ 5200, the C_f trend, PR #1151 | M, 6–10 GB10-h |
| #1107 | refinement and absorbing layers, PR #1150 | L, 6–8 GB10-h |
| #1057 | the A20 GPU table after #1117 | S, 4–6 GB10-h |
| #1155 | one definition for the WGSL and native constants | S |
| #1156 | the GB10 segfault in `concurrent_solves_on_the_shared_device_read_back` | S–M |
| #1157 | `encode` checks fluid results only | S |
| #1094 | one compressible case on OpenFOAM v2406, by a person | S |
| #1016, PR #1128 | the contributor's SA wall-layer stack (#1130–#1132) and trace relaxation | M; S |
| #1158 | M4's owner forks F1, F4, F5, F6C | per element |

## 8. References

- Jacob, Malaspinas & Sagaut, J. Turbul. 19:1051, 2018, doi:10.1080/14685248.2018.1540879.
- Kawai & Larsson, Phys. Fluids 24:015105, 2012.
- Kravchenko & Moin, Phys. Fluids 12:403, 2000, doi:10.1063/1.870318.
- Lee & Moser, J. Fluid Mech. 774:395, 2015, doi:10.1017/jfm.2015.268.
- Lyn, Einav, Rodi & Park, J. Fluid Mech. 304:285, 1995, doi:10.1017/S0022112095004435.
- Malaspinas & Sagaut, J. Comput. Phys. 275:25, 2014, doi:10.1016/j.jcp.2014.06.020.
- Qi, Klimach & Roller, Comput. Math. Appl. 2019 (arXiv 1601.06984).
- Schornbaum & Rüde, SIAM J. Sci. Comput. 38:C96, 2016 (arXiv 1508.07982).
- Spinelli, Horstmann, Masilamani, Soni, Klimach, Stück & Roller, *HPC performance study of different collision
  models using the Lattice Boltzmann solver Musubi*, Comput. Fluids 255:105833, 2023,
  doi:10.1016/j.compfluid.2023.105833 (HRR and cumulant on the cylinder at Re 3900).
- Xu & Sagaut, J. Comput. Phys. 245:14, 2013 (arXiv 1203.6350).
- Geier, Schönherr, Pasquali & Krafczyk, Comput. Math. Appl. 70:507, 2015 (cumulant collision).
- Trias, Gorobets & Oliva, Comput. Fluids 123:87, 2015, doi:10.1016/j.compfluid.2015.09.013.
