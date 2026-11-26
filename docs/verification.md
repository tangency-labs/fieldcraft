# Verification and validation

Every solver capability is checked against three kinds of reference, and every accepted result is
frozen so it cannot drift.

| Kind | What it checks | Sources |
|---|---|---|
| Theory | The code solves the equations right | Closed-form solutions (beams, plates, Lamé, Kirsch, waves, shocks, Poiseuille…), patch tests, convergence rates, mesh invariants |
| Experiment | The model predicts reality | Measurements published in papers and reports (Taylor bar, cavity and cylinder flows, FSI experiments…) |
| Standard benchmarks | The code agrees with the published targets of standard problems | NAFEMS benchmarks, MacNeal–Harder and other published standard problems and reference solutions. Values openly published in benchmark reports and guides (ESRD, vendors' benchmark guides) may be cited; we never run the product that made them |

Rules:
- Each case is **our own model** of the problem, built from the published description (geometry,
  material, loads, boundary conditions). We never run a commercial product and never copy its
  input files, manual text or figures; only the published numbers are used, with a citation.
- Closed-source or non-permissively licensed software is never used to verify a result,
  brute-force a value or decide behaviour. Results a contributor shares from their own
  commercial licence are a cross-check in the discussion, never a reference value (AGENTS.md).
- Each case lives in `tests/oracle/<case>/`:
  - `recipe.json` — the model as FieldCraft commands (or a SolidCraft recipe for geometry);
  - `expected.json` — the reference values with `kind` (theory | experiment | benchmark),
    `source` (citation: document, table or figure, page), `quantity`, `value`, `tolerance` and
    the reason for that tolerance (mesh, measurement scatter, model differences);
  - `fixture.json` — FieldCraft's accepted result, written by `cargo xtask oracle --bless <case>`
    only after the case passes against `expected.json`.
- `cargo xtask oracle` runs every case, compares with `expected.json` (pass/fail with error) and
  with `fixture.json` (any change beyond the fixture tolerance fails CI, even if the reference
  still passes), and writes `docs/oracle.md` with the table of cases, sources, errors and status.
- A fixture is updated only deliberately (`--bless`), in a PR that says why the result moved
  (a better element, a finer default mesh), never to hide a regression.
- The theory suite is the gate for every merge; experiment and benchmark cases run in CI too, with
  their cited tolerances.

## The case format

`cargo xtask oracle` builds `fieldcraft-testkit`'s `oracle` binary (release) once and runs every
case through the command engine; the results are in [oracle.md](oracle.md).

```
cargo xtask oracle                      # every case; writes docs/oracle.md
cargo xtask oracle --case nafems-le1    # one case (repeatable); docs untouched
cargo xtask oracle --bless nafems-le1   # freeze a passing case's results (or --bless all)
cargo xtask oracle --quick              # skip cases marked "slow" (and "nightly")
cargo xtask oracle --nightly            # also run the cases marked "nightly"
```

- `recipe.json` is a script `fieldcraft-cli run` also accepts: a JSON array of
  `{"command", "params"}`; a step may carry an `"id"` so a quantity can read its output. Meshes
  come from `mesh.block` (mapped hex, tet or wedge blocks with arcs and elliptical edges, named
  faces, edges and vertices that sets select), `mesh.tet.generate`, or nodes and elements written
  one by one; once structured and automatic meshers land, cases may use them.
- `design.json` (optional) is a SolidCraft script that builds the geometry, opened before the
  recipe (the test corpus's recipe form).
- `expected.json`:
  ```json
  {"title": "…", "kind": "theory|benchmark|experiment",
   "source": "document, section/table/figure, page",
   "status": "active|pending", "requires": "the capability a pending case waits for",
   "slow": false, "nightly": false,
   "quantities": [{"id": "syy_D", "quantity": "σ_yy at D (2, 0)",
     "measure": {"command": "solve.probe", "params": {"step": "Load", "at": [[2, 0, 0]]}, "pointer": "/nodes/0/stress/1"},
     "value": 92.7, "tolerance": {"relative": 0.01}, "reason": "why this tolerance",
     "known_failure": "optional: why FieldCraft misses the reference today"}]}
  ```
  A measure reads a JSON pointer either from a recipe step's output (`"from": "<id>"`) or from a
  command run after the recipe (`"command"`, `"params"`); an optional `"scale"` multiplies it, an
  optional `"plus"` (a measure of its own) adds another number of the same run (the mean of two
  nodes, each scaled by ½), and an optional `"over"` (a measure of its own) divides the result by
  another number of the same run, for ratios such as the stress in one element over the stress
  in another.
  A case too long for every change (minutes) is marked `"nightly": true`: it runs in the nightly,
  manual and release runs (`--nightly`), or when named with `--case` or `--bless`; any other run
  checks it by its fixture (a fixture that misses an edited reference fails), and docs/oracle.md
  always shows its fixture, so the page is the same either way. `"slow"` cases run on every
  change and only `--quick` skips them.
  A tolerance is `relative` (of the reference), `absolute`, both, or `"exact": true`. A pending
  case cites its reference values (or `null` until transcribed) and has no recipe yet.
- `fixture.json` holds each passing quantity's accepted `value` and its `band`
  (10⁻⁹ of the larger of value and reference, plus 10⁻³ of an absolute tolerance; a quantity whose
  reference is zero, which is rounding noise, keeps its absolute tolerance as its band; an
  ill-conditioned quantity may set a wider `fixture_relative` with a `fixture_reason`; a quantity
  that depends on a mesh generated from geometry, which is not bit-identical across platforms, sets
  `"fixture": false` with a `fixture_reason` and is gated by its reference tolerance alone): the result is
  deterministic, the band only absorbs last-digit differences of the system maths library.

The run fails when a quantity misses its reference (unless it is a `known_failure`, which is
listed in oracle.md with its reason), when a known failure starts passing (remove the mark),
when a result leaves its fixture band, when a passing quantity has no fixture, or when a case
cannot run.

## Adding cases for a new capability

Every capability that lands brings its oracle cases in the same pull request:

1. Turn the pending cases that waited for it into active ones (write `recipe.json`, drop
   `status`/`requires`, add `measure` and `tolerance` with its reason to each quantity), and add
   theory cases for its closed-form solutions (patch test and convergence for an element; the
   analytic problem for a load, constraint or procedure).
2. Choose meshes that are converged well inside the tolerance (show the refinement in the case's
   README), or the mesh a benchmark prescribes.
3. Run `cargo xtask oracle --case <name>`; when it passes, `--bless <name>`; commit the case, its
   fixture and the regenerated docs/oracle.md.
4. A case that misses its reference stays honest: keep the tolerance its reason supports and add
   `known_failure` with what is missing. Never widen a tolerance to make a case pass.

## The meshing corpus

The meshers have a regression gate of the same kind: every part of the meshing test corpus
(`tests/corpus`, docs/meshing/corpus.md) has a `fixture.json` with each stage's outcome and its
deterministic measures, in the oracle's format, and every tet mesh the corpus makes is checked
against validity invariants that do not depend on the mesher's own gate. `cargo xtask corpus` and
the heavy test `corpus_matches_its_fixtures` (CI's heavy job) fail when a part moves from its
fixture or a mesh breaks an invariant; `cargo xtask corpus --bless <part>` rewrites a fixture.
docs/meshing/corpus.md gives the bands and the rules for blessing.

## TRT velocity inlet against published closures (#820)

FieldCraft's TRT velocity inlet (`crates/lbm/src/inlet.rs`; literature C3.3e) against the
published closures, each implemented from its paper behind the same inlet path, on the same three
cases; all TRT with Λ = 3/16, no turbulence model, measured on CPU (2026-11-13). The comparison
harness swapped the closure in place and is not committed; the committed tests cover the adopted
inlet.

- **Slip box**: the #820 16³ box, slip sides, uniform inlet 0.05 (Mach 0.087), convective pressure
  outlet, a 1e-6 checkerboard or broadband seed; steps survived of 30 000 (the run stops or the
  density departs by more than 0.1). Survival is checked in 50-step chunks; a solver error reports
  the last fully completed chunk.
- **Duct**: walled 48 × 16 × 16, uniform inlet 0.05, convective outlet, τ⁺ = 0.52, 60 000 steps
  in 1000-step chunks (a solver error reports the last fully completed chunk);
  over the last 10 000 the total-mass drift per step, the inlet's mass flux against `ρ_b U A`
  (the prescribed velocity), and inlet against outlet plane.
- **Poiseuille**: the parabolic-inlet channel of `regularized.rs`, τ⁺ = 0.6, H = 8, 16, 32 cells;
  the largest profile error over `u_max` at mid-length (order from 16 to 32) and in the first
  cell inside the inlet, at H = 32.

| Closure | Slip box, checkerboard, τ⁺ 0.505 / 0.507 / 0.51 / 0.52 | Slip box, broadband | Duct: steps; mass drift/step; flux vs `ρ_b U A`; in vs out | Poiseuille: mid error, order; next to inlet |
|---|---|---|---|---|
| FieldCraft (Zou–He + trace extrapolation + Shapiro high-pass) | all 30 000 (final ρ′ ≤ 8e-9) | all 30 000 (final ρ′ ≤ 3e-8) | 60 000; 1e-20; 8e-16; 7.9e-5 | 7.7e-5, 1.99; 5.6e-6 |
| Zou & He 1997 / Hecht & Harting 2010 (`develop`) | 2 800 / 12 200 / 30 000 / 30 000 | 400 / 450 / 600 / 5 700 | fails at 11 000 | 7.7e-5, 2.00; 6.8e-6 |
| Guo, Zheng & Shi 2002 extrapolation | 300 / 350 / 500 / 23 950 | 300 / 350 / 500 / 27 050 | fails before 1 000 | 3.1e-4, 1.98; 2.7e-4 |
| Latt et al. 2008 regularised | 400 / 550 / 950 / 30 000 | 450 / 550 / 1 000 / 30 000 | 60 000; 4e-20; 4e-16; 8.0e-5 | 7.7e-5, 2.00; 5.1e-5 |
| LODI (Izquierdo & Fueyo 2008; Heubes et al. 2014), bounce-back | 1 200 / 1 750 / 6 650 / 20 400 | 1 150 / 1 700 / 6 850 / 22 000 | 60 000; 0; −3.7e-4; 7.7e-5 | 7.6e-5, 1.97; 9.2e-5 |
| LODI, regularised populations | 300 / 400 / 550 / 11 250 | 300 / 350 / 500 / 11 500 | fails before 1 000 | 9.3e-5, 2.08; 1.5e-4 |

The in-versus-out mismatch of about 8e-5 is the same for every closure that survives: it is the
outlet plane's, not the inlet's (the inlet plane and the first interior plane agree to 1e-15). The
LODI inlets use the outgoing amplitude from one-sided interior differences with no pressure
relaxation; with bounce-back the boundary density is not the cell's, and the inlet's mass flux
falls 3.7e-4 short of `ρ_b U A`. A sponge or PML layer (Xu & Sagaut 2013; Najafi-Yazdi & Mongeau 2012)
was not built. Below τ⁺ = 0.505, reachable only with a turbulence model, none of these holds at
this Mach number (literature C3.3c, C1.5).

## Fluid capabilities at a glance

The fluid capabilities verified against theory, benchmarks or experiment, with their numbers and proposed
showcase examples, are listed in [What works](manual/fluid/what-works.md). Large-eddy simulation is not among
them: its measured errors, what blocks published accuracy and the options are in
[RFC 0003](rfc/0003-lbm-les-roadmap.md) (the LBM-LES work is parked as of 2026-11-26).

## Long studies

Some studies take hours and run in no CI job; the heavy job (`cargo xtask test --heavy`, 60
minutes) runs a bounded regression of each instead, and each study is run by hand when the
physics under it changes. The study's test stays `#[ignore]`d and returns at once unless its
environment variable is set, so `--heavy` passes over it.

### Re_τ 180 channel against Moser, Kim & Mansour (#663, #664)

Smagorinsky's model (C = 0.1, regularised collision), the channel 2πδ × 2δ × πδ driven by
u_τ²/δ, 10 δ/u_τ of start-up discarded and 15 averaged. Reference: R. D. Moser, J. Kim,
N. N. Mansour, *Direct numerical simulation of turbulent channel flow up to Re_τ = 590*, Phys.
Fluids 11:943–945, 1999 (Re_τ = 178): U_b/u_τ ≈ 15.7, U_c/u_τ ≈ 18.3, C_f ≈ 8.18 × 10⁻³.

```
FIELDCRAFT_CHANNEL_STUDY=1 cargo test --release -p fieldcraft-lbm --test turbulence -- --ignored --nocapture study
```

About 1.8 h on 32 threads. Recorded on Windows x86_64 (PR #668):

| cells per δ | y⁺ per cell | U_b/u_τ | U_c/u_τ | C_f |
|---|---|---|---|---|
| 12 | 15.0 | 9.95 (−36.6 %) | 13.07 (−28.6 %) | 2.02e-2 |
| 18 | 10.0 | 11.55 (−26.4 %) | 14.42 (−21.2 %) | 1.50e-2 |
| 24 | 7.5 | 12.81 (−18.4 %) | 15.60 (−14.8 %) | 1.22e-2 |
| 32 | 5.6 | 14.17 (−9.8 %) | 17.00 (−7.1 %) | 9.96e-3 |

The heavy job runs the 12-cell grid alone with the same protocol
(`a_coarse_channel_at_re_tau_180_stays_turbulent_with_the_study_s_means`, about a minute and a
half on 14 cores): turbulent, and within 3 % of the 12-cell means above. The flow is chaotic, so
its means move with the platform (9.90 and 12.95 on macOS aarch64); the band is three times that
scatter.

**With Van Driest damping** (`van_driest`: 26, literature C5.3; probes, not validation). At
Re_τ 180 a grid of N cells per δ has Δy⁺ = 180 / N, so none of these grids resolves the wall:
damping is meant for a wall resolved down to y⁺ ≈ 1–2 at the first cell (about 90 cells per δ).
The probes are short (5 δ/u_τ of start-up discarded, 1 averaged) and the bulk speed was still
falling at their end. Recorded on Linux aarch64 (NVIDIA GB10, 20 cores, shared):

| cells per δ | Δy⁺ | damped, short run: U_b/u_τ | U_c/u_τ | undamped, full protocol: U_b/u_τ |
|---|---|---|---|---|
| 12 | 15.0 | unstable at 2.2 δ/u_τ | — | 9.95 |
| 18 | 10.0 | 12.42 (−20.9 %) | 15.63 (−14.6 %) | 11.55 |
| 24 | 7.5 | 14.00 (−10.8 %) | 17.26 (−5.7 %) | 12.81 |

All three are **not wall-resolving**. On the 12-cell grid the damping takes away the eddy
viscosity that held τ₀ = 0.5013 stable at a wall the grid cannot resolve; the run now says so
(the first cells' y⁺ in the instability's message, and a warning past y⁺ 5 when a run ends).
A wall-resolving study in a minimal flow unit (πδ × 2δ × 0.3πδ at 90 cells per δ, 4.3 M cells,
about 6 × 10¹⁰ cell updates per δ/u_τ) would take a day or more on the CPU and is not run.

What it does not establish yet:
- **No converged reference.** The bulk speed approaches MKM's monotonically as the first cell
  moves into the viscous sublayer, but no grid is within its tolerance of it: without wall
  damping the eddy viscosity is too large in the first cells off the wall, so the friction is too
  high, and the grids above are too coarse for damping to model the wall. A wall model
  (literature C5.2) or a wall-resolving grid is the next step.
- **The finest grid had not settled.** At 32 cells per δ the bulk speed was still falling slowly
  at the end of the averaging, so its mean is an upper bound of that run's, and the observed
  trend is not a convergence order.
- **One realisation per grid.** The means carry the sampling scatter of a single 15 δ/u_τ
  average (about 1 % at 12 cells per δ, from the two platforms above); no confidence interval is
  stated.

### Spalart–Allmaras: wall-resolved channel against its 1-D solution and Moser, Kim & Mansour

The Spalart–Allmaras model (no wall law) in a channel at Re_τ 178.12, one cell along x and z (the
flow is fully developed), halfway walls, the body force u_τ²/δ, the lattice viscosity fixed at
0.002 (τ⁺ 0.506) so that Δt ∝ Δx² (`crates/lbm/tests/sa_channel.rs`).

- **Theory (asserted).** The steady 1-D SA solution with the same constants, solved in the test by
  second-order finite volumes and pseudo-transient Newton, Richardson-extrapolated from 1728 and
  3456 cells (U_b⁺ 15.850273, U_c⁺ 18.447753, largest correction 7.9e-6 u_τ). The lattice u⁺
  profile must converge to it at second order: the observed order on the finest pair in [1.8, 2.2],
  `sa_mms.rs`'s band (heavy test). The trunk test checks the reference's own order and that the
  coarsest lattice run settles on the turbulent branch.
- **Benchmark (reported, not asserted).** R. D. Moser, J. Kim, N. N. Mansour, Phys. Fluids
  11:943–945, 1999, `chan180.means`: U_b⁺ 15.68, U_c⁺ 18.301, C_f 8.135 × 10⁻³. No published
  SA-against-DNS comparison at this Re_τ gives a tolerance, so the gap is reported as the model's
  error: the converged 1-D SA solution is +1.09 % in U_b, +0.80 % in U_c and −2.1 % in C_f.

| cells across 2δ | first cell y⁺ | L2 / L∞ of u⁺ against 1-D | U_b⁺ | U_c⁺ | C_f | U_b, U_c, C_f against MKM |
|---|---|---|---|---|---|---|
| 96 | 1.86 | 0.130 / 0.138 | 15.726 | 18.310 | 8.087e-3 | +0.30 %, +0.05 %, −0.59 % |
| 144 | 1.24 | 0.0592 / 0.0628 | 15.794 | 18.385 | 8.018e-3 | +0.73 %, +0.46 %, −1.43 % |
| 216 | 0.82 | 0.0265 / 0.0280 | 15.825 | 18.420 | 7.986e-3 | +0.92 %, +0.65 %, −1.82 % |
| 1-D (converged) | — | — | 15.850 | 18.448 | 7.961e-3 | +1.09 %, +0.80 %, −2.1 % |

Observed orders (L2, L∞): 1.937, 1.941 (96→144) and 1.982, 1.996 (144→216). The coarse grids sit
closer to MKM only because their discretisation error has the opposite sign to the model's.
The heavy study took 245 s.

### Spalart–Allmaras: its wall law in the channel against the wall-resolved 1-D solution and Lee & Moser

The Spalart–Allmaras model with its own wall law (`wall_model: "spalart_allmaras_bounce"`) in a column
1 × (2N + 1) × 1, δ = N cells, periodic in x and z, on-node walls, the body force u_τ²/δ, the regularized
collision (TRT keeps the 0.505 floor with this law and cannot run these grids, #1063), lattice bulk speed 0.1
(`crates/lbm/tests/sa_wall_channel.rs`). Re_τ 182.088 at N 3 and N 6, and Re_τ 1000.512 at N 10, 16 and 32, each
run to 1e-6 u_τ. The 1-D solutions are Richardson-extrapolated from the reference grid (1728 or 3200 cells) and
its double (the Re_τ 1000.512 solve starts from the Re_τ 182.088 solution scaled in wall units near the wall and
in outer units above, χ = min(κ y⁺, χ₁₈₂ Re_τ/182.088)), and the wall-modelled 1-D solution starts at y₁ = Re_τ / N with the law's stress share.

- **Theory (asserted).** (1) The 1-D solutions are turbulent and converge at second order at shared faces, wall-resolved
  and wall-modelled (trunk test, and heavy test for Re_τ 1000.512). (2) The momentum balance: the walls' force equals
  a × the fluid mass, within the run's last change Δ, since per δ/u_τ the momentum moves at most M Δ u_τ. (3) The
  walls' friction velocity is √(a (N − ½)), within Δ: the boundary nodes take no force, so each wall carries the
  share 1 − 1/(2N) (#1101). (4) The first node: U₁ = u_τ u⁺(u_τ y₁ / ν) and ν̃₁ = κ u_τ · 1 cell, exactly. (5) The
  trapezoid identity u_{j+1} − u_j = ½ (g_j + g_{j+1}), g_j = a (N − j) / (ν + ν_t,j), on every link, within Δ u_τ²/(ν + ν_t,1)
  (derived in the test's module docs). Each bar uses the run's own Δ.
- **Benchmark (reported, not asserted).** Lee & Moser, J. Fluid Mech. 774:395–415, 2015: U_b⁺ 1/0.0637309 = 15.69 (Re_τ 182.088),
  U_b⁺ 1/0.0500256 = 19.99 and U_c⁺ 22.593 (Re_τ 1000.512).

| Re_τ | N per δ | y₁⁺ | U_b⁺ | U_c⁺ | C_f | against wall-resolved 1-D: U_b, U_c, C_f | of which the law alone: U_b, C_f | lattice − law (as run), L∞ (u_τ) | against Lee & Moser: U_b, U_c, C_f |
|---|---|---|---|---|---|---|---|---|---|
| 182.088 | 3 | 60.7 | 13.390 | 15.458 | 1.1155e-2 | −15.79 %, −16.38 %, +41.0 % | −6.49 %, +14.4 % | 6.5e-2 | −14.66 %, —, +37.3 % |
| 182.088 | 6 | 30.3 | 14.567 | 17.006 | 9.4253e-3 | −8.38 %, −8.01 %, +19.1 % | −4.13 %, +8.8 % | 4.2e-2 | −7.16 %, —, +16.0 % |
| 1000.512 | 10 | 100.1 | 19.163 | 21.410 | 5.4463e-3 | −3.45 %, −3.21 %, +7.3 % | −1.06 %, +2.2 % | 2.2e-2 | −4.14 %, −5.24 %, +8.8 % |
| 1000.512 | 16 | 62.5 | 19.478 | 21.741 | 5.2718e-3 | −1.87 %, −1.71 %, +3.9 % | −0.61 %, +1.2 % | 5.7e-2 | −2.56 %, −3.77 %, +5.3 % |
| 1000.512 | 32 | 31.3 | 19.773 | 22.050 | 5.1155e-3 | −0.38 %, −0.32 %, +0.8 % | −0.28 %, +0.6 % | 1.3e-1 | −1.09 %, −2.41 %, +2.2 % |

The 1-D wall-resolved solutions: Re_τ 182.088 U_b⁺ 15.899968, U_c⁺ 18.485491 (correction 8.3e-6); Re_τ 1000.512
U_b⁺ 19.848774, U_c⁺ 22.119259 (correction 8.6e-5). Each lattice run settled after 108 to 147 δ/u_τ with
Δ ≈ 1e-7 u_τ and 1e-6 u_τ left.

The column "lattice − law (as run)" compares the run with the continuum wall-modelled solution at the same
stress. It is reported, not asserted, and it grows with N: in lattice units the log profile is ln(y/Δx) whatever N,
and the trapezoid sum of 1/(κ j) over the rows exceeds ln j/κ by (γ − ½)/κ = 0.19 u_τ as j grows (γ is Euler's
constant), so the more rows lie in the log layer, the closer the offset comes to that value (2.2e-2, 5.7e-2 and
1.3e-1 u_τ at N 10, 16 and 32). SA's own discrete ν̃ in the first rows moves it by a few hundredths.

- The walls carry a(δ − Δx/2) (#1101), not aδ: the law gets 1 − 1/(2N) of the stress. That half cell is most of
  the gap between the law alone and the run (at Re_τ 1000.512, C_f +1.2 % and +2.2 % for the law alone against
  +3.9 % and +7.3 % as run, first node at y⁺ 62 and 100).
- The law's own error is negative and large in these cases: in the wall-resolved solution ν̃ falls below κ u_τ y
  faster than the stress, so the law under-predicts U_b (−1.06 % at Re_τ 1000.512, N 10, the first node at y⁺ 100).
- **Limitation: the law assumes a log layer of constant stress.** At Re_τ 182.088 the first node at y⁺ 30 or 61 is at
  y/δ 1/6 or 1/3, where the stress has fallen by that share and there is no constant-stress layer: C_f is +19 % and
  +41 % against the wall-resolved model (the law alone +8.8 % and +14.4 %). The law is for first nodes well inside a
  log layer, y₁/δ small; at Re_τ 1000.512 with y₁/δ 1/32 to 1/10 the error is +0.8 % to +7.3 %.
- y⁺ 100 at Re_τ 182.088 is not run: it needs fewer than 2 cells per δ. Re_τ 2000 is not run in this study: its
  1-D reference was not solved (from the plain start the Newton iteration does not converge within its cap; the
  start scaled from the Re_τ 182.088 solution, which the Re_τ 1000.512 reference uses, was not tried there).
- **The first rows carry the stress.** Wall-modelled LES loses the stress of the first rows (the resolved part is
  missing there, so the law must be read three or more cells up; Kawai & Larsson 2012). In RANS the model's own
  ν_t carries it: the modelled share (1 − ω_j/2)(−Π_xy)/(ρ a (N − j)) of rows 1 to 3, ω_j with the cell's eddy
  viscosity, is 1.000000 in every case, as the link identity makes it. Reading the law higher would add error here,
  since the law's own error grows with y₁/δ. (The finite-difference readings (ν + ν_t) dU/dy of 0.73 to 1.22 are the
  difference quotient's error at one to three cells, not the stress.)
- Setup warns that τ⁺ 0.5001 to 0.5006 is below 0.508 (the regularized collision's stability note); the bars hold.

### Wall-modelled LES: Malaspinas & Sagaut at Re_τ 950 (void)

**Void (#1057).** The MRT collision is linearly unstable at rest below a lattice viscosity of 1.1801 × 10⁻³
(τ⁺ 0.50354; [Collision operators](manual/fluid/collision.md#the-stability-bound),
`crates/lbm/tests/mrt_stability.rs`). These runs were at 5.4 × 10⁻⁵, 1.08 × 10⁻⁴ and 1.62 × 10⁻⁴ (N 10, 20, 30):
22, 11 and 7 times below it. `solve.fluid` now refuses them. The table is kept as the record of what was run.

Malaspinas & Sagaut's channel run exactly as their paper sets it up (Re_τ 950, Re_m 37 042, N 10,
20 and 30, Smagorinsky C 0.08, MRT, their forcing), against the paper's own measures and Hoyas &
Jiménez's DNS. Set-up, the readings where the paper is silent, and figures:
[Large-eddy simulation](manual/fluid/les.md). Recorded on Linux aarch64 (NVIDIA GB10), GPU, f32.
Pass criteria fixed before the runs: E within 15 % of the paper's Table 1 and C_f by the model's
u_τ within its 4.22 % of Dean.

| N | E (paper) | C_f by the model's u_τ vs Dean | C_f by the carried drag vs Dean |
|---|---|---|---|
| 10 | 2.08 (2.11) | +1.8 % | −16.8 % |
| 20 | 0.99 (1.07) | +0.8 % | −13.7 % |
| 30 | 0.87 (0.83) | +3.5 % | −14.1 % |

Recorded as passing on the paper's measures at the time; void, as above. The 14–17 % low carried drag is not established as a property of the method.

### Requested, applied and carried wall stress

The shared test support in `crates/lbm/tests/support/wall_balance.rs` separates
requested Musker shear, applied momentum-exchange traction (including reconstruction
impulse exactly once), and carried traction inferred from body impulse minus the
change in **stored population momentum**, divided by the same wall area. Displayed
wall-target velocities cannot replace stored momentum in this budget.

For a closed, unheated box with uniform acceleration and no zones or on-node
boundary nodes, the per-step identity is
`P_after - P_before = body_impulse - wall_impulse`. The residual limit is
`relative_tolerance * scale + 64 * f64::EPSILON * max(lattice_mass, 1)`, where
`scale` is the maximum component magnitude of the two stored momenta and the
body/wall impulses. Relative tolerances are `1e-12` on CPU and `1e-5` on GPU;
the absolute floor covers the host sums on both backends, including rest.
Tests cover both streaming layouts and odd/even GPU steps. Open boundaries need
momentum fluxes and are outside this diagnostic's scope.

Carried drag is compared independently with DNS at 5%; requested/carried stress
has a proposed 1% agreement gate. The original Re_τ 950 N10 audit fails both:
requested/carried ≈1.245, carried C_f 0.0043803 versus DNS 0.00509 (−16.75% versus
Dean, reported separately). That audit's run used MRT below its stability bound (above), so the two numbers are frozen data of that run, not a measurement of the wall model; the test checks the classification of frozen data only. The fast frozen-data regression
`measured_n10_discrepancy_remains_a_failed_validation` asserts both failures and
fails if either starts passing; a changed validator must not silently accept
this known discrepancy.

Until the #785 wall repair, the long M&S study **reports** the additional
carried-drag/stress-agreement result without adding it to its failure list, as
requested in the #829 review. Its existing paper-reproduction and in-range DNS
acceptance checks remain active. This keeps those regression checks usable without
calling the known stress defect a validation pass. The frozen-data test checks
classification, not a fresh flow result; after the repair, rerun the study and
review its diagnostic gates before claiming physical validation.

### Statistical acceptance of study windows

The M&S and Han studies stop on `fieldcraft_lbm::statistics::Blocks` (#796). A burn-in in samples
precedes equally spaced samples; partial blocks never count. Acceptance is assessed on the largest
even prefix of complete blocks, so both block lengths cover the same window, and needs:

- at least 40 base blocks (`MIN_BLOCKS`), so the doubled-length series has 20 values;
- a Student-t 95 % mean half-width below 1 % of the absolute mean at both block lengths, and a
  change of less than 20 % in that half-width when the block length doubles;
- a half-window drift below 1 %;
- absolute block correlations below 0.5 at lags 1 to 5 on both block lengths.

The block minimum follows from the lags tested: a lag-k autocorrelation from N independent values
has a standard error near 1/√N (NIST/SEMATECH §1.3.5.12, after Bartlett 1946), and Box & Jenkins
(1976, §2.1.6) estimate lags only up to about N/4, so lag 5 on the doubled series needs 20 doubled
blocks. At the earlier minimum of 20 base blocks (10 doubled values, standard error about 0.3) a
simulation of 20 000 independent Gaussian windows rejected 47 %: 18 % on the doubled-series
correlations and 34 % on the half-width change. At 40 blocks 14 % are rejected and at 80 blocks
1.5 %; a rejected window only means the study keeps averaging. The unit test
`independent_windows_at_the_minimum_usually_pass` checks the same on 400 synthetic windows.

The interval uses the [NIST/SEMATECH Student-t definition](https://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm)
with its [0.975 critical values](https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm),
rounded upward; above 30 degrees of freedom the value for 30 is used (conservative). A zero mean
cannot pass a relative interval. The store is capped at 65 536 block means (512 KiB); a sample
beyond the cap is an error and leaves acceptance false.

Rejection controls: three FieldCraft cylinder drag histories from #796 (Re 3900, N12, HRR, GB10;
100–200 and 400–600 D/U, and 400–600 D/U from a potential-flow start) are replayed at blocks of
1.2 and 2.5 D/U (at least 40 blocks each) and must keep failing on the evidence: their relative
half-widths are 2.7–9.2 % and their drifts 1.8–8.8 %. They are failure controls, not validation
fixtures. The synthetic unit tests reject alternating, lag-2, drifting and pair-dependent signals.

### Wall-modelled LES: Han, Ooka & Kikumoto not reproduced at A20

Han, Ooka & Kikumoto's channel (Fluid Dyn. Res. 53:045506, 2021) at A20 (Re_τ 640), with their
wall-function bounce (eq. 14, on-node) and with plain bounce-back, on the paper's MRT collision
at C_s 0.1 (void: below the collision's stability bound) and on the hybrid recursive regularised collision. Gate fixed before the runs: ε_u⁺
over the whole channel and the log-law intercept b within 15 % of their Tables 2 and 3. Recorded
on Linux aarch64 (NVIDIA GB10), GPU, f32. Set-up, the flow's diagnosis (the MRT collision's instability at this viscosity),
the odd–even mode of the wall rows, figures and
how to rerun: [Large-eddy simulation](manual/fluid/les.md).

| A20, wall-function bounce | U_b⁺ | C_f vs Dean | ε whole (paper) | b (paper) |
|---|---|---|---|---|
| MRT, Smagorinsky C_s 0.1, Van Driest (the paper's): **void**, lattice viscosity 5.0 × 10⁻⁵, 24 times below the bound | 14.44 | +53 % | 0.248 (0.057) | 1.92 (6.4) |
| hybrid σ 0.99, Vreman C_s 0.18 | 28.77 | −54 % | 0.535 (0.057) | 15.23 (6.4) |
| hybrid σ 0.99, Smagorinsky C_s 0.1, Van Driest | 27.94 | −52 % | 0.491 (0.057) | 14.30 (6.4) |
| hybrid, dynamic σ, Vreman C_s 0.18 | 28.94 | −55 % | 0.607 (0.057) | 9.04 (6.4) |

**Not reproduced** on the hybrid collision, which carries about half the right drag: the first rows above the
wall cannot carry the stress (the log-layer mismatch). **The MRT row is void**: at A20's lattice viscosity the
MRT collision is linearly unstable at rest, without walls (growth 7.9 × 10⁻³ per step), and the channel's
u′ ≈ v′ ≈ w′ ≈ 5 u_τ is that instability's saturated state, which a periodic box at rest reaches too
(4.4–4.8 u_τ). The earlier reading, "driven at the walls" and "the documented failure in under-resolved
channels", is withdrawn. The paper states neither its MRT rates nor its damping constant, so its result is not
reproducible from what is published. A40 was not run.

### Bluff-body LES (#785)

> Probe runs of 2026-11-23/24 on develop 655d92fe (NVIDIA GB10, f32; hybrid collision σ 0.99, Smagorinsky C 0.1, no
> wall model; box 24 D long, body 6 D from a uniform inlet, convective outlet, slip sides, periodic span; averaged
> over 100–240 D/U in 8 blocks). They are not acceptance runs: no grid is converged and the intervals are block
> intervals of one realisation. References: circular cylinder at Re 3900, C_d 0.99 ± 0.05 (experiment) and 1.04
> (LES), A. G. Kravchenko, P. Moin, Phys. Fluids 12:403–417, 2000, doi:10.1063/1.870318; square cylinder, C_d 2.1 at
> Re 21 400, D. A. Lyn, S. Einav, W. Rodi, J.-H. Park, J. Fluid Mech. 304:285–319, 1995,
> doi:10.1017/S0022112095004435, and 2.18 (DNS at Re 22 000), F. X. Trias, A. Gorobets, A. Oliva, Comput. Fluids
> 123:87–98, 2015, doi:10.1016/j.compfluid.2015.09.013.

| case | start | obstacle walls | C_d | against the reference |
|---|---|---|---|---|
| cylinder Re 3900, 16 cells per D, blockage 8.3 % | free stream everywhere | halfway | 1.310 ± 0.043 | +32 % |
| same | at rest on the body | halfway | 1.260 ± 0.012 | +27 % |
| same | at rest on the body | interpolated | 1.113 ± 0.012 | +12 % |
| same, blockage 4.2 % | at rest on the body | interpolated | 1.098 ± 0.018 | +11 % |
| same, 24 cells per D, blockage 8.3 % | at rest on the body | interpolated | 1.067 ± 0.019 | +8 % |
| square Re 22 000, 16 cells per D, blockage 7.1 % | at rest on the body | halfway (faces on cell faces) | 1.81 ± 0.07 | −14 % (Lyn), −17 % (Trias) |
| same, blockage 3.6 % | at rest on the body | halfway | 1.73 ± 0.03 | −18 %, −21 % |

> What the earlier records of these cases (#787: 1.32–1.35; #789: 4.11) contained besides resolution:
> - **A single-step force.** Near τ⁺ = ½ a halfway wall's force alternates step by step. On the square cylinder the
>   odd steps alone read 4.09 and the even steps alone −0.48; their mean, 1.81, is the pressure integrated over the
>   body's faces (1.79). #789's 4.11 was one of the two.
> - **The start.** With the free stream filling the box, body included, the body appears impulsively and the
>   pressure wave it sends stays between the reflecting inlet and outlet: on the cylinder the force's rms is 3.7
>   against 0.41 with a start at rest on the body, and the mean is 4 % higher; the square case went unstable at
>   161 D/U.
> - **The staircase**, 13 % on the cylinder (rows 2 and 3).
>
> The remaining error is resolution: 16–24 cells per diameter leave the separating shear layers on one or two cells.
> Density stays within 0.84–1.18 (square) and 0.93–1.08 (cylinder) of the reference, local Mach numbers up to 0.28.

### Wall-modelled channel at Re_τ 1000 on the regularised collision (#785)

> Probe runs of 2026-11-23/24 on develop 655d92fe: 2πδ × 2δ × πδ, constant force u_τ²/δ, Smagorinsky C 0.1, the
> regularised collision unless stated; 15–20 δ/u_τ discarded and 15–20 averaged, shorter than the studies'
> protocol, so every row is inconclusive by the statistical acceptance. Reference: M. Lee, R. D. Moser, J. Fluid
> Mech. 774:395–415, 2015 (U_b⁺ 19.99 at Re_τ 1000.5).

| wall treatment | cells per δ | U_b⁺ | C_f against Lee & Moser |
|---|---|---|---|
| none (halfway bounce-back) | 16 | 10.5–12.1 | +173 to +261 % |
| Malaspinas & Sagaut (`musker`) | 16 | 28.6–33.1, still rising | −51 to −63 % |
| same | 32 | 35.5, still rising | −68 % |
| same, hybrid collision σ 0.99 | 32 | 38.9, still rising | −74 % |
| Han, Ooka & Kikumoto (`spalding_bounce`), Van Driest 26 | 16 / 32 | 27.3 / 28.1 | −46 / −49 % |
| Han's wall, the law sampled at the 2nd row (probe) | 16 | 21.5 | −13.9 % |
| sampled at the 3rd row, h/δ 0.19 (probe) | 16 | 19.4 | +5.8 % |
| sampled at the 4th row, h/δ 0.17 (probe) | 24 | 19.5 | +5.4 % |
| sampled at the 6th row, h/δ 0.19 (probe) | 32 | 19.5 | +5.1 % |

> Both published forms read the law of the wall in the first cells, where the grid resolves none of the stress
> (resolved stress 0 to 0.5 of the total in rows 1–2), so the wall takes too little and the flow speeds up. Reading
> it three or more cells up (S. Kawai, J. Larsson, Phys. Fluids 24:015105, 2012) brings the carried friction to
> about +5 %; the probe changed only where the law is sampled.

### Matched on-node wall at Re_τ 1000 against Lee & Moser (#1105, #785)

`wall_model: "spalding_bounce_matched"`: Han, Ooka & Kikumoto's on-node wall with Spalding's law read four cells from
the wall and Malaspinas & Sagaut's eq. 57 mixing-length eddy viscosity on the rows below. The channel is 2πδ × (2δ + Δx) ×
πδ, driven by u_τ²/δ (u_τ = δ = 1), regularised collision, Smagorinsky C 0.1, lattice bulk speed 0.1, started from a perturbed
blunt profile; 20 δ/u_τ of start-up discarded and at least 50 averaged, stopped when
[`statistics::Blocks`](#statistical-acceptance-of-study-windows) accepts the bulk speed (one block per δ/u_τ), a cap of 400 δ/u_τ
being inconclusive. Reference: M. Lee, R. D. Moser, J. Fluid Mech. 774:395–415, 2015 (Re_τ 1000.5, U_b⁺ 19.99; the profile at
16 points). Bar, fixed before the runs: C_f, (U_b⁺_LM / U_b⁺)² − 1, within 5 %, and U⁺ within 5 % on every row from the sample
row (the rows at least 4 cells from the nearer wall) to the centreline. U_b⁺ is the trapezoid rule over the full 2δ, the walls'
zero at y = 0 and 2δ included, as Lee & Moser's U_b covers 2δ.

The sample height is a fixed 4 cells, our choice, not the paper's number for this method. Kawai & Larsson (Phys. Fluids
24:015105, 2012) state two conditions: sample away from the cells next to the wall, which the LES resolves poorly (§III), and
h_wm fixed in outer units, about 0.1δ, inside the inner layer (§III.B); their h_wm ≈ 4Δy_w is for their sixth-order compact
scheme (Δy_w = 0.011δ). Four cells address the first. The second holds only when 4Δx ≤ 0.1δ, about 40 cells per δ or more;
four cells are 0.25δ, 0.167δ and 0.125δ at 16, 24 and 32 cells per δ, and coarser grids carry the error this study measures.
The study prints h_wm/δ for its grid and marks it above 0.1. The model carries no channel half-height, so setup does not warn. In the laminar limit the
law read at 4Δx gives a slip of a(4Δx)²/(2ν) and a bulk high by 96/(M(M+1)) with M = 2N wall-to-wall cells: 35 %, 9.1 % and
2.3 % at 8, 16 and 32 cells per δ (`the_matched_wall_slips_by_the_law_at_four_cells_in_a_laminar_channel`, 35.3 % measured at 8).

```
FIELDCRAFT_MATCHED_CHANNEL=16 cargo test --release -p fieldcraft-lbm --test turbulence -- --ignored --nocapture matched_wall_channel
```

The heavy job runs a bounded regression of it (`a_coarse_matched_wall_channel_at_re_tau_1000_stays_stable_over_a_short_run`,
about four minutes on the build server): 8 cells per δ, 2 δ/u_τ discarded and 2 averaged, the study's code. It asserts a
stable run (largest local Mach below 0.3), finite positive row means and U_b⁺ within 30 % of Lee & Moser's; it is not the
accuracy bar. The study's table (grids 16, 24, 32) is not recorded yet: the runs take hours.

### Wall-modelled LES: Malaspinas & Sagaut against Lee & Moser and Dean (#664, #749)

**Void (#1057).** Every run of this section used the MRT collision below its stability bound (lattice viscosity
1.1801 × 10⁻³): Re_τ 1000 at 5.1 × 10⁻⁵, 1.02 × 10⁻⁴ and 1.53 × 10⁻⁴ (10, 20, 30 cells per δ); Re_τ 2000 at
4.6 × 10⁻⁵ and 6.9 × 10⁻⁵ (20, 30); Re_τ 550 at 1.0 × 10⁻⁴ to 1.8 × 10⁻⁴. They are 7 to 26 times below it. The
tables are kept as the record of what was run; none of the numbers, the Re_τ 2000, 30-cell agreement included,
is a result, and the Kawai & Larsson reading of the 9–15 % is not established.

Malaspinas & Sagaut's wall model with Musker's profile and the multiple-relaxation-time collision
([Large-eddy simulation](manual/fluid/les.md)), on the GPU in f32. Periodic channel, half-height δ, with the flow rate
held at the bulk velocity of Dean's correlation by a body force updated every step; 20 δ/u_τ of
start-up discarded and at least 50 averaged; Smagorinsky's C 0.08 at Re_τ 1000 and 0.1 at
Re_τ 2000. References: M. Lee, R. D. Moser, *Direct numerical simulation of turbulent channel flow
up to Re_τ ≈ 5200*, J. Fluid Mech. 774:395–415, 2015 (doi:10.1017/jfm.2015.268), and R. B. Dean's
correlation (J. Fluids Eng. 100(2):215–223, 1978, doi:10.1115/1.3448633). Acceptance, fixed before
the runs: C_f and the log layer's mean velocity within 5 % of Lee & Moser.

```
FIELDCRAFT_PAPER_CHANNEL=1 cargo test --release -p fieldcraft-lbm-gpu --test verification \
  -- --ignored study_malaspinas_and_sagaut_s_channel_against_lee_and_moser --nocapture
```

This study uses the MRT collision at the paper's relaxation times, which `solve.fluid` now refuses (below its
[stability bound](manual/fluid/collision.md#the-stability-bound)); the command records how the void runs were made.

Two minutes to an hour per grid. Recorded on Linux aarch64 (NVIDIA GB10), PR #749. C_f by the
force balance is the drag the flow carries, which equals the body force at steady state and is
what `solve.fluid` reports as the wall force; C_f by the wall model's u_τ is the one the model
infers at the second cell, Malaspinas & Sagaut's measure. E is their profile error measure (their
values with Musker's profile in brackets: Re_τ 950 and 2000).

| Re_τ | cells per δ | first cell y⁺ | C_f, force balance: vs Lee & Moser / vs Dean | C_f, wall-model u_τ: vs Dean | E (paper) | log layer vs Lee & Moser |
|---|---|---|---|---|---|---|
| 1000 | 10 | 46 | −11.4 % / −14.5 % | +1.3 % | 2.65 (2.11) | 7.6 % |
| 1000 | 20 | 23 | −8.9 % / −12.1 % | −0.9 % | 1.21 (1.07) | 6.7 % |
| 1000 | 30 | 16 | −8.7 % / −11.9 % | +1.4 % | 0.99 (0.83) | 8.3 % |
| 2000 | 10 | 100 | refused: τ 0.50007 is below the 0.5001 floor | — | — | — |
| 2000 | 20 | 52 | +10.2 % / +8.9 % | −3.1 % | 1.01 (1.15) | 9.0 % |
| 2000 | 30 | 33 | **−0.7 % / −1.8 %** | −1.8 % | 0.29 (0.88) | **3.3 %** |

- **C_f by the wall model's u_τ** is within the paper's 0.3–4.2 % of Dean on every grid. The E
  values here took x₁'s velocity before its rebuild and are superseded by the Re_τ 950
  reproduction above. A CPU run in f64 of Re_τ 1000 at
  10 cells per δ gives −14.8 % against Dean, against the GPU's −14.5 %.
- **The drag the flow carries meets the acceptance only at Re_τ 2000 with 30 cells per δ.** On
  the other grids it is 9–15 % off, the log-layer mismatch of S. Kawai, J. Larsson, *Wall-modeling
  in large eddy simulation: length scales, grid resolution, and accuracy*, Phys. Fluids 24:015105,
  2012 (doi:10.1063/1.3678331). Recorded as it is, not counted as a pass (owner decision
  2026-10-26); the next step is sampling the wall model's input several cells off the wall, as
  Kawai & Larsson do, on the same grids.
- At Re_τ 1000 the 20- and 30-cell grids put the first cell below y⁺ 30, outside the range the
  method is meant for; they are listed for completeness.

**Re_τ 550: out of range.** Re_τ 550 was in the original protocol, but its grids put the first cell
centre at y⁺ 13–24, below the method's range. This was found after the runs, so it is recorded here
as out of range, with its results, rather than counted. The protocol chose those grids by cell size
(Δy⁺ 54, 39 and 30) where the method's range is for the first cell's centre, half a cell off the
wall, at y⁺ 30 or more: that was our protocol's error, recognised on 2026-10-25 after the results
below. Malaspinas & Sagaut validate only at Re_τ 950 and above, so 550 is outside the paper's
tested range too; the y⁺ placement decides it. Owner decision 2026-10-26: keep its numbers, not in
the acceptance.

These runs used the study as it stood on 2026-10-25, before it followed the paper exactly:
2πδ × 2δ × 2πδ, MRT, C 0.08, the wall model with **Spalding's** profile (not Musker's), the flow
rate held at **Lee & Moser's** bulk velocity and updated every 0.05 δ/U_m (not every step), and the
running-mean stop rule on the force. Recorded on Linux aarch64 (NVIDIA GB10), GPU, f32. The current
study lists Re_τ 550 at 10, 14 and 18 cells per δ outside the acceptance (`FIELDCRAFT_PAPER_CHANNEL=550`);
that version has not been run at 550. The C_f by the wall model's u_τ, E and C_f against Dean were
not measured for these runs.

| cells per δ | measured Re_τ | first cell y⁺ | U_b⁺ (Lee & Moser 18.40) | C_f, force balance vs Lee & Moser | log layer vs Lee & Moser | cells | GPU time |
|---|---|---|---|---|---|---|---|
| 10 | 470 | ≈ 23 | 21.30 | −25.3 % | 17.9 % | 79 k | 10 s |
| 14 | 461 | ≈ 16 | 21.72 | −28.3 % | 23.2 % | 217 k | 37 s |
| 18 | 475 | ≈ 13 | 21.08 | −23.8 % | 19.0 % | 460 k | 121 s |

The first cell's y⁺ is the measured Re_τ × 0.5 / N.

### What #664's acceptance still lacks

- **Developing pipe flow against Shah & London (1978)**: not built.
- **The graphics-card cooling case (#642) to statistically steady state**: not run.
- **The stop rule on a turbulent flow.** It is shown on a physically fluctuating laminar one:
  vortex shedding behind the Schäfer–Turek cylinder at Re 100
  (`vortex_shedding_stops_on_its_running_mean_with_the_long_run_s_mean_drag` in
  `crates/lbm/tests/steady.rs`, heavy job). The rule on the instantaneous force runs all 16 s;
  on the running mean from 4 s it stops at 9 s with the mean drag 0.24 % from the long run's
  (tolerance 1 %) and the lift's RMS 0.8 % from it. A turbulent case (the channel, or #642's
  card) stopping on its running mean is not run.

Done for #664: the averaged statistics in the solve reply (`stop.averages`), the stop rule above,
model-off results unchanged, and Poiseuille and Taylor–Green converging like Δ² with the model on
(`crates/lbm/tests/turbulence.rs`).
