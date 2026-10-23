# Verification and validation

Every solver capability is checked against three kinds of reference, and every accepted result is
frozen so it cannot drift.

| Kind | What it checks | Sources |
|---|---|---|
| Theory | The code solves the equations right | Closed-form solutions (beams, plates, Lamé, Kirsch, waves, shocks, Poiseuille…), patch tests, convergence rates, mesh invariants |
| Experiment | The model predicts reality | Measurements published in papers and reports (Taylor bar, cavity and cylinder flows, FSI experiments…) |
| Published benchmarks | Agreement with established codes on standard problems | NAFEMS benchmarks, vendor benchmark and verification manuals (e.g. Abaqus Benchmarks and Verification manuals), other published reference solutions |

Rules:
- Each case is **our own model** of the problem, built from the published description (geometry,
  material, loads, boundary conditions). We never run a commercial product and never copy its
  input files, manual text or figures; only the published numbers are used, with a citation.
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
  command run after the recipe (`"command"`, `"params"`); an optional `"scale"` multiplies it, and
  an optional `"over"` (a measure of its own) divides it by another number of the same run, for
  ratios such as the stress in one element over the stress in another.
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

What it does not establish yet:
- **No converged reference.** The bulk speed approaches MKM's monotonically as the first cell
  moves into the viscous sublayer, but no grid is within its tolerance of it: without wall
  damping the eddy viscosity is too large in the first cells off the wall, so the friction is too
  high. Van Driest damping or a wall model is the next step (literature C5.1, "Limits").
- **The finest grid had not settled.** At 32 cells per δ the bulk speed was still falling slowly
  at the end of the averaging, so its mean is an upper bound of that run's, and the observed
  trend is not a convergence order.
- **One realisation per grid.** The means carry the sampling scatter of a single 15 δ/u_τ
  average (about 1 % at 12 cells per δ, from the two platforms above); no confidence interval is
  stated.

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
