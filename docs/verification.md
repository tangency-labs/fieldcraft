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
