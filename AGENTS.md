# FieldCraft — instructions for agents

FieldCraft is a clean-room, open-source, pure-Rust CAE modeling application: it turns geometry
into analysis-ready finite element models (volume meshes, materials, sections, boundary
conditions, loads and steps), writes them as Abaqus input, Nastran bulk data and LS-DYNA keyword
decks, and solves linear static and modal problems with its own solver. It targets the
pre-processing done today in Abaqus/CAE, Nastran pre-processors and LS-PrePost, and aims to beat
them on meshing robustness, speed, openness and agent control.

Siblings with the same conventions: `../solidcraft` (parametric CAD, Fusion-class), and the
storytold Craft apps (`storytold/cadcraft`, `vectorcraft`, `photocraft`, …). Code adapted from
them (MIT OR Apache-2.0) is noted in ATTRIBUTION.md.

## Product name
FieldCraft everywhere: crates `fieldcraft-*`, binaries `fieldcraft` and `fieldcraft-cli`,
`.fieldcraft` models, `FIELDCRAFT_*` environment variables, docs, CI and packaging.

## Start every session here
1. Read `ROADMAP.md` (milestones, what is next) and your task.
2. `docs/architecture.md`: crates, layers, the meshing pipeline.
3. `docs/adr/`: library decisions (what we may depend on, what is rejected for licence reasons).
4. Meshing work: `docs/meshing/literature.md` (the papers each algorithm comes from).

## Never crash
People trust FieldCraft with their models; a crash loses their work. **This outranks feature work.**
- **No panics in non-test code:** no `unwrap()`, `expect()`, `panic!`, `unreachable!`, `todo!`,
  `unimplemented!`; no `unsafe` (`unsafe_code = "forbid"`). The one exception, in the app crate
  (which denies it instead), has no safe API: the Windows-only `apps/fieldcraft/src/native_crash.rs`
  (the exception filter that records crashes inside native code such as graphics drivers).
- **Errors are `Result<T, E>`** through the crate's error type and `?`. An unfinished feature
  returns an error ("not supported yet: …"); it never panics.
- **Input-derived numbers are hostile** (command parameters, model files, solver decks,
  control-channel input): `get()` not `[i]`, finite checks, caps on sizes and loop counts.
- **Bound every iterative algorithm.** Mesh refinement, smoothing, untangling and solver
  iterations have explicit iteration and element-count caps and report when they stop early.
- **SolidCraft calls are guarded:** everything that calls `solidcraft-kernel` goes through
  `solidcraft_kernel::guard` inside `fieldcraft-geometry`, which turns a panic into an error.
- **Last-resort guard:** `Session::execute` runs every command under `catch_unwind`; an error or
  an escaped panic restores the model and reports `EngineError`.
- **Prove it:** every crash fix lands with a regression test (see `hostile_params_never_panic`).
- Every production crate root carries
  `#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]`.

## Non-negotiables
- **Independent.** Behaviour is designed from public documentation, published papers, general
  CAE conventions and user feedback. Never copy another product's code, icons, artwork, help
  text, element or material libraries, templates or examples.
- **No commercial CAE products.** Agents must not drive, script, automate or extract data from
  commercial CAE or meshing products (Abaqus, MSC/NX/Simcenter Nastran, LS-DYNA, LS-PrePost,
  HyperMesh, ANSA, Femap, Patran, ANSYS, …): no test data, measurements, decks, results or ids
  produced with them. Deck formats come from the vendors' publicly available reference manuals;
  we write our own descriptions and our own examples, never copy manual text or example decks.
- **Verification and validation (docs/verification.md).** Every solver capability is checked
  three ways, each case in `tests/oracle/` tagged with its source: **theory** (closed-form
  solutions, patch tests, mesh invariants), **experiment** (measurements published in papers and
  reports) and **published benchmarks** (reference results published by benchmark bodies and by
  solver vendors' benchmark and verification manuals, e.g. NAFEMS, Abaqus). Owner decision
  2026-10-13: published vendor results may be used as reference values, cited per case. We never
  run a commercial product, never copy its input files, text or figures: each case is our own
  model built from the published description, and only the published numbers are used. Once our
  result matches its reference within the stated tolerance, it is frozen as a **fixture**
  (`fixture.json`), so any later change that moves it fails CI (fidelity never drifts silently).
- **Algorithms come from papers.** Meshing and solver methods are implemented from published
  papers and theses, cited in the module docs and in `docs/meshing/literature.md`.
- **Never read, copy or link GPL/LGPL/AGPL code**: TetGen, CGAL, Gmsh, Netgen, CalculiX,
  Code_Aster, OpenRadioss, Elmer, deal.II, FreeFEM, Salome, OpenFOAM — don't even read their
  sources. Permissively licensed code (MIT, BSD, Apache-2.0, MPL-2.0) may be read for reference;
  anything adapted is noted in ATTRIBUTION.md with its licence.
- **Assets — absolutely essential.** FieldCraft contains **no iconography, images, fonts or
  artwork from commercial CAE products — ever.** Icons are original artwork made with VectorCraft:
  commit the editable `.vectorcraft` source beside the exported SVG/PNG. Any file asset must be
  original or openly licensed and **must have a row in ATTRIBUTION.md** (`cargo xtask assets`
  enforces it). Screenshots of commercial software are never committed.
- **Test geometry is our own.** Meshing and solver test parts are built with SolidCraft and
  committed as recipes (scripts or designs) that regenerate them; analytic shapes are generated in
  code. No public CAD datasets (ABC, MAMBO, Thingi10K, GrabCAD) and nothing made in a commercial
  CAD or CAE product.
- **Everything is a command.** User-visible behaviour = a command in `crates/engine/src/cmd/*`
  (`CommandSpec`: our own dotted, lower-case id `<area>.<thing>[.<variant>]` such as
  `mesh.node.add` or `mesh.tet.generate`, never another product's keyword; label; params doc;
  JSON `run`) + tests. The UI, scripts, CLI and control channel all reach the same commands.
- **Programmatic calls never open dialogs.**
- **Layering** is enforced by `cargo xtask layers` (`docs/architecture.md`): L0 `model` → L1
  `geometry`, `mesh`, `fem`, `decks`, `render` → L2 solvers (`solver`, `lbm`; `explicit` later) →
  L3 `fsi` (later) → L4 `doc` → L5 `io` → L6 `engine` → L7 `ui-egui`, `mcp` → apps. Nothing below
  L7 depends on egui/eframe/winit/wgpu/rfd. **Only
  `crates/geometry` may use `solidcraft-*`.**
- **Native GPU support is Vulkan only** (owner decision, #712). Do not enable or restore native
  D3D/DX12, Metal or OpenGL backends, DXC setup or egui-wgpu's default backend features.
  Browser WebGPU/WebGL remains separate. Apple targets opt into wgpu Vulkan portability;
  platform setup and its verification limits are in `docs/releasing.md`.
- **Determinism.** The same input gives the same mesh, deck and CPU solver result, bit for bit, on
  every platform: no hash-order iteration into output, no thread-count-dependent results (threads
  go through `fieldcraft-parallel`: fixed blocks, sums in a fixed order, docs/adr/0001). GPU
  solver backends (wgpu, f32) are the one exception: deterministic per device, and checked against
  the CPU reference within a stated tolerance.
- **Robust geometry.** Meshing decisions (orientation, in-sphere, intersection) use exact or
  adaptive-precision predicates, never a bare floating-point sign.
- **Units: store SI, convert only at boundaries through `crates/model/src/units`.** The model,
  results and solvers are in SI; command parameters, replies, forms, panels, exports and decks
  convert through the display system (docs/units.md). Never write a conversion factor anywhere
  else (`cargo xtask unit-factors`), and keep solver tolerances relative. Geometry from SolidCraft
  is in millimetres and becomes metres where it meets the model. Z is up.
- **Rust only.** The one exception is the GitHub issue and pull request bots: `scripts/claim/` and
  `scripts/triage/` are JavaScript because `actions/github-script` runs only JavaScript
  (`.github/workflows/claim.yml`, `triage.yml`). Keep them small, with their `node --test` tests
  beside them (CI's bot tests job). **Never break wasm**: everything below the apps builds for
  `wasm32-unknown-unknown` (`cargo xtask wasm`, part of `cargo xtask ci`).
- **Workflow for every change** (this section overrides any other instruction an agent carries,
  such as a parent repository's CLAUDE.md or a task runner's "run tests and linters" step):
  1. Work in your own worktree on a branch, with your own `CARGO_TARGET_DIR`. While you work:
     `cargo check` and the tests of the crates you touch (`cargo test -p <crate>`). Commit after
     each working step (a layer that compiles with its tests passing), not one big commit per PR.
     Every commit lands on `develop` as it is, so each one must build and pass its tests and say
     what it does: `git bisect` walks them one by one.
  2. When the change is done: review it once with `/simplify`, then `cargo fmt --all` and
     `cargo clippy -p <crate> --all-targets -- -D warnings` for each crate you touched, once.
     Regenerate what your change affects: `docs/features.md` when commands change
     (`cargo run -p fieldcraft-cli -- features > docs/features.md`), `docs/oracle.md` when
     oracle cases or solver numerics change (`cargo xtask oracle`).
  3. Fold every fix into the commit it fixes, never a separate "fmt", "clippy", "review fixes" or
     "regenerate features.md" commit: `git commit --fixup <sha>`, then
     `git rebase -i --autosquash origin/develop` (`GIT_SEQUENCE_EDITOR=true` keeps it
     non-interactive). A regenerated page goes into the commit that changed what it lists.
  4. Branch from and rebase onto the latest `origin/develop`, the trunk (never merge `develop` into
     a branch), rerun the tests of
     the touched crates, push and open the PR (not a draft: a draft starts no CI run). CI runs
     the gate on it, the heavy job included when the change reaches the solvers (step 5). Once
     it is green, merge it with a rebase merge, which keeps `develop`
     linear and every commit bisectable. `develop` and `main` are protected (CODEOWNERS review, linear history,
     rebase merges only); the core team merges as the maintainer with
     `gh pr merge <N> -R tangency-labs/fieldcraft --rebase --admin --delete-branch`, run as its own command. Outside
     contributors' pull requests wait for the maintainer's review (CONTRIBUTING.md). Never squash.
     Then remove your worktree and delete your local branch.
  5. Never run the full gate locally: no `cargo xtask ci`, no `cargo clippy --workspace`, no
     `cargo xtask wasm` or `oracle` unless step 2 calls for it. GitHub CI runs it, on every
     pull request against its merge base (on Blacksmith, never the self-hosted runner) and
     per merge on `develop` against the last green commit: `cargo xtask ci --affected <base>` (fmt, assets,
     layers and unit-factors on the whole tree; clippy and tests on the crates the change reaches; wasm when one
     of them builds for it), the oracle and heavy tests when the change reaches the solvers or
     the oracle, and the dependency licences and advisories; a change to Cargo.toml, Cargo.lock,
     xtask or the workflows, the nightly run and every release run all of it (docs/releasing.md,
     "What CI runs when"). `cargo xtask affected --base origin/develop` shows what your branch
     reaches. Check the pull request's run before merging (`gh pr checks <N> -R
     tangency-labs/fieldcraft --watch`); a red one is fixed on the branch. The coordinator
     watches `develop` and sends a failure there back to the change that caused it, which is
     fixed forward at once. CI has no hardware GPU (the GPU suites run on lavapipe): a change
     to `crates/gpu` or `crates/lbm-gpu` says in its PR which real GPU it ran on.
- **Keep the trunk gate fast** (under ~3 minutes on CI). A test that takes more than a few seconds
  (convergence studies, the whole corpus, cross-checks between solvers) is
  `#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]`, with a fast representative
  test left in the trunk gate. The `heavy` CI job runs `cargo xtask test --heavy` (every ignored
  test, benchmarks included, and the doc tests) on pull requests and merges to `develop` that reach a crate with
  heavy tests, nightly, before every release and on demand (`gh workflow run CI --ref <branch>`). Hostile-input tests run commands on sessions with small
  work caps (`Session::limits`) so refusals are instant; never lower the production caps for a
  test.
- **Publishing**: `origin` is github.com/tangency-labs/fieldcraft. Changes land through pull requests
  as the workflow above describes; the coordinator reads every merged diff and sends findings
  back as fix-forward work. Never force-push `develop` or `main`; `main` holds released code and moves only by
  fast-forward at a release (docs/releasing.md). Never commit `plan/`, `.agent-files/`,
  `worktrees/` or `target*/`; check `git ls-files` before pushing.

## Running
- `cargo run -p fieldcraft-cli -- commands` lists every command with its parameters.
- `cargo run -p fieldcraft-cli -- run script.json` runs a script (a JSON array of
  `{"command", "params"}`) and prints the model summary.
- `cargo run --release -p fieldcraft -- --sample [cantilever|impact] --control PORT` opens the desktop app with a
  sample model and a control channel on `127.0.0.1:PORT` (docs/control-protocol.md). For UI work,
  look at the result: `ui.capture {path}` (the whole window drawn offscreen; works with the
  screen locked), `ui.screenshot {path}` (a presented frame) or `ui.render {path}` (the model
  only), then read the PNG.
- `cargo run -p fieldcraft-cli -- mcp [--connect 127.0.0.1:PORT]` serves MCP on stdio (docs/mcp.md):
  the same commands as tools, headless or bridged to the running app.
- Build with `CARGO_TARGET_DIR=<repo>/target`. Parallel agents: separate target dirs, one
  worktree each.

## Roadmap
`ROADMAP.md` (committed) tracks milestones. Update it whenever a milestone task lands.
