# Contributing to FieldCraft

FieldCraft is one open-source modeler for structures, heat and fluids, written in pure Rust: open a
parametric SolidCraft part, mesh it (tet or all-hex), define the analysis, solve it with our own
implicit, explicit, thermal or lattice Boltzmann solvers, and look at verified results — from the
desktop app, the command line, or an AI agent over MCP.

Help is welcome, whether that's a one-line fix, a verification benchmark, a file format or owning a
whole solver. This page covers how to build it, the few rules that aren't negotiable, where help is
needed most, and how to claim work. The rules every change follows, by a person or an agent, are in
[AGENTS.md](AGENTS.md).

## Get it running (about 10 minutes on a fast machine)

You need stable Rust 1.95 or newer. Native graphics and GPU solvers require Vulkan and a
compatible driver; no DXC setup is needed. See [platform requirements](docs/releasing.md),
including the separate Vulkan portability installation required on macOS.

```sh
git clone https://github.com/tangency-labs/fieldcraft && cd fieldcraft
cargo run --release -p fieldcraft -- --sample cantilever   # the desktop app on a sample model
cargo run -p fieldcraft-cli -- commands                    # every command and its parameters
cargo test -p fieldcraft-solver                            # the tests of one crate
cargo xtask oracle --quick                                 # the verification cases (docs/oracle.md)
```

`develop` is the trunk and `main` holds released code. The layers and crates are in
[docs/architecture.md](docs/architecture.md); every user-visible capability is a command in
`crates/engine/src/cmd/` ([docs/features.md](docs/features.md) lists them).

## Where help is wanted most

See the [`help wanted`](https://github.com/tangency-labs/fieldcraft/labels/help%20wanted) and
[`good first issue`](https://github.com/tangency-labs/fieldcraft/labels/good%20first%20issue) labels, and
the `area:` labels for your field. In order of how much they would help users:

1. **Reproducible meshes** ([#115](https://github.com/tangency-labs/fieldcraft/issues/115)). The same part
   must give the same tet mesh on every platform and every run, as the solvers already do.
2. **All-hex meshing of thin and curved parts** ([#85](https://github.com/tangency-labs/fieldcraft/issues/85)).
   The octree mesher over-refines thin walls; concave transition templates and an up-front depth
   are next. Hex meshing is the centre of gravity of the project.
3. **Getting results in and out**: ParaView export ([#90](https://github.com/tangency-labs/fieldcraft/issues/90)),
   Gmsh import ([#91](https://github.com/tangency-labs/fieldcraft/issues/91)), Abaqus deck reading
   ([#98](https://github.com/tangency-labs/fieldcraft/issues/98)), CSV and STL
   ([#92](https://github.com/tangency-labs/fieldcraft/issues/92), [#93](https://github.com/tangency-labs/fieldcraft/issues/93)).
4. **Verification benchmarks** ([#97](https://github.com/tangency-labs/fieldcraft/issues/97),
   [#106](https://github.com/tangency-labs/fieldcraft/issues/106)). An engineer with a published benchmark
   can add it as an oracle case without touching solver code.
5. **Nonlinear elements and materials**: incompatible modes at finite strain
   ([#82](https://github.com/tangency-labs/fieldcraft/issues/82)), kinematic hardening
   ([#81](https://github.com/tangency-labs/fieldcraft/issues/81)), and the material library
   ([#99](https://github.com/tangency-labs/fieldcraft/issues/99)).
6. **The app**: typed command forms ([#112](https://github.com/tangency-labs/fieldcraft/issues/112)), XY and
   path plots ([#94](https://github.com/tangency-labs/fieldcraft/issues/94), [#95](https://github.com/tangency-labs/fieldcraft/issues/95)).

Bigger areas on the [roadmap](ROADMAP.md) — contact, dynamics, composites, coupled fields, wave
propagation, turbulence — are open to anyone who wants to own one: start a thread in
[Discussions](https://github.com/tangency-labs/fieldcraft/discussions).

## Before you start: claim an issue

Every change starts from an issue, and only the maintainer merges.

1. Pick an open issue labelled **help wanted** or **good first issue** that is not labelled
   **claimed**. If what you want to do has no issue, open one first and wait for a reply.
2. Comment `/claim` on it. The bot assigns it to you, labels it **claimed** and replies; then it is
   yours. It refuses an issue someone else holds (naming them and when their claim lapses), an
   issue an open pull request already covers, and an issue without **help wanted** or **good
   first issue**: on those, say what you would like to do and the maintainer replies. Work on
   unassigned issues may not be accepted.
3. Open a draft pull request early, linked with `Closes #N`. A pull request that references an
   open-for-contributors issue nobody holds (`Closes #N` in the body, `(#N)` in the title) claims
   it the same way; if someone holds it already, or an older open pull request references it,
   the bot says so on your pull request instead.
4. An issue with no visible progress for 14 days goes back to the pool: a daily sweep frees a claim
   with no linked pull request 14 days after it was made, and one with a pull request after 14
   days without pushes. `/unclaim` releases it early (you or a maintainer).

Issues labelled **in progress** are being worked on by the core team; please don't start on them.
The bot never lets them be claimed, and never frees an issue a maintainer holds.

## Labels and triage

The triage bot labels every new issue and pull request (`.github/workflows/triage.yml`). It only
adds labels; correct them freely and it won't put them back.

- **Pull requests** get `area:*` from the paths they change: `crates/mesh` → `area:mesh`,
  `crates/lbm*` → `area:cfd`, `crates/solver`, `crates/fem`, `crates/explicit` → `area:solver`,
  `crates/ui-egui` and `apps/` → `area:ui`, `crates/io`, `crates/decks` → `area:io`, `docs/` and
  top-level `*.md` → `area:docs`, `.github/`, `xtask/`, `packaging/`, `scripts/`, `tools/` →
  `area:infra`.
- **Issues** get `bug` or `enhancement` from the form they were filed with (or a crash or
  "Feature request" title), and an `area:*` from words in the title or the command an in-app
  report names. Until it has a priority, an issue carries `needs-triage`; adding a `priority:`
  label removes it.
- **Maintainers** comment `/triage area:mesh high S` to sort an issue or pull request: any
  existing label (case, and `-` or `_` for a space, don't matter: `in-progress`), `high` for
  `priority: high`, `S`/`M`/`L` for a size (which replaces the old one), `-label` removes one.
  Any `/triage` with a change clears `needs-triage`.
- **Short bug reports without steps to reproduce** get `needs-info` and one comment asking for
  the steps (and the version, when it's missing). The reporter's reply, or an edit that adds the
  steps, removes the label; an issue that stays in `needs-info` for 14 days is closed, and can be
  reopened any time. Feature requests, reports filed from the app (**Help › Report a Bug…**),
  and maintainers' issues are never asked.
- `needs-triage` is the queue: [open issues with it](https://github.com/tangency-labs/fieldcraft/issues?q=is%3Aopen+label%3Aneeds-triage).

## Reporting a problem or an idea

Use the templates under **New issue**. Show it: for anything you can see — the view, a mesh, a
panel, a plot — add a screenshot or a short screen recording (drag it into the issue) and mark
the problem on it. For a mesh, a section (`fieldcraft-cli snapshot … --section`) shows what the
outside hides. A picture is often the fastest way for someone else to see what you see.

## Making the change

1. Branch from `develop` (the trunk; `main` holds released code only) and open the pull request
   against `develop`. Keep it to one topic.
2. Before pushing: `cargo fmt --all`, and `cargo clippy -p <crate> --all-targets -- -D warnings`
   and `cargo test -p <crate>` for the crates you touched. CI runs the gate on your pull
   request once the maintainer approves the run (a draft starts none): the crates your change
   reaches, and the heavy tests and the oracle when it reaches the solvers.
3. User-visible behaviour is a command (`crates/engine/src/cmd/`) with tests, including hostile
   parameters.
4. Meshing and solver algorithms cite the paper they come from in the module docs.
5. A change you can see (the view, meshes, panels, results, icons) shows a picture before and after
   in the pull request.
6. CI has no hardware GPU: it runs the GPU suites on lavapipe (software Vulkan). A change to
   `crates/gpu` or `crates/lbm-gpu` says in the pull request which real GPU it was run on
   (Vulkan, for example an NVIDIA or AMD card on Linux or Windows), with the tests it ran.
7. Never contribute code, text or data from GPL-family projects or commercial CAE products
   (see [docs/adr/0001-libraries.md](docs/adr/0001-libraries.md)).

By contributing you agree your work is licensed under MIT OR Apache-2.0.

`develop` and `main` are protected: changes land on `develop` only through pull requests that the
maintainer has reviewed,
as rebase merges (keep each commit building and passing its tests). Be kind: see the
[code of conduct](CODE_OF_CONDUCT.md).
