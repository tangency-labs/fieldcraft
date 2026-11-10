# ADR 0001: Libraries and licences

Status: accepted (2026-10-13)

## Context

FieldCraft is MIT OR Apache-2.0. Everything it links must allow that, and the meshing and solver
algorithms, where the obvious open-source implementations are GPL-family, must come from papers.

## Decision

**Geometry: SolidCraft**, as git dependencies pinned by rev (`solidcraft-kernel`,
`solidcraft-geom`, `solidcraft-doc`, `solidcraft-io`; MIT OR Apache-2.0, built on the Apache-2.0
truck crates until v0.30 and on solidcore, an Apache-2.0 Rust port of NVIDIA's SMLib, from
v0.40, with its `solidcore-record` geometry records from the copy SolidCraft vendors). SolidCraft
is the primary geometry source: FieldCraft opens SolidCraft designs and
meshes from SolidCraft's tessellation, and reads face-tagged meshes SolidCraft exports. Other CAD
exchange files (STEP, which SolidCraft's kernel already reads) can be added as secondary inputs
through the same crate. Only `fieldcraft-geometry` may depend on SolidCraft.

**Desktop UI: SolidCraft's UI kit** (`solidcraft-ui-kit` and the `solidcraft-theme` token contract
it draws from; same pinned rev, under the SolidCraft Brand and Interface Licence, see below).
*Decided 2026-11-10 (owner):* the desktop UI
draws its windows, confirmations, menus, form rows, fields and buttons with the kit instead of
FieldCraft's adapted copies of SolidCraft's widgets, so the two apps share one implementation and
FieldCraft follows SolidCraft's egui line. FieldCraft installs its own tokens (its teal accent on
SolidCraft's chrome) and keeps what is its own: the viewport, mesh, element, results and setup
colours. `fieldcraft-ui-egui` is the only crate that may use them (`cargo xtask layers`); it uses
nothing else of SolidCraft's. The kit builds for `wasm32-unknown-unknown`, so the browser app keeps
the same UI. craft-ui (storytold/craft-libs), used before for tabs, icon buttons and widget states,
is no longer a dependency.

**SolidCraft's brand-licensed crates: used by written permission.** From SolidCraft 2874e78 on,
`solidcraft-ui-kit`, `solidcraft-theme` and `solidcraft-look` are under the SolidCraft Brand and
Interface Licence (each crate's `LICENSE`), not MIT or Apache-2.0: they are **not permissively
licensed**, and §3a and §3b of that licence require the maintainers' written permission to link
them from, or distribute them with, software that is not SolidCraft. *Permission requested
2026-11-10* by FieldCraft's owner (GitHub `jyviko`) for FieldCraft to link these three crates and
to distribute them as part of FieldCraft. *Approval* by SolidCraft's maintainer bherbruck is
recorded by their approving review of the pull request that adds this section,
https://github.com/tangency-labs/fieldcraft/pull/872; FieldCraft does not ship these crates without
it. The permission covers those three crates only: any other SolidCraft crate under the brand
licence needs its own permission, recorded here, before FieldCraft depends on it. `deny.toml` accepts the licence for exactly these three crates (a
`LicenseRef-SolidCraft-Brand` clarification and exception each), never in its global allow list.
FieldCraft does not use `solidcraft-look` yet: whether it draws SolidCraft's icons and themes or its
own is still the owner's decision.

**Meshing: our own**, implemented from published papers (`docs/meshing/literature.md`).
Robust predicates follow Shewchuk's adaptive-precision arithmetic (the paper; his reference
`predicates.c` is public domain). A permissive crate may be adopted instead after review; record
it here. *Decided 2026-10-13:* our own implementation, `fieldcraft_mesh::predicates`, written from
Shewchuk 1997, Devillers–Teillaud 2011 and Attene 2020 (no `unsafe`, deterministic, wasm-safe, a
validated coordinate domain instead of Shewchuk's overflow/underflow caveat, and indirect
predicates that `robust` lacks). `robust` (georust, MIT OR Apache-2.0) and `num-bigint` /
`num-rational` / `num-traits` (MIT OR Apache-2.0) are dev-dependencies only: differential test
oracles, never linked into FieldCraft.

**Linear algebra: our own.** *Decided 2026-10-13 (implicit solver task):* `fieldcraft-solver`
has its own sparse symmetric storage, nested-dissection ordering (George 1973, by coordinate
bisection), up-looking sparse Cholesky (Davis 2006, elimination trees per Liu 1990) and
shift-invert block Lanczos with full reorthogonalization (Ericsson–Ruhe 1980, Grimes–Lewis–Simon
1994) on a Householder + implicit-QR dense eigensolver (Golub–Van Loan). `faer` (MIT) was the
alternative and is faster on large problems, but it selects SIMD kernels by CPU feature at run
time (AVX2, AVX-512, NEON), which changes floating-point reduction order between machines and
breaks the bit-for-bit determinism rule; it is also a large dependency tree for the browser
build. Our code is sequential with a fixed operation order, `unsafe`-free and wasm-safe, and its
factor size is known (and capped) before any value is computed. Cost if this is wrong: large
models (beyond a few hundred thousand equations) factor slower than they could; the remedy is a
supernodal factorization of our own, or `faer` behind the same `Factor` interface with a stated
tolerance instead of bit-identity, as the GPU exception already does. *Updated 2026-10-13 (CPU
parallelism task):* the factorization is now our own supernodal left-looking Cholesky (Liu–Ng–Peyton
1993 fundamental supernodes on a postordered elimination tree, Ng–Peyton 1993, Rothberg–Gupta
1994), multithreaded with a fixed operation order (below). Its rounding differs from the old
up-looking factor in the last digits; every oracle fixture stays within its band.

**Threads: rayon** (MIT OR Apache-2.0), in `fieldcraft-parallel` only (an L0 crate with no
workspace dependencies; *decided 2026-10-13*). The solvers keep the bit-for-bit rule on any
number of threads by construction: work is cut into blocks whose sizes are constants of the
algorithm, never derived from the thread count; each floating-point value is computed by one task
in an order fixed by the input; results come back in block order and every sum over blocks or
elements (nodal forces, assembly, wall forces, energies) runs on one thread in the order the
sequential algorithm used. No float atomics, no `reduce` whose tree depends on work stealing:
rayon decides only when a block runs. `std::thread::scope` was the alternative and needs no
dependency, but it spawns threads per parallel region; the explicit solver enters several regions
per time increment (10⁵–10⁷ increments), so a persistent pool is needed, and writing our own would
re-implement rayon's. rayon was already in the shipped dependency tree (THIRD-PARTY-LICENSES.txt
unchanged). In the browser build (`wasm32`) the crate does not depend on rayon and runs the same
blocks one after another, with the same results. The thread count is a session and `solve.*`
parameter (`threads`, 1–256) whose default is `FIELDCRAFT_THREADS` or the available cores. Cost if
this is wrong: rayon's pool is a process-wide resource and nested pools add latency; replacing it
touches only `fieldcraft-parallel`.

**Machine memory: sysinfo** (MIT), in `fieldcraft-engine` only, desktop and CLI builds (not the
browser), with only its `system` feature, for the physical memory (*decided 2026-10-21, #650*). It
bounds how far a fluid solve's memory caps may be raised: half the physical memory, or on Linux
half the root cgroup's limit when that is lower (a container's). std has no such query; the
alternatives were our own per-platform code (`/proc/meminfo` on Linux is plain Rust, but macOS and
Windows need `sysctl` and `GlobalMemoryStatusEx` through FFI, which `unsafe_code = "forbid"`
rules out) or a crate built on C code. Cost if this is wrong: the bound is a little off on some
platform, and the remedy stays inside `fieldcraft_engine::machine`.

**Allowed:** MIT, BSD-2/3, Apache-2.0, ISC, Zlib, Unicode, MPL-2.0 (as an unmodified dependency),
CC0 and public domain.

**Rejected — never read, copy or link:**

| Project | Licence | Would have offered |
|---|---|---|
| TetGen | AGPL-3.0 | constrained Delaunay tetrahedralization |
| CGAL | GPL-3.0 / LGPL-3.0 | Delaunay, mesh generation |
| Gmsh | GPL-2.0+ | tet and hex meshing |
| Netgen | LGPL-2.1 | tet meshing |
| CalculiX | GPL-2.0 | implicit solver, Abaqus-like input |
| Code_Aster | GPL-3.0 | implicit solver |
| OpenRadioss | AGPL-3.0 | explicit solver |
| Elmer | GPL-2.0 / LGPL-2.1 | multiphysics solver |
| FreeFEM, Salome | LGPL | FE libraries, pre-processing |
| Palabos, OpenLB, waLBerla, lbmpy, VirtualFluids, TCLB | AGPL / GPL | lattice Boltzmann (implement from the papers) |
| OpenFOAM, SU2, 4C, DualSPHysics, preCICE | GPL / LGPL | CFD, FSI coupling, SPH |
| FluidX3D | non-commercial only | GPU lattice Boltzmann |
| LULESH | no licence | explicit hydrodynamics proxy |
| OpenCASCADE | LGPL-2.1 | B-rep kernel |
| Attene's CDT, Indirect_Predicates, VolumeMesher | GPL / LGPL | robust CDT, indirect predicates (implement from the papers) |
| TetWild (original), Evocube, robustPolycube, ultimaille, invertible-maps | GPL-3.0 / AGPL-3.0 | tet meshing, polycube, injective maps |
| libHexEx, libHexHex, QGP3D, MC3D, CoMISo, Graphite | GPL | integer-grid-map hex extraction, quantization |
| Shewchuk's Triangle | non-commercial only | 2D meshing |
| geogram `third_party/tetgen`, `third_party/triangle` and their wrappers | AGPL / non-commercial | (the rest of geogram is BSD-3 and readable) |
| libigl `copyleft/` modules | GPL | (libigl core is MPL-2.0 and readable) |
| OpenVolumeMesh, CGoGN | believed LGPL, unverified | treat as forbidden until checked |

deal.II is now dual-licensed Apache-2.0 (with LLVM exception) or LGPL, so it is readable under its Apache-2.0 option. Readable solver references (re-implement from the papers): STLBM (MIT), XLB, lettuce, pylbm, Laghos, MFEM, IBAMR, the Kratos core, MUI, FEBio, SPHinXsys, SPlisHSPlasH, fenris; full table in `docs/solvers/literature.md`.

**Repositories without a licence are all rights reserved: never read them** (for example Gao
2019's octree hex code, LoopyCuts). Readable references for the meshing work (re-implement from
the paper, note anything adapted in ATTRIBUTION.md): Shewchuk's `predicates.c` (public domain),
`georust/robust` (Apache-2.0), geogram (BSD-3, exclusions above), fTetWild (MPL-2.0), CinoLib,
Gen-Adapt-Ref, HexaLab, robust_hex_dominant_meshing (MIT), quartet (ISC). The full table with
sources is in `docs/meshing/literature.md`.

**Test geometry** follows the same rule as test data: nothing produced with a commercial CAD or
CAE product. That excludes ABC (Onshape), MAMBO (files from commercial translators) and the HexMe
models that come from ABC, GrabCAD or Siemens NX; HexMe's Gmsh-made models are allowed (CC-BY-4.0).
Thingi10K is not used either (owner decision 2026-10-13: authoring tools unrecorded). **The test
corpus is our own:** parts built with SolidCraft (scripts or designs authored in SolidCraft,
committed as their recipes so the geometry is regenerated, not stored) and analytic shapes.
Recipes are SolidCraft command scripts run through `solidcraft-engine` (same pinned rev), a
dependency of `fieldcraft-geometry` only behind its non-default `scripts` feature, which the test
tooling (`fieldcraft-testkit`) enables; the shipped library opens designs without it.

**Commercial products** (Abaqus, Nastran, LS-DYNA, LS-PrePost, HyperMesh, ANSA, Femap, Patran,
ANSYS) are never run or scripted, and their input files, text and figures are never copied. Deck
formats come from the vendors' public reference manuals, described in our own words with our own
examples. Published numerical results (benchmark and verification manuals) may serve as cited
reference values for our own models of the same problems (owner decision 2026-10-13;
docs/verification.md).

## Consequences

- The tet and hex meshers are significant original work; the literature survey decides the
  order and the methods.
- SolidCraft and FieldCraft move together: bumping the pinned rev bumps `solidcore-record` (from
  SolidCraft's repository) and the repeated `[patch.crates-io]` entries in the same commit. The
  pin is a SolidCraft release tag's commit or, when FieldCraft needs a change not yet released, a
  commit on SolidCraft's main branch,
  and a scheduled check opens designs saved by SolidCraft's latest release, so a release that
  saves a design format the pin can't read is caught before users meet it
  (docs/architecture.md, "SolidCraft pin").
