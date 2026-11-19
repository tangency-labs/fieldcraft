# RFC 0002: Port WARP3D's fracture capabilities

- **Status:** proposed, open for comment. **Decided:** FieldCraft ports; WARP3D does not run beside it as a sidecar solver (§6).
- **Discussion:** first raised in [#946](https://github.com/tangency-labs/fieldcraft/issues/946)
- **Scope:** port WARP3D's fracture and nonlinear-material capabilities to Rust, with its LGPL and non-free parts replaced by our own code

## 1. Summary

FieldCraft ports WARP3D's fracture and nonlinear-material capabilities to Rust. Every part of WARP3D that is LGPL, non-free or proprietary is left behind and replaced with FieldCraft's own code, written from the published papers. This RFC records the proposal, the licence audit behind it, and what each replaced part becomes.

## 2. Before any port lands

This port is an exception to two rules, and both must change first:

- **AGENTS.md, "No code from other apps".** It forbids porting code from any other CAE or solver app, whatever its licence. AGENTS.md says the rule changes in the Tangency Labs FieldCraft maintainer guide, section 1. That guide and AGENTS.md need an exception naming WARP3D's NCSA-licensed code.
- **`deny.toml` and ATTRIBUTION.md.** The NCSA licence is permissive, but `deny.toml` does not list it yet. ATTRIBUTION.md needs a row for the ported code: what it is, the NCSA notice, and the WARP3D version it came from.

Rules the port keeps:

- **The source is Fortran; the port is Rust.** It is translated into FieldCraft's own crates and conventions. It follows the "Never crash" rules, SI units through `crates/model/src/units`, and determinism through `fieldcraft-parallel`. No Fortran is linked.
- **Agents never read GPL or LGPL code.** Before any file is ported, a licence pass reads only the header of every file in the WARP3D tree. It lists the NCSA files to port and the GPL, LGPL or unclear files to replace. An agent stops at the header of any non-NCSA file and never reads past it.
- **Results are verified against standards** (docs/verification.md). Each ported capability gets its oracle cases from theory, published experiments and standard benchmarks, never by comparing against WARP3D's output.

## 3. What gets replaced

| What WARP3D uses | Licence | Replaced with | Effort |
|---|---|---|---|
| Kopp 3x3 symmetric eigen-solver (`gtmat1.f`) | LGPL | `fem::dense::symmetric_eigen` (cyclic Jacobi, Golub & Van Loan §8.5) | none, exists |
| hypre 2.7.0b (iterative solvers, AMG) | LGPL-2.1 (Apache-2.0 OR MIT only from 2.17.0, 2019-07-10) | our own PCG and AMG, from the papers | PCG about 1 week, AMG 3–4 weeks |
| MUMPS (direct solver) | CeCILL-C (LGPL-like) | our own supernodal Cholesky (Ng & Peyton 1994) | none (a speed check above 1M DOF) |
| METIS 4.0.3 (graph partitioning and ordering) | free only for education and research at non-profits and US government agencies; no redistribution without approval | our own nested-dissection ordering (`crates/solver/src/ordering.rs`) | none |
| Intel MKL (Pardiso) | proprietary | our own Cholesky | none |
| NASA VSS (sparse solver) | needs NASA's written approval; export-restricted | our own Cholesky | none |
| any other non-NCSA file found by the licence pass | — | written from the paper it implements | per file |

METIS 5.1.0 is Apache-2.0, but some GKlib files are LGPL-2.1+, so it is not linked either. Multilevel partitioning is published (Karypis & Kumar 1998) and can be written in Rust.

## 4. Graph ordering and parallel factorisation (what we already have)

- `crates/solver/src/ordering.rs` orders by nested dissection.
  - Separators come from coordinate bisection, or from the graph alone (George & Liu 1978), with deterministic tie-breaks.
  - Fill and time per ordering are measured.
- `crates/solver/src/cholesky.rs` builds the elimination tree (Liu) and fundamental supernodes (Liu, Ng & Peyton 1993), and factors left-looking (Ng & Peyton 1994).
  - Disjoint subtrees of the supernodal tree factor concurrently.
  - Large supernodes split into row tiles.
  - All of it runs through `fieldcraft-parallel`, so results are bit-identical whatever the thread count.
- Open contribution: multilevel graph separators, to cut fill and balance the subtrees.
  - Coarsen, partition, then refine with Kernighan–Lin or Fiduccia–Mattheyses.
  - Later, partitioning for distributed or iterative solves.
  - Both must stay deterministic and come from papers.

## 5. Port order

1. Licence pass over the WARP3D tree (headers only), plus the AGENTS.md, maintainer-guide, `deny.toml` and ATTRIBUTION.md changes above.
2. The crack definition, `model.crack.create`, with its crack-front ordering (§6.2).
3. Domain-integral J on linear-elastic solids, checked on CT and SENB specimens and the semi-elliptical surface crack (§6.3).
4. Deformation plasticity, and J with it.
5. Interaction-integral K and T-stress.
6. Cohesive elements.
7. Gurson with element extinction.
8. Chaboche kinematic hardening and Norton creep.

Each step ships with its manual page and its oracle cases, and cites the papers the WARP3D routines implement.

## 6. The deck and sidecar proposal, and what the port takes from it

A contribution proposed running WARP3D beside FieldCraft: FieldCraft writes a WARP3D input deck, runs the WARP3D binary (`solve.static` with `backend: "warp3d"`), and reads its result files back (`results.read`). FieldCraft does not take that route.

### 6.1 Why the port, not a sidecar

- **One solver path.** Every result FieldCraft shows comes from its own solver, verified against standards and frozen as oracle fixtures. A sidecar's results depend on an outside binary and its version, which our CI cannot pin or check.
- **Determinism and the browser.** Our solvers give the same result bit for bit on every platform and build for wasm. Starting a process does neither.
- **Nothing to install.** Users get fracture mechanics in the app they downloaded, on macOS, Linux, Windows and in the browser.
- **`backend` already has a meaning.** On solve commands it chooses CPU, threads or GPU for our own solvers, not another program.

So there is no `backend: "warp3d"`, no runner that starts WARP3D, and no reader for WARP3D's result files.

### 6.2 What carries over

The deck work solved problems the port has too, and these parts become FieldCraft's own:

- **The crack is a model definition**, not a parameter of one command: `model.crack.create` with
  - the front: a node set, ordered or not, which FieldCraft chains into an open path or a closed loop, taking the hex20 mid-side nodes along;
  - the crack-plane normal, or inferred from the symmetry conditions on the front;
  - `symmetric` for a half model (the energy on both crack faces);
  - the number of rings of the integration domains;
  - the weighting along the front: `pointwise` (a value at every station, tent functions inside and ramps at open ends) or `front_average` (one value for the whole front).

  The domain-integral step, the UI, the model check and any deck export all read this one definition.
- **The domain integral and its weighting** come from the papers, never from WARP3D's routines: C. F. Shih, B. Moran, T. Nakamura, *Energy release rate along a three-dimensional crack front in a thermally stressed body*, Int. J. Fracture 30:79–102, 1986; B. Moran, C. F. Shih, *A general treatment of crack tip contour integrals*, Int. J. Fracture 35:295–310, 1987. For the interaction integral and T-stress: T. Nakamura, D. M. Parks, *Determination of elastic T-stress along three-dimensional crack fronts using an interaction integral*, Int. J. Solids Struct. 29(13):1597–1611, 1992; M. Gosz, J. Dolbow, B. Moran, *Domain integral formulation for stress intensity factor computation along curved three-dimensional interface cracks*, Int. J. Solids Struct. 35(15):1763–1783, 1998.
- **Crack-face pressure** is an ordinary surface load on the crack faces. Inside FieldCraft no face renumbering is needed.
- **Tet fronts** are refused with a clear message until the integral supports them.
- **Results**: J, K_I, K_II, K_III and T at every station, per ring, in the session's results; `solve.probe` reads them, and a results panel plots them along the front (as the path plot does).

### 6.3 Verification

The semi-elliptical surface crack in a plate under crack-face pressure is a good first 3D case. It becomes an oracle case on our own model, built in SolidCraft from the published description:

- the reference values are the published ones: X. Wang, R. Bell, Eng. Fract. Mech. 71(9):1477–1496, 2004 (its tables for the surface crack), and I. S. Raju, J. C. Newman, *Stress-intensity factors for a wide range of semi-elliptical surface cracks in finite-thickness plates*, Eng. Fract. Mech. 11(4):817–829, 1979;
- WARP3D's verification inputs and meshes (`verification/test87` and the like) are not copied, and the numbers WARP3D's own scripts check (`driver.inp`) are not references: they are another solver's output.

### 6.4 A WARP3D deck export

Writing a WARP3D input deck, as FieldCraft writes Abaqus, Nastran and LS-DYNA decks, stays possible as an export for people who run WARP3D themselves: `deck.write` with `format: "warp3d"`, in `crates/decks`, written from WARP3D's public user manual, reading the crack from `model.crack.create`. It is not part of the port and FieldCraft does not run what it writes. Whether to take it is open (§7).

## 7. Open questions

- Is this the order the owner wants?
- Should the iterative solvers (PCG, AMG) come before fracture?
- Take the WARP3D deck export (§6.4) now, later, or not at all?

