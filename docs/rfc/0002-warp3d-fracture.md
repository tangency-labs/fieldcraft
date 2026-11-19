# RFC 0002: Port WARP3D's fracture capabilities

- **Status:** proposed, open for comment
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
2. Deformation plasticity.
3. Domain-integral J, checked on CT and SENB specimens.
4. Interaction-integral K and T-stress.
5. Cohesive elements.
6. Gurson with element extinction.
7. Chaboche kinematic hardening and Norton creep.

Each step ships with its manual page and its oracle cases, and cites the papers the WARP3D routines implement.

## 6. Open questions

- Is this the order the owner wants?
- Should the iterative solvers (PCG, AMG) come before fracture?

