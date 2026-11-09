# Linear solvers

The sparse direct solvers behind every implicit analysis: a fill-reducing nested dissection
ordering of the mesh nodes, a supernodal left-looking Cholesky factorization A = L Lᵀ for symmetric
positive definite systems, and a symmetric indefinite factorization P A Pᵀ = L D Lᵀ with bounded
Bunch–Kaufman pivoting on the same symbolic structure, which also reports the matrix's inertia.
Both factor in parallel with a fixed operation order, so the factor is the same bit for bit on any
number of threads. The factor's size is computed and checked against a cap before any of its
memory is committed.

| | |
|---|---|
| **Commands** | the `max_factor_entries`, `max_dofs`, `max_iterations` and `threads` parameters of `solve.static`, `solve.heat`, `solve.frequency`, `solve.buckle`, `solve.nonlinear` |
| **Step procedure** | `static`, `coupled_temperature_displacement`, `heat_transfer`, `frequency`, `buckle`, `nonlinear_static` |
| **Solver** | `crates/solver` (CPU, f64): `ordering`, `cholesky`, `ldlt`, `constraints`, `sparse` |
| **Method** | nested dissection (George 1973; Simon 1991; George & Liu 1978); elimination tree and supernodes (Liu 1990; Liu, Ng & Peyton 1993); left-looking supernodal Cholesky (Ng & Peyton 1993; Rothberg & Gupta 1994); bounded Bunch–Kaufman pivoting (Ashcraft, Grimes & Lewis 1998) within supernodes (Schenk & Gärtner 2006) |
| **Verified by** | `cholesky::tests`, `ldlt::tests`, `ordering::tests`, `constraints::tests`, `sparse::tests`; every solver test and oracle case runs on them |

## Synopsis

```json
{"command": "solve.static", "params": {"step": "Load", "max_factor_entries": 500000000, "threads": 16}}
```

This raises the factor cap for this solve only (later solves keep the default) and runs on 16
threads with the same result as on one.

## Parameters

| Field | Value | Default | Meaning |
|---|---|---|---|
| `max_factor_entries` | 1 to 2 000 000 000 | 400 000 000 | cap on the entries of the factor L (8 bytes each); may be raised deliberately for one solve |
| `max_dofs` | 1 to 6 000 000 | 6 000 000 | cap on the model's unknowns; may only be lowered |
| `max_iterations` | 1 to 3000 | 3000 | frequency and buckling: cap on the Lanczos vectors; may only be lowered |
| `threads` | 1 to 256 | the session's (`FIELDCRAFT_THREADS` or the cores) | changes the time, never the result |

Internal constants (`crates/solver/src/cholesky.rs`, `ldlt.rs`, `ordering.rs`):

| Constant | Value | Role |
|---|---|---|
| `PIVOT_TOLERANCE` | 10⁻¹⁰ | relative pivot size below which a matrix is singular to working precision |
| `PANEL` | 64 | columns per panel of a supernode's dense factorization |
| `TILE` | 64 | rows per tile of the dense kernels |
| `PARALLEL_WORK` | 2²¹ multiply–adds | dense work from which a kernel's tiles run in parallel |
| `SUBTREE_SHARE` | 256 | subtrees are handed to threads in pieces of about 1/256 of the total work |
| `ALPHA` | (1 + √17)/8 | Bunch–Kaufman growth bound |
| `ROOK_STEPS` | 64 | most steps of one rook search |
| `LEAF` | 16 | parts this small are ordered as they come |
| `PERIPHERAL_SEARCHES` | 8 | searches for a pseudo-peripheral node per part |

## Theory

### Assembly and pattern

Matrices are stored in compressed rows with both triangles (`sparse::Csr`). The pattern is built
once from the element connectivity on the node graph: a node's row holds every node it shares an
element with, and each node pair is a dense block of their unknowns (a node's unknowns are
contiguous in the `DofLayout`), so stiffness and mass share it. Sums run in a fixed order
(elements by label, then the element's own node order).

### Ordering

Nested dissection (George 1973) orders a separator after the two parts it splits, so the factor
of a 3D mesh of n nodes fills O(n^{4/3}) instead of the O(n^{5/3}) of a band or envelope ordering.
Two separators are tried and the order whose factor has fewer entries is kept
(`fill_reducing_order`):

- **coordinate bisection** (`nested_dissection`): split at the median of the longest axis of the
  part's bounding box; the nodes of the smaller side that touch the other side are the separator
  (Simon 1991, with vertex separators);
- **graph dissection** (`graph_dissection`): the separator is a middle level of the breadth-first
  level structure rooted at a pseudo-peripheral node (George & Liu 1978; Gibbs, Poole &
  Stockmeyer 1976), which cuts a thin curved wall (a shell meshed with solids, a bellows, a
  revolved part) across its surface where a plane would cut it many times.

The symbolic count that compares them costs O(|L|) time and O(n + |A|) memory. Ties are broken by
node index, so the ordering is deterministic. The free unknowns are numbered node by node in this
order; reference point equations follow the mesh nodes.

![Nested dissection: separators ordered after the parts they split](../../media/manual/linear-solvers-1.svg)

### Sparse Cholesky

`fieldcraft_solver::cholesky` factors a symmetric positive definite A = L Lᵀ.

**Symbolic.** The elimination tree (Liu 1990) and the column counts of L come from the row
subtrees (Davis 2006, §4.1–4.4), so the factor's size is known and checked against
`max_factor_entries` before any memory is committed. The columns are renumbered in a postorder of
the tree (an equivalent ordering with the same fill), so every subtree is a contiguous range of
columns, and consecutive columns with nested structures merge into fundamental supernodes (Liu,
Ng & Peyton 1993): a supernode's columns share one row structure and are stored as one dense
lower trapezoid, so the factor holds exactly the entries of L.

**Numeric** (left-looking supernodal; Ng & Peyton 1993; Rothberg & Gupta 1994). A supernode
starts as its columns of A, receives the update L_K L_Kᵀ of every descendant supernode K whose
rows reach it (K in increasing order), and is factored densely by panels of 64 columns: a panel is
updated with the supernode's columns before it, its diagonal block is factored column by column,
and the rows below are solved against that block.

**Threads and determinism.** Each entry of L is computed by one task, by a sequence of operations
fixed by the matrix and the constants `TILE` and `PANEL`: sums over a supernode's columns run in
column order from zero, and updates are subtracted in supernode order. Disjoint subtrees of the
supernodal tree are factored concurrently (a supernode reads only its descendants), and the dense
work of a large supernode is split into row tiles; neither changes any entry's sequence of
operations, and multiplications and additions are never fused, so the factor is the same bit for
bit on any number of threads and on every platform.

**Singularity.** A pivot at or below 10⁻¹⁰ times the matrix's own diagonal entry means the matrix
is singular to working precision (an unconstrained rigid-body motion, a disconnected part, a
zero-stiffness direction): the factorization stops and names the column, the first such column
in the postorder whichever thread meets it.

![A supernodal factor and its tree](../../media/manual/linear-solvers-2.svg)

### Symmetric indefinite LDLᵀ

`fieldcraft_solver::ldlt` factors P A Pᵀ = L D Lᵀ, L unit lower triangular and D block diagonal
with 1 × 1 and 2 × 2 pivots, on the Cholesky factor's symbolic factorization: the same ordering,
elimination tree, supernodes, fill, factor-size cap and parallel subtree tasks.

**Pivoting.** The bounded Bunch–Kaufman (rook) search (Bunch & Kaufman 1977; Ashcraft, Grimes &
Lewis 1998): with α = (1 + √17)/8 and ω_i the largest off-diagonal magnitude of column i, column k
is a 1 × 1 pivot when |a_kk| ≥ α ω_k; otherwise the search walks from column i to the row r of its
largest entry: r is a 1 × 1 pivot when |a_rr| ≥ α ω_r, (i, r) a 2 × 2 pivot when ω_r = ω_i, else
the search continues from r (at most 64 steps; each strictly raises ω, and at the cap the last two
columns are the 2 × 2 pivot). Every entry of L is then bounded by 1/(1 − α) ≈ 2.78 for 1 × 1
pivots. A 2 × 2 pivot has |a_ii a_rr| < α² a_ir², a negative determinant: one positive and one
negative eigenvalue.

**Pivots within supernodes** (Schenk & Gärtner 2006). The columns of a supernode are fully summed
once its descendants' updates are in, so any of them can be the next pivot: a symmetric
interchange inside the supernode leaves every row structure unchanged, so nothing is delayed and
the static symbolic factorization holds (multifrontal codes delay a failed column into the parent
instead, which changes the structure as it goes; Duff & Reid 1983; Liu 1992). The search looks at
the supernode's own rows only. Inside a panel each column is brought up to date with the panel's
columns before it when the search looks at it (the partial factorization of LAPACK's `dsytrf`,
Anderson et al. 1999); after the panel the rest of the supernode receives its update by the tiled,
parallel kernel of the Cholesky factor. Every pivot decision is made by one task, ties broken by
the lowest index, so the factor, its pivots and its inertia are the same on any number of threads.

**Inertia** (Sylvester's law): the counts of positive and negative eigenvalues of A are those of
D, a 1 × 1 pivot by its sign, a 2 × 2 one by its eigenvalues. A factor of K − σM counts the
eigenvalues below σ; a tangent's negative pivots count its unstable directions.

**Refinement.** `Symmetric::solve_refined` solves A x = b, then refines with the residual against
A while it shrinks, at most a given number of times (iterative refinement in working precision;
Higham 2002, §12.2). Cholesky solutions are returned as solved.

**Choice of factorization** (`Symmetric::new(a, cap, kind)`): `Definite` is Cholesky, refusing a
matrix that is not positive definite (the default wherever a matrix is, so linear results are
unchanged bit for bit); `Indefinite` is LDLᵀ; `Auto` is Cholesky and, when the matrix turns out
not to be positive definite, LDLᵀ on the same symbolic factorization.

| Analysis | Matrix | Factorization |
|---|---|---|
| Linear static, coupled step | K_ff | Cholesky |
| Heat transfer | K_ff (steady), C + θΔt K (transient) | Cholesky |
| Frequency | K − σM | Cholesky |
| Buckling without a shift | K | Cholesky (`Definite`) |
| Buckling with a shift | K + σ K_σ | LDLᵀ (`Indefinite`), refined twice per solve |
| Nonlinear static | the tangent | `Auto`: Cholesky, LDLᵀ when indefinite ([Nonlinear static](nonlinear-static.md)) |

![The rook pivot search](../../media/manual/linear-solvers-3.svg)

### Memory limits

`solve.static`, `solve.nonlinear`, `solve.heat`, `solve.frequency` and `solve.buckle` accept
`max_factor_entries`: 400 000 000 by default, or an explicit whole number from 1 to
2 000 000 000. `{"step": "Load", "max_factor_entries": 500000000}` deliberately raises the cap
for that solve; it does not change the default of later solves. The other work and iteration caps
can only be lowered.

Each factor value is 8 bytes: the default allows 3.2 GB of values, and the ceiling 16 GB. These
are factor values alone, not total RAM: indices, assembled matrices, numeric workspace and other
solver storage need additional memory. Raising the cap neither reserves memory nor guarantees
that it is available. A factor over the cap is reported before its values are allocated, with the
required entry count and its byte estimate: "the factorization needs N entries (B bytes of factor
values alone; indices and workspace are additional); the limit is M entries (refine less or raise
max_factor_entries)". An allocation that fails is "not enough memory for the factorization (N
entries)". 32-bit builds remain limited by addressable allocation sizes. The ceiling is defined
once, as `FACTOR_ENTRY_CEILING` in `crates/engine/src/cmd/solve.rs`.

The same entry budget guards two more stores:

- **Constraint reduction** (`constraints::Budget`, 8 × `max_factor_entries` bytes over the matrices
  being reduced): a held distributing coupling makes every coupled node depend on every other, a
  dense (3N)² block, estimated from the rows and refused before anything is built ("distributing
  coupling `name` couples N nodes and its reference point is held …; split the surface into
  smaller couplings, hold the surface itself, or use a kinematic coupling"); every other reduction
  counts its entries and stops at the budget before allocating ("reduced matrices need more
  than …").
- **The Lanczos basis** (vectors and their mass products, 2 n × vectors): at most
  `max_factor_entries` / (2 n) vectors, within `max_iterations`. At 6 000 000 equations this is 33
  vectors with the default cap and 166 at the ceiling.

## Algorithm

1. **Pattern.** `NodeGraph::new` from the elements', trusses', shells', beams' and discrete
   elements' node lists.
2. **Order.** `fill_reducing_order` on the mesh nodes: coordinate and graph dissection, the one
   with the smaller symbolic fill kept.
3. **Reduce.** `constraints::Reduction` eliminates prescribed values and constraint rows in that
   order (free unknowns numbered by the node order), within its memory budget.
4. **Symbolic.** Elimination tree, postorder, column counts, supernodes; the entry count against
   `max_factor_entries` (refused before allocation).
5. **Allocate** the supernodes' values (`try_reserve_exact`; a failure is the out-of-memory
   error).
6. **Numeric.** Subtrees in parallel, supernodes left-looking, dense panels and tiles; the run is
   cancellable between supernodes. A non-positive pivot (Cholesky) or a supernode without an
   acceptable pivot (LDLᵀ) stops with the column, which the solver turns into the node and
   direction: "the stiffness is singular at node N in direction …" (Cholesky) or "the matrix is
   singular at …: no pivot of its supernode is more than round-off (a mechanism, or a shift or
   load at an eigenvalue)" (LDLᵀ).
7. **Solve.** Forward and back substitution over the supernodes, the permutation applied; LDLᵀ
   solves with D's 1 × 1 and 2 × 2 blocks, optionally refined.

There is no iteration in the factorization itself; the bounded loops are the rook search (64
steps), the refinement (a fixed count, stopping when the residual stops shrinking) and the
Lanczos basis (`max_iterations`).

## Output

The solve replies report `dofs` (all unknowns), `equations` (the free ones factored) and
`factor_entries` (entries of L); frequency and buckling replies add `lanczos_vectors`. Nonlinear
steps record each iteration's negative pivot count, and a shifted buckling step its
`below_shift` count, both from the LDLᵀ inertia.

## Worked example

The `bracket-hex` sample (5 778 hex20 elements, the model of the [Linear static](linear-static.md#worked-example)
example) solved on 1, 4 and 14 threads, its results written to CSV after each solve and compared
byte for byte, and solved once more with a factor cap below what its factor needs.

| | |
|---|---|
| **Model** | `file.sample` `bracket-hex`: filleted steel L-bracket, hex20 of 3 mm, 27 708 nodes; upright's end clamped, 1 kN on the foot's end; step `Load` (static) |
| **Solves** | `solve.static` with `threads` 1, 4 and 14; then with `max_factor_entries` 50 000 000 |
| **Compared** | `results.export.csv` after each solve: probe results of nodes 1 to 10 000 (coordinates, displacement, reaction, stress, principal and von Mises stress, 17 significant digits) and the reactions of set `top` |
| **Machine** | Apple M4 Max (14 cores), macOS 26, CPU |

```json
[
  {"command": "file.sample", "params": {"name": "bracket-hex"}},
  {"command": "solve.static", "params": {"step": "Load", "threads": 1}},
  {"command": "results.export.csv", "params": {"what": "probe", "step": "Load", "nodes": [1, 2, 3, "…", 10000], "path": "t1.csv"}},
  {"command": "results.export.csv", "params": {"what": "reactions", "step": "Load", "sets": ["top"], "path": "r1.csv"}},
  {"command": "solve.static", "params": {"step": "Load", "threads": 4}},
  {"command": "results.export.csv", "params": {"what": "probe", "step": "Load", "nodes": [1, 2, 3, "…", 10000], "path": "t4.csv"}},
  {"command": "results.export.csv", "params": {"what": "reactions", "step": "Load", "sets": ["top"], "path": "r4.csv"}},
  {"command": "solve.static", "params": {"step": "Load", "threads": 14}},
  {"command": "results.export.csv", "params": {"what": "probe", "step": "Load", "nodes": [1, 2, 3, "…", 10000], "path": "t14.csv"}},
  {"command": "results.export.csv", "params": {"what": "reactions", "step": "Load", "sets": ["top"], "path": "r14.csv"}},
  {"command": "solve.static", "params": {"step": "Load", "max_factor_entries": 50000000}}
]
```

(`"…"` stands for the node labels 4 to 9 999, written out in the run; a probe takes at most
10 000 nodes.) Results from the `solve.static` replies and the files. Wall times are of
`solve.static` (assembly, ordering, factorization, solve and stress recovery) in two runs, with
other jobs running on the machine (load average about 18).

| Quantity | 1 thread | 4 threads | 14 threads |
|---|---|---|---|
| `dofs` | 83 124 | 83 124 | 83 124 |
| `equations` | 82 065 | 82 065 | 82 065 |
| `factor_entries` | 91 826 205 | 91 826 205 | 91 826 205 |
| `max_displacement` | 0.5511123143618812 mm | the same | the same |
| `max_von_mises` | 295.0983647905056 MPa | the same | the same |
| `reaction_total` y | 999.999999981553 N | the same | the same |
| Probe CSV, 10 000 nodes (3 481 210 bytes) | — | identical to 1 thread | identical to 1 thread |
| Reaction CSV, set `top` | — | identical to 1 thread | identical to 1 thread |
| Wall time, two runs | 11.2 s, 11.8 s | 3.9 s, 3.6 s | 3.8 s, 2.5 s |

The solve under a cap of 50 000 000 entries is refused in 0.7 s, before the numeric
factorization, with:

```text
limit reached: the factorization needs 91826205 entries (734609640 bytes of factor values alone;
indices and workspace are additional); the limit is 50000000 entries (refine less or raise
max_factor_entries)
```

![Displacement magnitude on the deformed hex20 bracket, the model whose factor of 91.8 million entries the threads share](../../media/manual/examples/structural-linear-solvers.svg)

**What to look at.** The factor holds 91.8 million entries, 0.73 GB of values, for 82 065
equations, about 1 120 entries per equation. Every entry of L comes from one task with a fixed sequence of operations,
so the three solves write the same files to the last digit while the time falls about threefold
from one thread to four. The refused solve names the entries the factor needs and their bytes,
computed from the symbolic factorization alone; `max_factor_entries` of at least 91 826 205
lets it through.

## Verification

| Case | Check | Test |
|---|---|---|
| Sparse SPD factor and solve; cap refused before allocation; singular matrix located; cancellation; empty and 1 × 1 matrices | exact solution, located errors | `cholesky::tests::factor_and_solve_a_sparse_spd_matrix` |
| Supernodes wider than two panels, 1, 2 and 7 threads | residual < 10⁻¹², factor bit-identical on every thread count | `supernodes_panels_and_threads` |
| Indefinite shifted systems | inertia equals the eigenvalue count below σ; residual < 10⁻¹⁰; on an SPD matrix the same fill as Cholesky and no negative pivot | `ldlt::tests::indefinite_systems_solve_with_their_inertia` |
| Saddle points with zero diagonals | 2 × 2 pivots; inertia (K's size, the constraints') | `two_by_two_pivots_for_zero_diagonals` |
| Singular Laplacian, NaN, caps | located `Singular` errors; cap before the numeric factorization | `singular_and_hostile_matrices_are_located_errors` |
| LDLᵀ on 1 and several threads | values and pivot order bit-identical | `same_factor_on_any_number_of_threads` |
| `Auto` selection and refinement | Cholesky for SPD, LDLᵀ otherwise with the right inertia; refined solutions | `auto_selection_and_refinement` |
| LDLᵀ cost on a hex8 stiffness pattern of 140 000 equations (1.5 · 10⁸ entries) | 1.1 times the Cholesky factor; 1.2 times on a shift with 10 600 negative eigenvalues and 7 700 2 × 2 pivots | `speed_against_cholesky` (heavy) |
| Orderings | a permutation with separators last; graph dissection deterministic; less fill than planes on a thin curved wall; a block keeps the better ordering | `ordering::tests` |
| Constraint budget | a held distributing coupling of 12 321 nodes refused with its size before allocating; byte counts that overflow refused | `oversized_reductions_are_refused_before_allocating`, `reduced_matrix_budget_refuses_byte_counts_that_overflow` |
| Lanczos basis under a raised cap | 33 and 166 vectors at 6 000 000 equations | `raised_factor_budget_keeps_the_lanczos_vector_limit` |
| Cap error text | entry count and byte estimate before allocation | `factor_limit_estimates_bytes_without_allocating_the_factor`, `heat_factor_cap_reports_value_bytes_before_allocation` |

Every implicit solver test and oracle case runs on these factorizations; the oracle fixtures fail
CI on any change of a frozen result beyond its band (10⁻⁹ of the value).

## Limitations

- Direct solvers only: memory grows with the factor (O(n^{4/3}) entries for a 3D mesh under
  nested dissection); there is no iterative solver and no out-of-core factorization.
- LDLᵀ keeps a static structure: L's entries below a supernode are not bounded by the search
  (refinement covers accuracy where it matters), and a nonsingular saddle point whose zero
  diagonal couples only to rows of a later supernode would stop as singular. No FieldCraft system
  has one: hybrid pressures are condensed per element.
- A zero eigenvalue is never counted in the inertia: a zero pivot is a located error.
- Unsymmetric matrices are not factored: the skew part of a follower pressure's load stiffness and
  of friction is left out of the tangent ([Nonlinear static](nonlinear-static.md)).

## Implementation

| File | Contents |
|---|---|
| `crates/solver/src/sparse.rs` | `Csr`, `NodeGraph`, `add_block`, `add_scalar`, `add_diagonal` |
| `crates/solver/src/ordering.rs` | `nested_dissection`, `graph_dissection`, `fill_reducing_order`, `fill` |
| `crates/solver/src/cholesky.rs` | `UpperColumns`, `Symbolic`, `Factor`, `FactorError`, `PIVOT_TOLERANCE`, `PANEL`, `TILE` |
| `crates/solver/src/ldlt.rs` | `LdltFactor`, `Inertia`, `Kind`, `Symmetric` (`new`, `solve`, `solve_refined`, `inertia`) |
| `crates/solver/src/constraints.rs` | `Reduction`, `Budget` |
| `crates/solver/src/analysis.rs` | `pattern`, `factor_raw`, `factor_symmetric`, `factor_error`, `budget`, `lanczos_basis` |
| `crates/solver/src/lib.rs` | `Limits`, `Error::factor_limit` |
| `crates/engine/src/cmd/solve.rs` | `limits`, `FACTOR_ENTRY_CEILING`, `threads` |

## References

- A. George. *Nested dissection of a regular finite element mesh.* SIAM J. Numer. Anal. 10(2):345–363, 1973.
- A. George, J. W. H. Liu. *An automatic nested dissection algorithm for irregular finite element problems.* SIAM J. Numer. Anal. 15(5):1053–1069, 1978.
- N. E. Gibbs, W. G. Poole, P. K. Stockmeyer. *An algorithm for reducing the bandwidth and profile of a sparse matrix.* SIAM J. Numer. Anal. 13(2):236–250, 1976.
- H. D. Simon. *Partitioning of unstructured problems for parallel processing.* Computing Systems in Engineering 2:135–148, 1991.
- J. W. H. Liu. *The role of elimination trees in sparse factorization.* SIAM J. Matrix Anal. Appl. 11(1):134–172, 1990.
- J. W. H. Liu, E. G. Ng, B. W. Peyton. *On finding supernodes for sparse matrix computations.* SIAM J. Matrix Anal. Appl. 14(1):242–252, 1993.
- T. A. Davis. *Direct Methods for Sparse Linear Systems.* SIAM, 2006.
- E. Ng, B. W. Peyton. *Block sparse Cholesky algorithms on advanced uniprocessor computers.* SIAM J. Sci. Comput. 14(5):1034–1056, 1993.
- E. Rothberg, A. Gupta. *An efficient block-oriented approach to parallel sparse Cholesky factorization.* SIAM J. Sci. Comput. 15(6):1413–1439, 1994.
- J. R. Bunch, L. Kaufman. *Some stable methods for calculating inertia and solving symmetric linear systems.* Math. Comp. 31(137):163–179, 1977.
- C. Ashcraft, R. G. Grimes, J. G. Lewis. *Accurate symmetric indefinite linear equation solvers.* SIAM J. Matrix Anal. Appl. 20(2):513–561, 1998.
- O. Schenk, K. Gärtner. *On fast factorization pivoting methods for sparse symmetric indefinite systems.* Electron. Trans. Numer. Anal. 23:158–179, 2006.
- I. S. Duff, J. K. Reid. *The multifrontal solution of indefinite sparse symmetric linear equations.* ACM TOMS 9(3):302–325, 1983.
- J. W. H. Liu. *The multifrontal method for sparse matrix solution: theory and practice.* SIAM Review 34(1):82–109, 1992.
- E. Anderson et al. *LAPACK Users' Guide.* SIAM, 1999.
- N. J. Higham. *Accuracy and Stability of Numerical Algorithms*, 2nd ed. SIAM, 2002.

### Factor-only byte budgets (library API)

`cholesky::Factor::estimate_memory` analyzes the actual symbolic pattern and reports factor
values, index storage, upper/reordered matrix storage and workspace separately.
`Factor::with_memory_limit` and `ldlt::Symmetric::with_memory_limit` refuse an insufficient
factor budget before allocating numeric values, while retaining the independent entry cap.
The real LDLᵀ pivot buffers are included; this API does not estimate native complex factors.

Workspace estimates count the largest products/solve buffers per scheduled subtree, sum
subtrees that may run concurrently, and compare with the sequential top and symbolic phases.
Indices use actual vector capacities, including the structure copied into the returned factor.
Saturating byte arithmetic is rejected rather than accepted as a small wrapped estimate.
The existing constructors and their floating-point operation order are unchanged.

This budget covers the factorization only, **not the process**. The caller must separately
reserve its full and reduced CSR matrices, model, retained results and Lanczos basis before
assigning the remaining budget. Symbolic storage is still allocated for the estimate; the
numeric budget check is not a pre-assembly memory guard. Allocator overhead and RSS are not
measured by this API. Whole-solve `max_memory_gib` command admission remains pending; the
existing command defaults and ceilings still apply.

## Diagnostic timing

Set `FIELDCRAFT_SOLVER_TIMING=1` before starting a native process to emit tab-separated
`FIELDCRAFT_TIMING`, stage name and wall seconds on stderr. It is disabled by default and
sampled once per process. The diagnostic covers ordering, symbolic analysis, factor allocation,
Cholesky numeric factorization and back substitution, plus static/frequency calls and the
Lanczos start, block operator, projection, reorthogonalization, dense eigenproblem, convergence
check and reconstruction. It does not alter arithmetic, tolerances or partition sizes.

Spans are **inclusive**: Lanczos contains its sub-stages and concurrent `backsolve_inclusive`
records overlap inside `lanczos_block_operator`. Do not sum those records as wall time. A span
also appears when its scope returns an error, so check the solve result. Logging overhead and
other machine workloads affect diagnostic times; use separate controlled runs for performance
acceptance. The timer does not measure peak memory or loading/saving of model files.

For `solve.static`, `solve.frequency` and `solve.harmonic`, the engine also emits
`FIELDCRAFT_COMMAND_TIMING`, command id and wall seconds for the complete `Session::execute`
call, including assembly, factorization, solution, summary construction and in-session result
storage. This is the primary diagnostic time. Separate `file.open` and `file.save` commands and
external RPC transport are excluded. Solver spans explain this total; they do not replace it.
An error still emits a partial elapsed span; only successful commands are comparable.

For supernode diagnostics, also set `FIELDCRAFT_SUPERNODE_TIMING=1`. Cholesky emits each
supernode's index, row count, width, inclusive wall time, assembly time, dense-product time,
padded product FLOP count, below-panel solve time, and the product time within assembly as
`FIELDCRAFT_SUPERNODE` records.
`FIELDCRAFT_SUPERNODE_TOP` lists the nodes processed serially above the subtree tasks. Assembly
includes its dense products; concurrent node times must not be summed as factorization wall
time. Product throughput includes packing, allocation and scheduling, and counts each multiply
and add separately, including padded micro-tiles. These diagnostics do not alter arithmetic or
work partitions; logging itself can affect measured time.

Triangular solves also emit `FIELDCRAFT_TRIANGULAR_WORK` when solver timing is enabled:
right-hand-side columns, executed lanes (including any padding), equations, logical factor
coefficient visits, scalar multiply/subtract updates, and permuted scalar values. Each record
is emitted after a successful traversal. Counts use that traversal's actual dimensions; they
are not hardware memory-traffic counters. Multiply coefficient visits by eight for the logical
coefficient byte volume, not measured DRAM traffic.

The `backsolve_lane_group_inclusive` and `backsolve_inclusive` spans may overlap across
workers. Their sum is accumulated elapsed time, not CPU time or command wall time.

Set `FIELDCRAFT_STATIC_RESIDUAL=1` to report `FIELDCRAFT_STATIC_RESIDUAL` on stderr:
the reduced-system residual norm, right-hand-side norm, and their ratio, computed
from the actual constrained stiffness and solved displacement. A zero right-hand-side
norm makes the ratio undefined (`NaN`). This opt-in diagnostic retains the reduced
matrix and right-hand side; its time and memory belong to the complete static command.
It does not add fields to the result archive.

With solver timing enabled, eligible wide Cholesky updates emit
`FIELDCRAFT_DENSE_WORK`: cache-panel variant flag, output rows, output columns,
reduction width, padded multiply/add FLOPs, and an upper bound on extra panel-merge
additions including padded columns. Counts describe logical work, not hardware
traffic. `dense_update_reference` spans can overlap and their sum is accumulated
elapsed time, not factorization wall time or CPU time. The cache-panel experiment
is retained only in solver tests; production uses the original full-width update.

For a controlled dense-dispatch diagnostic, set `FIELDCRAFT_DENSE_PROFILE=1`.
`FIELDCRAFT_DENSE_PROFILE` records output rows/columns, reduction width, tile count,
Y packing wall time, dispatch-to-join wall time, and accumulated X packing, kernel,
and whole-task elapsed times, followed by the longest task elapsed time. The kernel
interval includes output allocation and storing its dense tile. Task intervals can
include preemption and overlap each other: they are neither CPU time nor pure waiting.
Dispatch-to-join includes actual work, scheduling and joining; subtracting accumulated
task times from it does not measure scheduler overhead.

With that profiling enabled, `FIELDCRAFT_DENSE_SERIAL=1` executes the identical dense
output tiles sequentially. It changes only their dispatch, retaining tile boundaries,
per-element arithmetic and the supernode dependency order. This is a diagnostic
control, not a production optimization. Subtree parallelism remains active, so the
comparison measures the net effect of dense-task concurrency in that context.

Supernode timing also emits `FIELDCRAFT_SUPERNODE_DETAIL`: node index, accumulated
product-to-factor writeback wall time, written scalar entries, and diagonal panel
factorization time. The writeback count excludes the separate contiguous below-panel
copy. Existing assembly and product intervals include their nested operations; do
not add them together. `numeric_subtrees` and `numeric_top` partition the factor's
subtree-dispatch and serial top-node phases, apart from setup and error handling.

## Experimental multilevel ordering

There is an opt-in experimental third candidate for `fill_reducing_order`, in
`ordering_multilevel.rs`. Default settings retain the coordinate and graph-level
candidates. The experiment contracts heavy edges, grows a balanced coarse partition,
then refines it during uncoarsening with locked single-vertex moves and best-prefix
rollback. Crossing edges yield a vertex separator through bipartite matching and
minimum vertex cover; child parts precede their separator. Leaves of at most 64
vertices use quotient-graph approximate-degree bounds. This is a bounded variant,
without supervariable amalgamation or a claim to reproduce another AMD package.

Visits, adjacency, coarse ids, eight starting vertices and all ties are deterministic.
There are at most 32 coarsening levels, four refinement passes per level and 64
separator levels. Disconnected components are handled separately; dense or
unsplittable parts fall back to their stable input order. The third candidate is
enabled for standalone linear static and unprestressed frequency solves, with numerical, convergence and thread-determinism checks.
Complete solve-command measurements determine the benefit. Include ordering and diagnostics
in that total; a local ordering budget alone cannot determine the net benefit.
Small separators alone are not enough.

The deterministic least-fill choice compares all three candidates and keeps the
existing coordinate/graph candidate on a tie. `FIELDCRAFT_ORDERING_EXPERIMENT=legacy`
is a diagnostic override for comparisons with the original two-candidate selection;
no environment setting can enable the third candidate in an excluded analysis.
`FIELDCRAFT_ORDERING_TRACE=1` reports separator depth, part/child/separator sizes,
coarsening levels and cut weight. The ignored `measure_multilevel_graph` harness
reports exact node-graph column counts, three-DOF expansion and actual supernode
row fronts separately. None of these diagnostics changes the result-file format.

Algorithm sources: [Karypis and Kumar (1998)](https://people.eecs.berkeley.edu/~jordan/sail/readings/karypis-kumar.pdf),
[Fiduccia and Mattheyses (1982)](https://ics-websites.science.uu.nl/docs/vakken/ea/pdf/FidMath_GBP.pdf),
and [Amestoy, Davis and Duff (1996)](https://people.engr.tamu.edu/davis/publications_files/An_Approximate_Minimum_Degree_Ordering_Algorithm.pdf).
This is an independent implementation; no external partitioner is linked or vendored.

Frequency solves automatically assign the existing fixed-seed
random stream to the same physical free degrees of freedom as the original
coordinate/graph selection, including random columns replacing lost rank. The
factor may use a different equation order, but the physical starting block remains
fixed. It does not relabel returned modes, aggregate mode clusters or loosen any
convergence check. Building this coordinate map costs time and must remain inside
the complete command measurement. No opt-in flag is required for the supported
production paths.

Diagnostic controls for isolating finite-iteration effects are
`FIELDCRAFT_LANCZOS_TOLERANCE` (frequency only, bounded to 1e-13 through the original
1e-9, default unchanged) and `FIELDCRAFT_RITZ_TRACE=1` (basis size, Ritz index,
transformed relative residual bound and tolerance). A transformed residual bound
is not the physical residual `||K phi - lambda M phi|| / ||K phi||`. Tighter
convergence is a diagnostic reference, not permission to replace a frozen baseline;
record its extra iterations and full command time. These controls do not change
result identifiers or payloads.

The production choice is scoped to standalone linear static and unprestressed
frequency assembly. Buckling and its preload, nonlinear, heat, coupled thermal,
prestressed analyses and their preloads, harmonic/random/transient operators,
response spectra and sensitivities retain legacy ordering, including when the
experimental environment variables are set. The first unrestricted diagnostic moved two unrelated oracle fixtures; those failures
are not waived, and the third ordering is not validated for those analysis paths.

If a trial has adjacent squared frequencies separated by at most `sqrt(epsilon)`
relative to their scale, or a mode's summed translational effective mass is at
most `128 * epsilon` of total mass, the entire frequency solve is repeated with
legacy ordering. These are conservative fallback triggers, not relaxed result
acceptance bands. No modes are relabeled, rotated or combined. The first solve's
storage is released before the retry. Warnings report its actual iterations and
basis size; convergence fields describe the returned legacy solve. Full command
time includes both attempts. A comparison override requesting legacy ordering
never retries. This preserves individual modes in numerically ambiguous cases
at the cost of extra work for those cases.

A model containing any pretension section retains legacy ordering for every step,
including a standalone static or frequency step without a named base state.

The FM refinement checks the remaining partition capacity against the smallest
vertex weight before searching a gain heap. If no vertex can fit, it returns no
move without draining and restoring the heap. This preserves the weighted gain
ordering, stale-entry rules and fixed tie breaks. The diagnostic override
`FIELDCRAFT_ORDERING_REFINEMENT=baseline` retains the original heap scan for
same-program comparisons; it does not select a different permutation.

`FIELDCRAFT_ORDERING_PROFILE=1` reports exclusive phase wall-time totals for
canonical/induced graph construction, coarsening, initial growth, refinement,
separator extraction, leaf ordering and uncoarsening. `recursion_other` is the
remaining control/allocation/component work; `total` covers the complete
multilevel ordering. Input file decoding and later symbolic factor counts are
outside that total. FM query counts and impossible-query time are diagnostic
subsets of refinement, not additional wall time to add to the phase sum.

Multilevel graph storage uses CSR rows. Heavy-edge contraction accumulates each
coarse row with vertex markers before sorting its unique neighbors; matching
visits, aggregate lists, markers and accumulation buffers are reused between
levels and recursive subgraphs. Fine-to-coarse maps and layer graphs remain
owned until uncoarsening completes. Induced subgraphs share a resettable global
to local index map. The canonical graph uses two-pass symmetric counting and
filling, independent row sorting/deduplication and ordered compaction. Parallel
row processing has one writer per row; recursive subgraphs are still serial.

The contraction ratio is the number of coarse vertices divided by the number
of fine vertices. A ratio of 0.9 or higher stops coarsening; graphs of at most
128 vertices also stop, with at most 32 coarsening levels. Recursive leaves of
at most 64 vertices use the quotient-graph AMD approximation. The diagnostic
`FIELDCRAFT_ORDERING_STORAGE=baseline` retains the row-vector implementation.

The CSR heavy-edge matcher visits a fixed integer Fisher–Yates permutation
(SplitMix64, seed `0x243f6a8885a308d3`), with the minimum original vertex id
breaking equal edge-weight ties. This is independent of thread count and input
adjacency order. `FIELDCRAFT_ORDERING_VISITS=degree` retains degree/id visits
with CSR storage; the storage baseline always retains the original degree/id
visits. Neither diagnostic broadens the analysis types eligible for multilevel
ordering; models with pretension sections retain legacy ordering.
