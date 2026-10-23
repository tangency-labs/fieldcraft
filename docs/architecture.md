# Architecture

FieldCraft turns geometry into finite element models, writes them for the solvers engineers
already use, and solves linear problems itself.

```
SolidCraft B-rep ──tessellate──▶ boundary surface ──▶ volume mesh ──▶ FE model ──▶ Abaqus .inp
 (solidcraft-kernel)           (welded, watertight,   tet (Delaunay)   materials,     Nastran .bdf
                                face-tagged)          hex (core goal)  sections,      LS-DYNA .k
                                                                       BCs, loads,  ──▶ own linear
                                                                       steps            solver
```

## Crates and layers

`cargo xtask layers` enforces the table in `xtask/src/layers.rs`: a crate may depend only on
crates in strictly lower layers.

| Layer | Crate | Role | Status |
|---|---|---|---|
| L0 | `fieldcraft-model` | Solver-neutral FE model: nodes and elements by label, their geometry tags, sets, materials, sections, couplings, BCs, loads and steps, the model check and the resolved view the writers read. Knows no file format and no geometry type. | started |
| L0 | `fieldcraft-parallel` | Deterministic data parallelism for the solvers: fixed work partitions on a rayon pool, results in block order, the thread-count setting (`FIELDCRAFT_THREADS`, `threads`); sequential in the browser. No workspace dependencies. | started |
| L1 | `fieldcraft-geometry` | The only door to SolidCraft: opens `.solidcraft` designs (rebuilt from their timeline, bodies with persistent face names), B-rep bodies, measures, tessellation into a face-tagged `SurfaceMesh`; face-tagged 3MF files read into the same type without the kernel (docs/formats/tagged-mesh.md). No SolidCraft type in its API. | started |
| L1 | `fieldcraft-mesh` | Surface preparation (weld, watertightness, feature edges), robust predicates, Delaunay tetrahedralization and refinement, Voronoi-based optimization, hex meshing, quality metrics. | started |
| L1 | `fieldcraft-fem` | Element kernels shared by every structural solver: the element table (tet4, tet10, wedge6, hex8 full and incompatible modes, hex20 full and reduced), shape functions, quadrature, element stiffness, consistent and lumped mass, body and pressure loads, stresses and their extrapolation, isotropic and orthotropic elasticity, the uniform-strain (one-point) hex8 and tet4 kinematics, hourglass shape vectors and element frequency bounds of the explicit solver, objective stress rates and J2/Johnson–Cook return mapping. Pure functions, no solver loop. | started |
| L1 | `fieldcraft-decks` | Abaqus input, Nastran bulk data and LS-DYNA keyword readers and writers, to and from `model`. Text out, no file system; `deck.write` saves it. | started |
| L1 | `fieldcraft-render` | Display data without a GPU dependency: the visible faces, element edges and nodes of a model, the orbit camera, a CPU rasterizer for headless snapshots, geometry surfaces with their face boundaries, deformed shapes and banded contour plots of nodal results, slices, surface values and streamlines of flow fields on grids of cells, picking. | started |
| L2 | `fieldcraft-solver` | Implicit solver: assembly from the resolved model, boundary conditions, sparse direct solution, linear static and modal analysis, stresses and reactions. | started |
| L2 | `fieldcraft-explicit` | Explicit dynamics: central-difference time integration with lumped mass, one-point hex8 with hourglass control, nodal-pressure tet4, bulk viscosity, rigid walls, large deformation with objective rates, plasticity, damage and erosion, energy account, histories and field frames. | started |
| L2 | `fieldcraft-lbm` | Lattice Boltzmann flow solver on the CPU: D3Q19 two-relaxation-time or recursive regularized scheme, bounce-back and slip walls, Zou–He inlets and outlets, convective outlets, obstacles voxelised from the face-tagged surface with exact predicates, momentum-exchange forces per face; a D3Q7 temperature lattice carried by the flow (held temperatures, heat fluxes, adiabatic faces, Boussinesq buoyancy, the heat into every face) and conjugate heat transfer: conduction through the obstacles that conduct, coupled to it at their faces. Uses `mesh` for the predicates and the boundary surface. | started |
| L3 | `fieldcraft-fsi` | Fluid–structure coupling between the L2 solvers (later). | planned |
| L4 | `fieldcraft-doc` | The document: geometry, mesh, model, results, history. | planned |
| L5 | `fieldcraft-io` | `.fieldcraft` files (docs/formats/fieldcraft.md): text model definition, binary mesh, geometry references with content hashes, the journal; atomic saves; autosave storage for crash recovery; read and write roots that confine an automation client's file access (`roots`). Bytes in, bytes out; the file system only in `fs`, `recovery` and `roots` (in the browser, files kept in memory: `memory`). Result files later. | started |
| L6 | `fieldcraft-engine` | Commands (`CommandSpec`), `Session::execute` with rollback and the panic guard, scripts; every file a command names checked against the session's roots (`Session::file_access`); headless pictures of the session (`view`: model, geometry, deformed shape and contours on the CPU rasterizer). | started |
| L7 | `fieldcraft-gpu` | The viewport's GPU pipelines (lit faces, wide lines, MSAA, depth) without a window toolkit, the adapter policy (`FIELDCRAFT_GPU`: hardware GPU, software GPU, CPU rasterizer) and offscreen pictures (`fieldcraft-cli snapshot`). The window draws through the same pipelines. | started |
| L7 | `fieldcraft-lbm-gpu` | The lattice Boltzmann solver's GPU backend: the CPU scheme's step (TRT or regularized collision with Guo forcing, the per-link wall rules, Zou–He boundaries and convective outlets, momentum exchange, and the temperature lattice with its wall conditions, buoyancy and heat exchange, and the conducting solids with their implicit interface links) as wgpu compute shaders, driven through `lbm`'s `Backend` trait from the flow the CPU prepared; single precision, one population buffer streamed in place (the AA pattern, `docs/lbm-gpu-performance.md`), deterministic per device, checked against the CPU within its tests' tolerances; opt-in per solve (`solve.fluid` with `backend: "gpu"`). Opens its device with `gpu`'s adapter policy. | started |
| L7 | `fieldcraft-mesh-gpu` | Projection smoothing's GPU device: the coloured Winslow sweeps' interior nodes (Winslow's update or the neighbour mean, the guarded step on the scaled Jacobian and edge ratio around the node) as a wgpu compute shader, driven colour by colour through `mesh`'s sweep device while the host moves the boundary nodes on the geometry; single precision about the mesh's centre, deterministic per device, the same verdicts as the CPU's coloured sweeps with a smallest scaled Jacobian within 0.02 of theirs; opt-in per projection (`mesh.hex.project` with `backend: "gpu"`). Opens its device with `gpu`'s adapter policy. | started |
| L7 | `fieldcraft-ui-egui`, `fieldcraft-mcp` | Desktop UI in the Craft family's layout (title bar with quick actions, workspace ribbon from the catalogue, browser, viewport with navigation bar and orientation gizmo, context panel, command search and command line, journal timeline; native file pickers as host services, background opening and solving; results display with legend and probe; offscreen capture of the whole window) and control channel requests (docs/control-protocol.md); MCP server over stdio (docs/mcp.md): tools over engine commands, a headless session or the running app through its control channel, read and write roots. | ui, mcp started |
| testkit | `fieldcraft-testkit` | Test tooling, a dev-dependency only: the meshing test corpus (`tests/corpus/*/recipe.json`, analytic shapes) and its runner (`cargo xtask corpus`). | started |
| apps | `fieldcraft`, `fieldcraft-cli`, `fieldcraft-web` | Desktop, command line, browser (docs/web.md). | started |

Rules beyond the layer order:
- Nothing below L7 uses egui, eframe, winit, wgpu or rfd.
- Only `fieldcraft-geometry` uses `solidcraft-*`, pinned by rev in the workspace `Cargo.toml`
  ("SolidCraft pin" below).
- Everything below the apps builds for `wasm32-unknown-unknown`.
- `fieldcraft-geometry`'s non-default `scripts` feature builds designs from SolidCraft command
  scripts (`script_design`, through `solidcraft-engine`) for the test corpus; the library opens
  designs without it.

### SolidCraft pin

`solidcraft-*` come from https://github.com/tangency-labs/solidcraft (formerly
bherbruck/solidcraft, which redirects there), pinned in the workspace `Cargo.toml` to the commit
of a SolidCraft release tag (now v0.18.2). SolidCraft's `[patch.crates-io]` (its fixed truck
crates) does not reach dependents, so the workspace repeats it from the same rev; bump both
together, and keep the list in step with SolidCraft's (v0.16 added `truck-geometry`). The same
goes for SolidCraft's patched `egui-wgpu` (MSAA clamped to what the device supports, 1x on
OpenGL adapters, so the desktop app's 4x never fails pipeline validation at start).

The pin decides which designs FieldCraft opens: SolidCraft's reader migrates every older design
format and refuses a newer one, and SolidCraft releases often raise the format. A design newer
than the pin reads (`fieldcraft_geometry::DESIGN_FORMAT`) fails with "this design was saved by a
newer SolidCraft (format N); this FieldCraft reads up to M. Update FieldCraft to open it" (#605:
FieldCraft 0.7.0 read up to format 21, SolidCraft 0.16.2 saves 27; #655: a FieldCraft on 0.16.2
reads up to 27, SolidCraft 0.17 and 0.18 save 28).

A newer pin can also name the same geometry otherwise: SolidCraft 0.16 names the pieces of a
split face from their history (`F9:face:4` became `F9:section`), so a model saved before names
faces, edges and vertices the design no longer has (0.17, 0.18.1 and 0.18.2 renamed none of the
corpus's and samples' faces, edges or vertices: #655). Opening a model or a design never lets that
pass in silence (`crates/engine/src/cmd/names.rs`): a stale face the model's mesh lies on is
found again where the mesh's facets tagged with it land on the new tessellation (all on one face
whose area agrees, no other tag's mesh there), its edges and vertices through the faces they
separate, and every definition and mesh tag follows, after the renamed definitions are checked
to read back from a file; anything else is refused. A name that still resolves but whose mesh
lies on another face is reported (`mismatched`). The open lists
each rename (`renamed`, and the warnings the status line shows) and marks the model changed, so
it is saved deliberately; each definition still naming what the design lacks is listed
(`stale.references`) and asked to be picked again, and a remesh says which seeds and size
controls it could not apply.

Three checks keep the pin current:
- `crates/geometry/tests/design_formats.rs` opens every design in
  `crates/geometry/tests/designs` (the older fixtures and `solidcraft-<version>/`, each saved by
  that release) and checks that a later release's copy of a recipe gives the same bodies,
  volumes and face names as the older one.
- `crates/geometry/tests/design_names.rs` freezes every face, edge and vertex name of the split
  corpus bracket (`tests/names/`): a pin bump that renames anything shows the renames there, to
  bless once opening older models finds them again.
- The scheduled `SolidCraft compatibility` workflow (`.github/workflows/solidcraft-compat.yml`,
  daily and on demand) downloads SolidCraft's latest release, saves the same recipes with it
  into `solidcraft-latest/` and runs those tests. Red means a SolidCraft release saves designs
  this FieldCraft can't open (or opens differently): move the pin to that release, port the API
  changes in `fieldcraft-geometry`, add its designs as `solidcraft-<version>/`, and check the
  meshing corpus (`cargo xtask corpus`), the heavy tests and the oracle stay the same.

## Meshing pipeline

1. **Boundary surface.** SolidCraft tessellates each body; every triangle carries its B-rep face,
   so feature edges and face-based loads and BCs survive meshing. Faces of an opened design also
   carry SolidCraft's persistent names (`F2:side:l4`, from the feature and sketch curve that made
   them), so a load on a named face follows it when the design is edited and opened again.
   Bodies of components are placed in the world through their occurrences (one body per
   placement), as SolidCraft displays and exports them.
2. **Surface preparation** (`fieldcraft-mesh::surface`, `mesh.surface.check`). Weld the tessellation's per-face
   vertex copies (at twice SolidCraft's linear tolerance), check the surface is closed, manifold, consistently and
   outwardly oriented (voids inwards) and free of self-intersections (exact triangle–triangle tests), and find
   feature edges (face boundaries, sharp dihedral angles). Defects come back located by triangle, vertex and
   position.
3. **Tetrahedral meshing.** Bowyer–Watson Delaunay on Shewchuk's adaptive predicates (symbolic
   perturbation for degeneracies), boundary recovery by the robust constrained Delaunay method of
   Diazzi et al. 2023 (implicit Steiner points on indirect predicates, gift-wrapping facet
   recovery; `fieldcraft-mesh::cdt`, command `mesh.tet.generate`), so every boundary triangle and
   its B-rep face tag survives; rounding is repaired by flips (§5.1). Then the quality stage
   (`fieldcraft-mesh::tet`): a curvature-graded size field, edge-collapse coarsening, constrained
   Delaunay refinement with boundary points projected onto the B-rep faces and curves (exact
   surfaces from SolidCraft where analytic, the tessellation otherwise), ODT smoothing
   interleaved as in Tournois et al. 2009, the Klingner–Shewchuk 2007 improvement schedule
   against slivers, and tet10 with mid-edge nodes on the geometry. Every step is a checked local
   replacement, so the mesh stays valid throughout.
4. **Hexahedral meshing (the core).** Grid-based all-hex first, the one family that meshes any
   input validly. `fieldcraft-mesh::octree` (command `mesh.hex.octree`): a 27-tree over the part,
   refined by curvature, thickness and narrow faces, balanced and cut by the 3-refinement
   node-marking templates (Schneiders 1996, Ito et al. 2009) into a conforming grid; the core of
   grid elements inside the part (exact inside test) is cleaned to a manifold boundary; the part's
   patches, sharp curves and corners are mapped onto it (Gao et al. 2019); one padding layer and a
   pillow under each patch carry the boundary and capture every sharp curve; the foldover-free
   energy (Garanzha et al. 2021) with a growing fitting term moves the boundary onto the geometry
   while every element stays valid; the worst elements are improved, and every element is proven
   valid by the certified Bernstein Jacobian check (Johnen et al. 2013; Marschner et al. 2020)
   before the boundary quads are tagged with their B-rep faces. Where a stage fails the grid is
   refined there and the mesh rebuilt. The dual schemes of Livesu et al. 2021 (fewer elements,
   concave transitions) are the next upgrade of the transitions. Then hex-dominant by tet
   recombination (Pellerin et al. 2018), then polycube (DualCube labelling) for blocky parts;
   frame-field/integer-grid-map methods stay research.
5. **Structured hex and face meshing.** Faces get all-quad meshes on shared, prescribed edge
   nodes (`fieldcraft-mesh::quad`): transfinite interpolation for four-sided faces (Gordon and
   Hall 1973), and for free faces a half-resolution constrained Delaunay triangulation,
   recombination (Lee and Lo 1994) and one uniform split that lands on the boundary nodes, then
   smoothing and untangling; every quad is checked with exact predicates. One `SizeField` trait
   (`fieldcraft-mesh::size`) sets interior sizes and edge divisions, with even counts per loop.
   Mapped and swept hex build on these face meshes (`fieldcraft-mesh::hex`, `mesh.hex.generate`):
   the B-rep topology is recovered from the face-tagged surface; a body is extruded (straight
   multi-level sweep: the levels' outlines imprinted into one cross-section, cut automatically
   into blocks (holes apart, O-grids, cuts from concave corners and fillets chosen by the quality
   of the blocks they leave), regions mapped, triangle primitives, O-gridded or free-meshed, then
   stacked), revolved (the same in a rotating
   frame), swept one-to-one through rings of linking faces (least-squares affine layers, Roca and
   Sarrate 2010), or mapped as a logical cube (Gordon and Thiel 1982); interval assignment
   couples the counts of shared curves. Nodes on analytic faces (`Shape::face_surfaces`) are
   exact, hex20 mid-edge nodes go onto the faces, boundary quads carry B-rep face tags, and the
   result passes the certified validity gate plus tiling and fidelity checks. A body no method
   fits gets each method's reason.
   Without geometry, `mesh.block`
   (`fieldcraft-mesh::block`) meshes one six-sided block from its eight corners and optional
   circular or elliptical edges by the same edge-based transfinite interpolation, as hex8, hex20,
   tet4/tet10 (Kuhn's six tetrahedra per cell) or wedge6, with grading; blocks join where their
   nodes coincide, face to face (mismatched divisions or kinds are refused), and their named
   faces, edges and vertices become geometry tags that sets select. It meshes the analytic
   shapes of the verification cases.
6. **Quality.** Scaled Jacobian, aspect ratio, dihedral angles; every generated mesh is checked
   for valid Jacobians and conserved volume before it reaches the model. Validity is one
   criterion shared with both solvers (`fieldcraft_model::admissible_jacobian`: normalized
   Jacobian above 1e-12): the gate proves it on the whole element, the solvers check it at their
   integration points, so a mesh the gate accepts is never refused by a solver.

## Model definition

`fieldcraft-model` holds the analysis model in no solver's vocabulary (the writers map it to
element names, cards and keywords):

- **Geometry tags.** Geometry is referred to by name only: a body and, inside it, the persistent
  name of a face, edge or vertex (`{"body": "Body1", "face": "F2:top"}`). A mesher tags what it
  makes: each element with the body (solids), face (shells) or edge (beams) it fills, each
  boundary element face with its geometry face, element edges with geometry edges, nodes with
  vertices (`Model::geometry_tag`, `tag_element`, `tag_element_face`, …; `mesh.tag` by hand).
  Replacing or removing an element drops its tags; a remesh clears them and writes new ones.
- **Sets** are defined by geometry names, explicit labels or both, and resolve against the
  current mesh into sorted nodes, elements and element faces. A face name also selects the pieces
  a kernel splits it into (`name#0`, `name#1`). Because the definition names geometry, a remesh,
  or an edit of the design followed by a remesh, keeps every set, section and load in place.
- **Element faces** are numbered per kind (`ElementKind::faces`): solids' faces are oriented with
  the right-hand normal pointing out; shells have a positive and a negative side.
- **Definitions** (materials, sections with their element formulation, reference points,
  couplings, boundary conditions, loads, steps) refer to each other by name. A reference must
  exist when it is made, and a definition still referred to (a load listed by a step included)
  cannot be deleted; what depends on
  the mesh is checked when the model is resolved.
- **Steps** list the boundary conditions and loads active in them. Nothing propagates from one
  step to the next.
- **Units** are a label (`mm-t-s-N-MPa`, …); numbers are never converted, and the check warns
  when a modulus, density or gravity is implausible for the chosen system.
- **Check and resolve.** `Model::check` lists errors and warnings, each naming its item (a set
  member not on the mesh, an element with no section or in two, a moment on solid nodes, a
  missing density for a frequency step, …). `Model::resolve` returns what writers and the solver
  read besides the definitions: every set resolved (each element then lies in exactly one
  section's set), the node label of every reference point, and every force, moment and traction
  as nodal values. It fails when the check finds an error.
- **Loads that keep their meaning on any mesh.** A force or moment on a set is a total by
  default: split by the consistent nodal weights `∫ N_i dA` of the set's faces (Gauss quadrature
  on each face's own linear or quadratic geometry), or equally over a set of nodes only. Applying
  it at every node is an explicit, mesh-dependent opt-in (`per_node`). A traction is force per
  area on faces. The expansion lives in `resolve`, so the deck writers and the solver share it.
- **Reference point labels** are stored with the point: chosen from the top of the label range
  when it is created, so remeshing or adding nodes and points never renumbers them.
- **Files.** A model deserializes by replaying its contents through the same checked methods the
  commands use, so a file cannot hold what commands would refuse.

## Implicit solver

`fieldcraft-solver` solves a model's linear static and frequency steps (`solve.static`,
`solve.frequency`; `solve.probe` reads the stored results). Methods and papers are cited in the
module docs; ADR 0001 records why the linear algebra is our own.

1. **System.** `Model::resolve()`; each solid section's elements take their material's
   elasticity matrix and the section's formulation (`fieldcraft_fem::element::ELEMENT_TYPES`, one
   row per element kind: tet4, tet10, wedge6, hex8, hex20); 2D sections (plane stress, plane
   strain, axisymmetric) take the 2D table (`fieldcraft_fem::plane::PLANE_TYPES`: tri3, tri6,
   quad4, quad8), chosen by the section kind, and the model then has two degrees of freedom per
   node: the 2D elements declare `ElementDofs::new(DofSet::PLANE)` (Ux, Uy) to the `DofLayout`. A 2D model
   lies in the XY plane (axisymmetric: x the radius, Y the axis; nodes on the axis are held
   radially), is loaded on element edges, and cannot mix with solids or couplings yet. Rod
   sections' rod2 and rod3 elements are trusses (`fieldcraft_fem::truss`). Point masses, rotary
   inertia and springs (`Resolved::masses`, `springs`, expanded per node by the model) assemble on
   the unknowns their nodes are given (`discrete::needs`, `discrete::dof_of`); dashpots are
   ignored with a warning (no velocities in static and frequency steps). Every reference point a
   coupling or a discrete element uses gets six degrees of freedom (translations, rotations) after
   the mesh's. Shell sections make quad4 (MITC4) and tri3 (MITC3) flat shells with a drilling
   rotation (`fieldcraft_fem::shell`; quad8 and tri6 shells are "not supported yet"), beam
   sections two-node Timoshenko beams with the profile's section properties and the section's
   orientation vector (`fieldcraft_fem::beam`); both declare `ElementDofs::new(DofSet::STRUCTURAL)`
   (`system::Structural`), so their nodes carry translations and rotations and the solids and
   trusses sharing those nodes add into the translations. Shells and beams need isotropic
   materials until material orientations exist, and are 3D only. The step's boundary conditions
   prescribe the kinds a node carries (rotations of solid-only nodes are ignored) and a coupled
   point's six directions; its loads become a load vector: forces, moments and tractions from
   `Resolved::nodal_loads`, forces and moments on coupled reference points, pressure as
   consistent face loads (shell faces too), gravity and body forces through the element shape
   functions (consistent forces and moments on shells and beams).
2. **Assembly** on the node graph (one dense block per node pair: 3 × 3 for solids, 2 × 2 in 2D
   models, up to 6 × 6 at shell and beam nodes; elements in label order, so sums are
   reproducible bit for bit).
3. **Ordering and elimination** (`constraints`). Nested dissection of the nodes by coordinate
   bisection. Couplings become linear constraint rows: kinematic = rigid link of each coupled
   translation to the point's translation and rotation; distributing = the point follows the
   least-squares rigid fit of the set's translations (equal weights, the deck writers' choice).
   Prescribed values and rows are eliminated master–slave in a fixed order, one dependent degree
   of freedom per row (the coupled node or the distributing point while it is free, else the
   largest coefficient), giving u = T_f q + T_p u_p; the reduced matrix T_fᵀ K T_f stays
   symmetric positive definite (no penalty, no multipliers). Free mesh equations follow the
   node order, reference point equations come last. Without couplings T is a selection and the
   numbers are those of plain elimination, bit for bit.
4. **Static solution.** Supernodal left-looking sparse Cholesky on the postordered elimination
   tree (the factor size is computed and capped before any value; subtrees and dense row tiles on
   threads with a fixed operation order, so the factor is the same on any thread count), reactions T_pᵀ (K u − f) (at nodes and at prescribed reference points, forces and
   moments), integration-point stresses extrapolated to the nodes and averaged, von Mises. Shells
   report section forces and bottom, middle and top stresses at their nodes (their top surface
   joins the nodal average), beams their end section forces; nodes report rotations and reaction
   moments. A zero pivot names the node or reference point and direction that can move freely.
5. **Frequency.** Shift-invert block Lanczos (block 6, full reorthogonalization) on
   K − σM, σ = 0 for a supported structure and slightly negative for a free one, so rigid-body
   modes come out at zero; mass-normalized shapes, participation factors and effective masses.
6. **Results** stay in the session per step, stamped with the session revision they were solved
   at; `Session::result` hands out only those of the current revision, so an edit makes them
   stale and undoing back to the solved state makes them current again.

Verification (`crates/solver/tests/verification.rs`, analytic only): displacement and traction
patch tests for every element and formulation (exact to round-off), cantilever tip deflection
against Timoshenko with convergence rates of the linear elements, Lamé's cylinder, Kirsch's
plate with a hole (Kt → 3), cantilever and free–free beam frequencies with six rigid-body modes.
`crates/solver/tests/plane.rs`: the 2D elements' displacement patch tests in plane stress, plane
strain and axisymmetry, uniform traction, a plane stress cantilever against Timoshenko with
convergence rates, Kirsch's plate, Lamé's cylinder in plane strain and axisymmetric, a rotating
disc (plane stress) and a rotating cylinder (axisymmetric) against Timoshenko & Goodier, masses
and rigid-body modes.
`crates/solver/tests/couplings.rs`: a block given a rigid motion through a coupled point
(kinematic and distributing), clamping through a coupled point equals clamping (static and
modes), tension through a distributing coupling against the equivalent pressure away from the
end, pure bending by a moment through a distributing coupling, torsion of a tube by a moment on
a kinematically coupled point (twist and shear stress). `crates/solver/tests/one_d.rs`: member
forces of determinate trusses (a plane Pratt truss, a space tripod) against the method of joints,
the three-bar truss, a rod3 constant-stress patch and a bar hanging under its weight, spring–mass,
skew springs, the two-mass chain, a bar with a tip mass and the torsional pendulum.
`crates/solver/tests/structural.rs`: membrane and bending patch tests of quad4 and tri3 shells
in a tilted plane, the simply supported plate under pressure against Navier's series (deflection
and first frequency), cantilever beams against Timoshenko (exact) and Euler–Bernoulli
frequencies, a portal frame against slope-deflection, solids continued by rods and skinned by
shells. `crates/solver/tests/benchmarks.rs` (published reference values, tagged with their
source, separate from the analytic gate): the Scordelis–Lo roof and the pinched cylinder.

## Explicit solver

`fieldcraft-explicit` integrates a model's `explicit_dynamic` steps in time (`solve.explicit`;
`solve.history` and `solve.probe` read the stored histories and frames). Methods and papers are
cited in the module docs and in docs/solvers/literature.md (Parts H and J). Both structural solvers
number their degrees of freedom with `fieldcraft_fem::DofLayout` (docs/fem.md).

1. **System.** `Model::resolve()`; hex8 sections with reduced integration become one-point
   uniform-strain hexahedra with Flanagan–Belytschko hourglass control, tet4 sections Bonet–Burton
   average nodal pressure tetrahedra (`fieldcraft_fem::uniform` holds their kinematics), rod2 and
   rod3 trusses (`discrete`). Point masses add to the lumped mass; springs and dashpots (on the
   previous half increment's velocity) add internal and viscous forces and their own stable
   increment. Other element kinds and formulations, couplings, reference points, rotary inertia,
   rotational springs and moments return "not supported yet". Masses are row-sum lumped. A zero-valued boundary condition holds a degree of freedom; a
   displacement value is reached at constant velocity over the step. Loads act in full from the
   start. The step lists initial velocities and fixed rigid planes (with friction) by set.
2. **Time step.** Each element's highest frequency is bounded in closed form (exact for a cube);
   the increment is the step's scale (0.9 by default) times the smallest element increment with
   Wilkins' bulk-viscosity correction, recomputed every increment. The ten controlling elements
   are reported; a step that would need more increments than its cap is refused before it runs.
3. **Integration.** Central differences in the half-step form, elements in label order, no
   parallel reductions (bit-for-bit repeatable). Rigid walls and prescribed degrees of freedom are
   kinematic velocity constraints whose impulses are the reactions.
4. **Energy account.** Kinetic, internal (stress work), hourglass, bulk viscosity, external and
   constraint work, and their balance; a run warns when the balance is off by more than 1% of the
   largest energy or the hourglass energy exceeds 10% of the internal energy.
5. **Results.** Histories (energies, momentum, constraint force, chosen nodes) at the history
   interval and field frames (displacement, velocity, element stress and strain) at the output
   interval, capped in size; kept in the session like the implicit results.

6. **Large deformation and materials (X2).** With `nonlinear_geometry` (the default) the
   elements are updated Lagrangian: strain and spin increments on the midpoint configuration,
   forces on the current one (hex8 gradients in closed form, `UniformStrain::hex8_deformed`). The
   material state (stress, back stress, strain, stiffness hourglass forces) is rotated with the
   Jaumann rate (Hughes–Winget) or the Green–Naghdi rate (Flanagan–Taylor), `fem::rotation`;
   materials are hypoelastic or J2 plastic (bilinear isotropic/kinematic, hardening table,
   Johnson–Cook with adiabatic heating) by radial return with Johnson–Cook damage,
   `fem::plasticity`. Elements whose damage reaches 1 or whose plastic strain reaches the
   material's erosion strain are eroded (reported with their energy; their nodes keep mass and
   momentum). The energy account adds the plastic work and the eroded energy.
7. **Contact (X3).** The model's contacts (`interaction.contact.*`) become interfaces before the
   elements are built (nodes that start inside a surface are moved onto it or refused):
   segment-to-segment penalty contact both ways (projection, polygon clipping and quadrature over
   each pair of facing faces), a bucket-grid search every increment in parallel by faces, a
   penalty per area from the masses and the stable increment that never shortens it, and Coulomb
   friction by a return map at the secondary nodes; rigid walls take Coulomb impulses. The contact
   forces join the equations of motion beside the external ones, their work and the friction
   dissipation join the energy account, and each contact's peak force, first and last contact,
   deepest penetration and friction are reported (docs/fem.md).

Verification (`crates/explicit/tests/verification.rs` and `large_deformation.rs`, analytic
only): stress waves in hex8 and tet4 bars (front at c·t, σ = ρ c v, work input), a bar striking a
rigid wall (contact time 2L/c, rebound speed, plateau force ρ c v₀ A), the energy balance,
bit-for-bit repeatability; Dienes' simple shear for both rates, a stressed element turned rigidly
by 90°, the J2 uniaxial curve and unloading at large strain, damage failure at the fracture
strain, and the Taylor bar against Taylor's 1948 estimate; `crates/explicit/tests/one_d.rs`:
spring–mass period, damped decay rate and period, a dashpot alone, two masses on a line spring,
waves in rod2 and rod3 bars; `crates/explicit/tests/contact.rs`: colliding bars, sliding
blocks under constant and decaying friction on a rigid wall and on a plate, the contact patch
test, Hertz's sphere impact, a free body untouched by general contact, initial penetrations and
thread-count independence; in the engine (`solve.rs`), a
cantilever's free-vibration period against `solve.frequency`.

## Lattice Boltzmann solver

`fieldcraft-lbm` solves a model's fluid steps (`solve.fluid`; `solve.probe` reads velocity,
pressure and temperature at points and wall force and heat histories). Papers are cited in the
module docs and in `docs/solvers/literature.md` (Part I).

1. **Model.** A fluid material (density, Newtonian viscosity; conductivity and specific heat
   for heat), a fluid domain (box, cubic cell spacing, a condition on each box face: wall,
   moving wall, slip, velocity inlet with a uniform or parabolic profile, pressure outlet,
   periodic pair; obstacle bodies and moving-wall conditions on their named faces; initial
   velocity; body acceleration; optional `thermal` settings: the initial temperature, a thermal
   condition per box face and per named obstacle face — held temperature, heat flux or
   adiabatic — gravity with the expansion coefficient for buoyancy, and the obstacles that
   conduct heat, each with its material and power) and a fluid step (end
   time or steps, optional time step, frames, the collision operator, the TRT parameter Λ).
   Deck writers refuse fluid steps.
2. **Setup.** Lattice units from the spacing, the time step and the density. Without a given
   time step the solver keeps the characteristic lattice velocity at 0.05 and τ ≤ 1. A lattice
   Mach number above 0.3, τ below 0.505 or a lattice acceleration above 0.01 is refused with the
   numbers and what to change; caps on cells, steps, cell updates and frame memory are checked
   before anything is allocated. With heat, the scalar relaxation time follows from the
   diffusivity k / (ρ c) in lattice units (refused below 0.505, warned above 2) and a flux on an
   inlet, outlet or periodic face is refused.
3. **Grid.** The engine triangulates each obstacle body (or takes an imported surface) and
   prepares its `BoundarySurface`; cells whose centres have a positive winding number are solid,
   decided with exact `orient2d`/`orient3d` and one symbolic perturbation, so faces through cell
   centres and lines through edges are classified consistently.
4. **Stepping.** TRT collision (or, per step, the recursive regularized collision that stays
   stable as τ approaches ½, `docs/lbm-collision.md`) with Guo forcing in place; pull streaming as one shifted copy per
   direction with per-link rules for the cells next to walls and box faces; Zou–He on inlet and
   outlet layers, or on a convective outlet the unknown populations updated by a convective
   outflow condition and each cell held at the face's pressure (Zou–He where the flow comes back
   in); a cell on the edge of two open faces belongs to one of them (velocity faces first).
   Wall links give the momentum exchanged with every box face and B-rep face (the face of the
   triangle nearest the link). A non-positive density or a lattice speed above
   0.5 stops the run with an error. With heat, a D3Q7 scalar lattice (weights ¼ and ⅛,
   c_s² = ¼) steps after the flow on the same links: it is advected by the cell velocity and
   diffuses with the antisymmetric relaxation time of its two-relaxation-time collision, whose
   Λ_T is fixed so that its walls and steady temperatures do not depend on that time (#629); held temperatures are anti-bounce-back, adiabatic
   faces bounce-back, a flux is added to the bounced population; velocity inlets carry the
   held temperature in at the cell's velocity, pressure outlets copy their neighbour. The heat
   crossing every wall link is summed per face as the flow's force is. The temperature feeds
   back through a per-cell Boussinesq force −β (T − T₀) g. Conducting obstacles carry their own
   temperature on their cells (`conduction.rs`): finite-volume conduction with harmonic-mean
   face conductivities per axis, explicit and sub-cycled, their power spread over their cells,
   the heat they give a held or flux box face counted as that face's;
   each D3Q7 link from the fluid into one is an anti-bounce-back at the interface temperature of
   the two half cells in series, and the heat it takes from the fluid goes to the solid cell,
   so heat is conserved to round-off; that heat is affine in the solid's temperature and taken
   point-implicitly, which keeps the coupling stable for any heat capacity. A run to steady state scales the solids' heat capacity
   down to the smallest the lattice allows (the steady temperatures do not depend on it). The
   GPU backend steps the conducting cells too (single precision, checked against the CPU).
   Fan and porous zones also run on the device: the same implicit half-force velocity,
   curve feedback before each collision, and optional zone flow/pressure histories. Zone
   headers and curves share the flags buffer, keeping the existing storage binding count
   ([implementation and measured cost](lbm-gpu-performance.md#fans-and-porous-resistance)).
5. **Results.** Frames of velocity, gauge pressure and (with heat) temperature per cell (single
   precision), force and heat histories per wall and the fluid mass, in model units, stored
   per step like the implicit results. Each conducting solid's peak temperature is read from
   the frames (the last one when the run reached steady state, else the hottest of them) and
   checked against its optional temperature limit (`FluidResult::components`).

The regularized collision has its own analytic and benchmark checks
(`crates/lbm/tests/regularized.rs`, and both collisions in the porous, fan, natural convection,
lid-driven cavity and Schäfer–Turek cases; `docs/lbm-collision.md`).

Verification (`crates/lbm/tests/verification.rs`, analytic): plane Poiseuille exact to round-off
with Λ = 3/16 and second order with Λ = 1/4; Couette exact with the wall shear μU/H; Taylor–Green
decay second order; Poiseuille in a staircase pipe; Stokes drag of a simple cubic sphere array
against Hasimoto (1959); a parabolic-inlet, pressure-outlet channel against Poiseuille; mass
conserved to round-off; bit-for-bit reproducibility. `crates/lbm/tests/thermal.rs` (heat):
linear conduction between held plates to round-off; a flux wall against q H / k and the wall's
heat balance; an inlet temperature carried through a plug flow; a periodic blob that keeps its
total and returns after one lap; (ignored, minutes) the de Vahl Davis cavity at Ra 10³ and 10⁴
against the published Nusselt numbers (both within 0.1 %). `crates/lbm/tests/conjugate.rs`
(conjugate heat transfer): a composite wall in series to round-off (orthotropic too); a
heat-generating slab at second order with its whole power reaching the plates; a heated wall
under a channel flow against the conduction profile; a heated block in a channel conserving
heat every step to round-off; a slow solid scaled for steady state landing on a long run's
temperatures; the heat-generating slab's peak passing and failing limits either side of it; House, Beckermann and Smith's cavity with a centred conducting body (conduction
limit every run; Ra 10⁵, ignored: +0.5 % on 48 cells, converging at second order to −0.22 %
extrapolated); Rayleigh's square array of cylinders against the Perrins series (conducting
rods within 0.3 %, insulating ones a first-order staircase error, −0.85 % on 64 cells); a
low-capacity solid with strong links stable at τ_T = 0.52; a two-layer wall with a thermal contact
dropping q R across it to round-off, and the heat across a contact conserved every step. `crates/lbm/tests/benchmark.rs` (ignored, published
reference values): Schäfer–Turek 2D-1 and 2D-2, the lid-driven cavity at Re 100 (Ghia et al. 1982).

## Commands

Every user-visible action is a `CommandSpec` in `crates/engine/src/cmd/<area>.rs` with a JSON
form. `Session::execute` clones the model, runs the command under `catch_unwind`, and restores
the model on an error or a panic; commands that changed the model are journalled, so the
journal replays as a script. Each journalled command is also an undo step (the state before it):
`edit.undo` restores that state and moves the command off the journal onto the redo stack,
`edit.redo` puts it back, so the journal always rebuilds the current model. History is bounded by
step count and by an estimate of the memory the kept states hold (oldest steps go first).
Commands of kind `Session` (`file.*`, `edit.*`) are never journalled. `file.save` and `file.open`
write and read `.fieldcraft` files; the desktop app autosaves unsaved work for `file.recover`.

Each command lists what it cannot run without (`requires`: an open design, geometry, a mesh, a
material, a step, a result, …; `crates/engine/src/requirement.rs`), next to its id. A requirement
is a real prerequisite, never a workflow order. `Session::execute` checks them first and fails
without doing anything, naming each missing piece and the toolbar entry that provides it ("Needs
a material — Properties › Material"). The desktop UI reads the same check
(`Session::availability`, a few counts, evaluated every frame) to grey ribbon buttons, command
search results and Solve, with the same text as their tooltip.

## Test corpus

Meshers are measured on our own parts (`docs/adr/0001-libraries.md`, test geometry). Each part is
a SolidCraft command script in `tests/corpus/<name>/recipe.json` (the file `solidcraft-cli run`
takes, plus a `corpus` object with its categories and, when it has one, its closed-form volume)
with a one-line `README.md`; analytic shapes (sphere, torus, cylinder, thick-walled pipe, plate
with a hole) are generated in code. `fieldcraft-testkit` runs a recipe in SolidCraft's engine,
saves the design and opens it through `open_design`, the path a person's design takes, so no
geometry is stored. `cargo xtask corpus` (a release build of the testkit's `corpus` binary)
builds every part, tessellates it and runs each stage of `fieldcraft_testkit::runner`, writing
docs/meshing/corpus.md. A mesher joins by implementing `runner::Stage` (its columns, and a run
over the part's tessellated surfaces) and adding itself to `runner::stages()`.

## Verification oracle

`tests/oracle/<case>/` holds the verification and validation cases (docs/verification.md): our
own model of each problem as a command recipe (`recipe.json`, optionally a SolidCraft
`design.json`), its cited reference values with tolerances and reasons (`expected.json`), and our
frozen results (`fixture.json`). `cargo xtask oracle` runs the testkit's `oracle` binary
(`fieldcraft_testkit::oracle`): one `Session` per case, quantities read by JSON pointer from step
outputs or commands such as `solve.probe {"at": [[x, y, z]]}`, compared with the reference and the
fixture; it writes docs/oracle.md. CI runs it as its own job.
