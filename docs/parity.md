# Capability checklist (reference: Abaqus/CAE)

What an analyst does in a full CAE pre-processor, grouped by milestone, with FieldCraft's status.
Scope comes from the public documentation of Abaqus/CAE and public analyst discussion,
described in our own words; nothing was taken from the product itself (AGENTS.md: no commercial
products are run, no documentation text is copied). Status: **done**, **started**, **planned**.
**HEX** marks items centred on hexahedral meshing, the core goal.

## Where we aim to be better

| # | Opportunity | Why analysts want it | FieldCraft approach |
|---|---|---|---|
| 1 | Hex meshes without manual partitioning | Partitioning a part into sweepable cells by hand is the most cited chore | Sweep and mapped meshing with automatic decomposition that proposes SolidCraft splits (M5); octree all-hex for anything else (M6) |
| 2 | Undo, redo, replay, crash recovery | No general undo; replay files are hard to reuse; corrupted model files | Every action is a journalled command with rollback, undo and redo; autosave with recovery after a crash |
| 3 | Sets that survive changes | Sets and loads break when geometry regenerates or is remeshed | Sets bind to SolidCraft's persistent face, edge and body names and resolve after meshing (M3) |
| 4 | An API that is not a bundled scripting runtime | Version breakage, slow scripts, geometry picked by position | JSON commands over CLI, control channel and MCP |
| 5 | Located diagnostics | Generic error messages | Errors name the element, face or keyword; `mesh.quality` lists failing elements |
| 6 | Lossless deck round trip | Unsupported keywords dropped on import | Deck readers keep unknown keywords and write them back (M9) |
| 7 | Safe, diffable model files | Binary model databases that corrupt or lock to a version | `.fieldcraft`: text, line-diffable definition and journal, checksummed binary mesh, versioned format tag, atomic saves with backup, autosave (M3, done) |
| 8 | Geometry checked before meshing | Slivers and tiny faces fail late | Surface preparation reports defects with their location (M1) |
| 9 | No size caps, any solver | Node limits in learning editions; one solver per tool | Free, open source, three deck formats |
| 10 | Clear step history | Loads and BCs propagate between steps implicitly | An explicit step × load/BC table (M3) |

## Checklist

| Area | Capability | Milestone | HEX | Status |
|---|---|---|---|---|
| Geometry | Open SolidCraft designs, bodies, persistent face names, assemblies placed | M1 | | done |
| Geometry | Face-tagged mesh import (3MF) | M1 | | done |
| Geometry | Weld, closed/manifold/orientation and self-intersection checks, feature edges | M1 | yes | started |
| Geometry | Partitions and datums as SolidCraft features (split body/face, construction planes) | M5 | yes | planned |
| Geometry | Midsurface extraction for shells | M8 | | planned |
| Mesh | Size field: global, per face and edge seeds, bias, curvature and proximity refinement | M2 | yes | planned |
| Mesh | Free tet meshing (tet4, tet10 with mid-side nodes on the geometry) | M2 | | started |
| Mesh | Quality metrics, certified validity, verify per set | M1–M2 | yes | done (metrics, gate) |
| Mesh | 2D quad meshing of faces | M5 | yes | started (`fieldcraft_mesh::quad`; `mesh.planar`: 2D models of a polygon or a planar face, quad4/quad8/tri3/tri6, edges tagged by their neighbouring face) |
| Mesh | Mapped hex of six-sided regions | M5 | yes | done (`mesh.hex.mapped`) |
| Mesh | Sweep hex (extrusions, revolutions) | M5 | yes | done (`mesh.hex.sweep`: multi-level extrusions, revolutions, one-to-one sweeps) |
| Mesh | Automatic decomposition into sweepable and mappable regions | M5 | yes | planned |
| Mesh | Conforming interfaces, tet–hex transition with pyramids | M5 | yes | planned |
| Mesh | Octree/grid-dual all-hex for any part | M6 | yes | planned |
| Mesh | Hex-dominant by recombination; polycube | M6 | yes | planned |
| Mesh | Smoothing, untangling, boundary padding, layers through thin walls | M6 | yes | planned |
| Mesh | Mesh editing (move, merge, split, collapse, remesh region), imported orphan meshes | M6 | yes | started (node/element add/remove, merge coincident nodes, duplicate and crack check, renumber) |
| Model | Sets and surfaces bound to geometry names | M3 | yes | done |
| Model | Materials: elastic (isotropic, orthotropic), density, unit-system check | M3 | | done |
| Model | Materials: plasticity, hyperelasticity, damage, thermal | M8 | | started (von Mises plasticity: bilinear, table, Johnson–Cook; Johnson–Cook damage, erosion; plasticity in explicit and nonlinear static steps; written to Abaqus `*PLASTIC`, Nastran `MATS1`, LS-DYNA plastic materials where expressible; compressible neo-Hookean, Mooney–Rivlin and Yeoh for nonlinear steps, Abaqus `*HYPERELASTIC, NEO HOOKE`, `MOONEY-RIVLIN`, `YEOH`; thermal strain in nonlinear steps) |
| Model | Sections: solid, shell, beam profiles; assignment by set | M3 | | done |
| Model | 2D sections: plane stress, plane strain, axisymmetric (XY plane, Y the axis); element edges in sets; centrifugal load | M4 | | done |
| Model | Composite layups, material orientations, continuum shells | M8 | yes | planned |
| Model | Element formulation per region (integration, incompatible modes, hybrid) mapped per solver | M3 | yes | done (model; per-solver mapping with the deck writers) |
| Model | Boundary conditions: fixed, displacement, symmetry | M3 | | done |
| Model | Loads: force and moment, pressure, gravity, body force | M3 | | done |
| Model | Reference points, kinematic and distributing coupling | M3 | | done |
| Model | Steps: linear static, frequency; output presets | M3 | | done |
| Model | Steps: nonlinear static, explicit dynamics, heat transfer; amplitudes | M7–M8, M15 | | started (nonlinear static with or without nonlinear geometry, increments, tolerances and a load path (cyclic loading), written to Abaqus `*STEP, NLGEOM` and LS-DYNA implicit (the load path as LS-DYNA's curve; refused by Abaqus for now), refused by Nastran; explicit dynamics; heat transfer; amplitudes per load planned) |
| Interaction | Tie between non-matching meshes | M7 | | done (`interaction.tie`: dual mortar constraints in the implicit analyses, the tie patch test; Abaqus `*TIE`, LS-DYNA tied contact) |
| Interaction | Contact (general and pairs), friction, automatic pair detection | M7, M10 | | started (explicit: pairs, general and self contact, segment to segment with a penalty that keeps the time increment, Coulomb friction with exponential decay, friction on rigid walls, `interaction.contact.find`, Abaqus and LS-DYNA contact keywords; implicit: mortar contact pairs in nonlinear static steps by an augmented Lagrangian with Coulomb friction, self contact next) |
| Interaction | Springs and dashpots (ground or two-point; DOF, vector or line), point masses with rotary inertia | M4 | | done (implicit; explicit translational; all three decks with what each cannot take reported) |
| Interaction | Bolt pretension, connectors | M7 | | planned |
| Decks | Abaqus input writer (`deck.write`, flat deck; reports what Abaqus cannot express) | M3 | | done |
| Decks | Nastran bulk data writer (`deck.write`, SOL 101 or 103 per deck; reports what Nastran cannot express) | M3 | | done |
| Decks | LS-DYNA keyword writer (implicit; explicit in M8) | M3 | | done (implicit static and eigenvalue; beam sections planned) |
| Decks | Readers keeping unknown keywords | M9 | | planned |
| Solve | Own linear static and modal solver, verified against analytic solutions | M4 | yes | started (solid elements: tet4, tet10, wedge6, hex8, hex20; 2D plane stress, plane strain and axisymmetric tri3, tri6, quad4 (full, incompatible modes), quad8 (full, reduced); quad4/tri3 MITC shells and Timoshenko beams with profiles, mixed with solids; reference points, kinematic and distributing couplings, moments; centrifugal loads) |
| Solve | Own nonlinear implicit solver: Newton–Raphson with line search and automatic increments, total Lagrangian solids and 2D elements, mixed u–p at finite strain, follower pressure | M15 | yes | started (solids, plane strain, axisymmetric; St Venant–Kirchhoff, neo-Hookean and J2 plasticity (Simo's exponential return map at finite strain, radial return at small strain, consistent tangents, the return shared with explicit dynamics); verified by the elastica, tube inflation, patch and objectivity tests, uniaxial tension to large strain, Hill's thick cylinder, a cantilever's collapse and the Bauschinger effect; mortar contact with friction (patch test, Hertz, sliding blocks); Mooney–Rivlin and Yeoh; thermal strain at finite and small strain; a symmetric indefinite LDLᵀ factorization with inertia, Newton through indefinite tangents (hybrid elements at ν = 0.4999 under prescribed displacement); linear buckling of solids (Euler columns, Bryan's plate) and imperfections from its modes; arc-length with critical points from the inertia (shallow-arch snap-through); shells next) |
| Solve | Own explicit dynamics solver (central differences, one-point hex8 with hourglass control, nodal-pressure tet4, bulk viscosity, energy balance), verified against analytic wave and impact solutions | M10 | yes | started (X1 core and X2: large deformation, objective rates, J2 and Johnson–Cook plasticity, damage, erosion; X3: deformable contact and friction, verified by bar collision, sliding blocks, the contact patch test and Hertz; hyperelastic materials on the shared material point; opt-in mass scaling) |
| Results | Deformed shape, contours, mode shapes, probes, reaction sums | M4 | | started (results in the session, `solve.probe`, reaction sums; display planned) |
| Results | XY and path plots, cuts, display groups, free-body forces, animation | M9 | | started ([solid-field path plots](path-plots.md): typed polyline query, genuine surface picking, reference/deformed overlay, gaps and CSV; remaining postprocessing scope is separate) |
| Results | Readers for public result formats (Nastran OP2, LS-DYNA d3plot, Abaqus `.fil`) | M9 | | planned |
| Results | Mesh convergence studies | M9 | | planned |
| App | Desktop viewport, model tree, commands, control channel | continuous | | done |
| App | Undo/redo on the command journal, autosave, `.fieldcraft` file | M3 | | done |
| App | MCP server, browser build | continuous | | planned |

## Deck notes for the writers

- Node order: the model stores Abaqus/VTK order. Nastran orders the mid-side nodes of 20-node hexes
  and 15-node wedges differently; LS-DYNA takes tets and wedges as degenerate hexes and triangles
  as degenerate quads. Every mapping is tested by the sign of the element Jacobian after writing.
- Each format's field widths, continuation rules, number forms and label limits are respected by
  the writer and checked by golden-file tests on our own examples.
- Abaqus: the model's node order and face numbers are Abaqus's, so connectivity passes through;
  element types come from the element kind, the section and its formulation (table in
  `crates/decks/src/abaqus.rs`). The step table is bridged to Abaqus's propagation: zero-valued
  boundary conditions active in every step go to the initial state, and steps restate their
  boundary conditions and loads with `OP=NEW` when the active set changes.
- Nastran (`fieldcraft_decks::nastran`): entries with a real use large (16-column) fields, the rest
  small fields; a real is exact when its shortest round-trip form fits 16 columns, else rounded to
  at least ten significant digits. Linear solids map full integration to `PSOLID` `TWO`/`THREE` +
  `FULL`, hex20 reduced to `TWO` + `FULL`, hex8 incompatible modes to `BUBBLE`; one-point reduced
  linear solids, hybrid solids, reduced-integration or one-point shells, body forces, gravity on part
  of the model, orthotropic beams and rods, and orthotropic shells off planes of constant Z are
  reported as unsupported. Static and frequency steps need separate decks (`deck.write` `steps`).
- LS-DYNA (`fieldcraft_decks::lsdyna`, `deck.write {"format": "lsdyna"}`): long format (20-column
  fields; every value with up to 13 significant digits, and 17-digit values from 0.01 to 1e17, are
  written exactly; the rest is rounded to at least 13 digits and counted). The 10-node tetrahedron
  follows the manual's figure: mid-side nodes 8, 9, 10 on edges 2-4, 3-4, 1-4 (the model has them
  on 1-4, 2-4, 3-4). Incompatible modes map to ELFORM −18 (eight-point enhanced strain). Several
  steps are written as `*CASE` analyses. Refused with the item named: beam sections (not yet),
  body forces per volume, partial-DOF kinematic couplings, constraints on a coupled reference
  point, gravity on different sets or on parts of sections, hybrid other than tet4, reduced
  quadratic elements, materials without density, more than 10 shell integration points.
- 2D models: Abaqus writes CPS*/CPE*/CAX* (3, 6, 4, 4R, 4I, 8, 8R) with the thickness on the
  `*SOLID SECTION` line, edge surfaces `S1`… and `BR`/`BZ` body forces for axisymmetric elements;
  Nastran writes plane stress as membranes (`PSHELL` with MID1 only, grids holding 3456 through
  `GRID` PS) and edge pressures as consistent `FORCE`s, and refuses plane strain and
  axisymmetric elements (their cards differ between Nastran dialects); LS-DYNA writes plane stress
  and plane strain (unit thickness) as `*SECTION_SHELL` ELFORM 12/13, one-point quad4 (reduced
  integration, with hourglass control) and tri3, and refuses quadratic 2D elements, full or
  incompatible-mode quad4 and axisymmetric sections. Centrifugal loads: Abaqus `CENTRIF`
  (ω², point, axis), Nastran `RFORCE` about an axis through the origin, LS-DYNA not yet.
- Open item (LS-DYNA tet10 node order): the writer follows the manual's figure (mid-side nodes 8,
  9, 10 on edges 2-4, 3-4, 1-4), which disagrees with an earlier reading as the Abaqus order. It
  needs one confirmation run by a licensed user; agents never run commercial solvers.
