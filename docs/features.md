# Features

Generated from the capability catalogue (`crates/engine/src/catalogue.rs`) by
`cargo run -p fieldcraft-cli -- features > docs/features.md`; a test fails when it is out of date.
Every entry is a toolbar button. **Ready** entries run their command; **coming soon** entries are
shown greyed with the roadmap milestone that brings them (ROADMAP.md). Scripts, the command line,
the control channel and AI agents over MCP (docs/mcp.md) run the same commands.

87 ready, 35 coming soon.

## Geometry

| Panel | Capability | Status | Command |
|---|---|---|---|
| Import | Open Design — Open a SolidCraft design: its bodies with persistent face names | ready | `geometry.open` |
| Import | Tagged Mesh — Read a face-tagged surface mesh (3MF) | ready | `geometry.import.mesh` |
| Inspect | Bodies — List the open design's bodies, volumes and named faces, edges and vertices | ready | `geometry.bodies`, `geometry.centre`, `geometry.features` |
| Inspect | Tessellate — Tessellate a body into a face-tagged boundary surface | ready | `geometry.tessellate` |
| Inspect | Check Geometry — Short edges, small and sliver faces and near faces against the element size, located | ready | `geometry.check` |
| Inspect | Check Surface — Weld, closed, manifold, orientation and self-intersection checks; feature edges | ready | `mesh.surface.check` |
| Prepare | Datum Plane — A datum plane: principal plane and offset, three points or datum points, point and normal, normal to an axis | ready | `geometry.datum`, (`geometry.datum.edit`), (`geometry.datum.delete`), (`geometry.datum.suppress`) |
| Prepare | Datum Axis — A datum axis: principal axis, two points or datum points, point and direction | ready | `geometry.datum.axis`, (`geometry.datum.axis.edit`) |
| Prepare | Datum Point — A datum point: coordinates, a mesh vertex, along an edge, a face's centre, an axis through a plane | ready | `geometry.datum.point`, (`geometry.datum.point.edit`) |
| Prepare | Partition — Cut bodies into cells by a datum plane; the two faces of each cut are paired by name | ready | `geometry.partition`, (`geometry.partition.delete`), (`geometry.partition.suppress`) |
| Prepare | Merge Faces — Mesh adjacent faces as one (virtual topology): a fillet in pieces, a face split by an imprint; the B-rep is unchanged | ready | `geometry.faces.merge`, (`geometry.faces.merge.delete`) |
| Prepare | Absorb Faces — Absorb tiny and sliver faces into a neighbour for meshing, one picked or every one Check Geometry found | ready | `geometry.face.absorb` |
| Prepare | Split Face — Split a face along a datum plane for meshing only: two faces on its surface, the cut an edge between them | ready | `geometry.face.split`, (`geometry.face.split.delete`) |
| Prepare | Defeature — Take small holes, fillets, chamfers and bosses out for meshing, from a proposal list with each one on or off | ready | `geometry.defeature`, (`geometry.defeature.propose`), (`geometry.defeature.delete`) |
| Prepare | Stitch — Close an imported surface's gaps under a tolerance: open points welded, T-junctions split | ready | `geometry.stitch`, (`geometry.stitch.delete`) |
| Prepare | Midsurface — Extract midsurfaces of thin parts for shell meshes | coming soon (M8) |  |

## Mesh

| Panel | Capability | Status | Command |
|---|---|---|---|
| Size | Global Size — One element size for the whole part | coming soon (M2) |  |
| Size | Seed Edges — Element counts or sizes per design edge, with bias and an exact, minimum or free constraint; the meshers honour them | ready | `mesh.seed.edge`, `mesh.seed.clear`, `mesh.seed.list` |
| Size | Local Size — Element size on faces or bodies, in a box or in a sphere, growing away at a set rate; every mesher honours it | ready | `mesh.size.set`, `mesh.size.growth`, `mesh.size.clear`, `mesh.size.list` |
| Size | Curvature — Refine by curvature and proximity | coming soon (M2) |  |
| Tet | Tet Mesh — Quality tet4 or tet10 mesh of a body to a size: refined, smoothed and improved, boundary on the B-rep faces, faces tagged | ready | `mesh.tet.generate` |
| Tet | Quadratic — Second-order elements with mid-side nodes on the geometry | coming soon (M2) |  |
| Tet | Delaunay Points — Delaunay tetrahedralization of a point list (debugging aid) | ready | `mesh.delaunay.points` |
| 2D | 2D Mesh — Quad or triangle mesh of a planar face or polygon in the XY plane for plane stress, plane strain and axisymmetric models; edges tagged | ready | `mesh.planar` |
| Structured Hex | Quad Faces — All-quad meshes of faces on shared edge nodes | coming soon (M5) |  |
| Structured Hex | Hex Mesh — Structured hex mesh of a body: automatic, or extruded, swept, revolved or mapped (the Method drop-down) | ready | `mesh.hex.generate`, (`mesh.hex.sweep`), (`mesh.hex.mapped`) |
| Structured Hex | Partitioned Mesh — One mesh of a partitioned part, a method per cell (hex, tet, or all hex as conforming blocks); cells joined across partitions, tet cells to hex cells through conforming pyramids | ready | `mesh.cells`, (`mesh.hex.blocks`) |
| Structured Hex | Extrude to Surface — Hex layers from a meshed or planar section up to a curved target surface, following the side faces | ready | `mesh.hex.extrude` |
| Structured Hex | Mesh Plan — A hex-meshing plan from the features that built a body (sweeps, blocks, O-grids, layers, imprints, copies); saved, editable, applied with the structured meshers | ready | `mesh.plan`, `mesh.plan.apply` |
| Structured Hex | Mesh Block — A mapped block of hex, wedge or quad elements between corner points; projection blocking (a block topology projected onto the part) | ready | `mesh.block`, `mesh.hex.project` |
| Structured Hex | Decompose — Find sweepable and mappable regions of 3D parts (extruded cross-sections are already cut into blocks by Hex Mesh) | coming soon (M5) |  |
| Automatic Hex | Octree Hex — Unstructured all-hex mesh of any part: adaptive grid, boundary layers on every sharp edge, untangled and proven valid | ready | `mesh.hex.octree` |
| Automatic Hex | Hex-Dominant — Hex-dominant mesh by recombining tets | coming soon (M6) |  |
| Automatic Hex | Polycube — Polycube hex mesh of blocky parts | coming soon (M6) |  |
| Improve | Smooth — Smooth node positions, keeping every element valid | coming soon (M6) |  |
| Improve | Untangle — Repair inverted elements | coming soon (M6) |  |
| Improve | Layers — Boundary padding and layers through thin walls | coming soon (M6) |  |
| Improve | Relax Nodes — Relax the nodes inside a sphere to fix spikes and slivers; surface nodes slide on the geometry; kept as a field to apply again after a remesh | ready | `mesh.node.relax`, `mesh.node.relax.apply` |
| Improve | Move Nodes — Drag nodes of the mesh, the nodes round them following by a falloff or harmonically, on the geometry; the original drawn as a ghost and the elements coloured by how their quality changes | ready | `mesh.node.move` |
| Edit | Node — Add a node | ready | `mesh.node.add`, `mesh.node.remove` |
| Edit | Element — Add an element on existing nodes | ready | `mesh.element.add`, `mesh.element.remove`, `mesh.clear` |
| Edit | Mirror Mesh — Mirror a body's mesh across a plane, joined on the plane (mesh half a symmetric part) | ready | `mesh.mirror` |
| Edit | Merge Nodes — Merge coincident nodes, keeping sets, loads and boundary conditions on the survivors; check for duplicates and cracks, renumber | ready | `mesh.node.merge`, `mesh.check.duplicates`, `mesh.renumber` |
| Edit | Remesh Region — Remesh a region in place; split and collapse elements | coming soon (M6) |  |
| Edit | Orphan Mesh — Import a mesh without geometry | coming soon (M6) |  |
| Verify | Quality — Quality metrics and the certified validity check, failing elements listed | ready | `mesh.quality` |
| Verify | Mark Bad Elements — Mark the elements below quality limits as a set, located and drawn in the marker colour | ready | `mesh.quality.mark` |
| Verify | Verify Set — Quality and validity per set | coming soon (M2) |  |
| Verify | Export STL — Write the mesh's outer surface (or a set's faces) as binary STL | ready | `mesh.export.stl` |
| Verify | Tag Geometry — Tag elements, faces, edges and nodes with the geometry they mesh | ready | `mesh.tag`, `mesh.tag.clear` |

## Properties

| Panel | Capability | Status | Command |
|---|---|---|---|
| Materials | Material — Elastic (isotropic, orthotropic) material with density; neo-Hookean hyperelasticity for nonlinear steps | ready | `model.material.create`, `model.material.edit`, `model.material.delete` |
| Materials | Plasticity — Von Mises plasticity: bilinear, hardening table or Johnson–Cook | ready | `model.material.plasticity` |
| Materials | Hyperelastic — Hyperelastic material models | coming soon (M8) |  |
| Materials | Damage — Johnson–Cook damage and element erosion | ready | `model.material.damage` |
| Materials | Thermal — Conductivity, specific heat and expansion | coming soon (M8) |  |
| Materials | Equation of State — Polynomial, ideal gas, Mie–Grüneisen, JWL and tabulated equations of state | coming soon (M11) |  |
| Sections | Section — Solid, shell, beam, rod or 2D (plane stress, plane strain, axisymmetric) section assigned to a set, with its element formulation | ready | `model.section.create`, `model.section.edit`, `model.section.delete` |
| Sections | Composite — Composite layups and continuum shells | coming soon (M8) |  |
| Sections | Orientation — Material orientations | coming soon (M8) |  |
| Sets | Set — A set of nodes, elements and faces, or of named geometry that survives remeshing | ready | `model.set.create`, `model.set.edit`, `model.set.delete`, `model.set.resolve` |
| Model | Units — The units the model is shown and typed in | ready | `model.units` |

## Assembly

| Panel | Capability | Status | Command |
|---|---|---|---|
| Instances | Instances — Part instances placed from SolidCraft occurrences | coming soon (M7) |  |

## Steps

| Panel | Capability | Status | Command |
|---|---|---|---|
| Steps | Step — A linear or nonlinear static, frequency, heat or explicit step with its output | ready | `model.step.create`, `model.step.edit`, `model.step.delete` |
| Steps | Activate — Choose which loads and boundary conditions act in a step | ready | `model.step.activate` |
| Steps | Step Table — The step by load and boundary condition table | ready | `model.activity` |
| Procedures | Amplitude — Time and frequency amplitudes | coming soon (M7) |  |
| Procedures | Nonlinear — Nonlinear static steps | coming soon (M7) |  |
| Procedures | Explicit — Explicit dynamic steps | coming soon (M8) |  |
| Procedures | Heat Transfer — Heat transfer steps | coming soon (M8) |  |

## Loads & BCs

| Panel | Capability | Status | Command |
|---|---|---|---|
| Boundary Conditions | Boundary Condition — Fixed, displacement or symmetry condition on a set or reference point | ready | `model.boundary_condition.create`, `model.boundary_condition.edit`, `model.boundary_condition.delete`, `model.boundary_condition.symmetry` |
| Loads | Load — Force, moment, pressure, traction, gravity or body force | ready | `model.load.create`, `model.load.edit`, `model.load.delete` |
| Flow | Fluid Domain — The box of fluid cells, its fluid, and the obstacle bodies voxelised from their face-tagged surfaces | ready | `model.fluid_domain.create`, `model.fluid_domain.edit`, `model.fluid_domain.delete` |
| Flow | Inlet & Outlet — Velocity inlets, pressure outlets, walls, moving and slip walls, periodic faces; held temperatures, heat fluxes and adiabatic faces; conducting solids' temperature limits | ready | `model.fluid_domain.boundary`, `model.fluid_domain.thermal_boundary`, `model.fluid_domain.limit` |

## Interactions

| Panel | Capability | Status | Command |
|---|---|---|---|
| Constraints | Reference Point — A reference point for couplings, loads and boundary conditions | ready | `model.reference_point.create`, `model.reference_point.edit`, `model.reference_point.delete` |
| Constraints | Coupling — Kinematic or distributing coupling of a set to a reference point | ready | `model.coupling.create`, `model.coupling.edit`, `model.coupling.delete` |
| Constraints | Tie — Tie non-matching meshes together (dual mortar constraints) | ready | `interaction.tie`, `interaction.tie.edit`, `interaction.tie.delete` |
| Contact | General Contact — Contact between everything that may touch, or self contact of a surface, with friction (explicit dynamics and nonlinear static steps) | ready | `interaction.contact.general`, `interaction.contact.edit`, `interaction.contact.delete` |
| Contact | Contact Pair — Contact between two surfaces: with friction in explicit dynamics, mortar contact with friction in nonlinear static steps | ready | `interaction.contact.pair` |
| Contact | Find Contacts — Find the parts of the mesh that touch, and make contact pairs of them | ready | `interaction.contact.find` |
| Contact | Rigid Wall — A fixed rigid plane a set's nodes cannot pass, with friction, in an explicit dynamics step | ready | `interaction.rigid_wall`, `interaction.rigid_wall.delete` |
| Connectors | Bolt Load — Bolt pretension | coming soon (M7) |  |
| Connectors | Spring — Spring or dashpot to the ground or between two points, along a direction or their line | ready | `model.spring.create`, `model.spring.edit`, `model.spring.delete` |
| Connectors | Point Mass — Point mass and rotary inertia at nodes or a reference point | ready | `model.mass.create`, `model.mass.edit`, `model.mass.delete` |
| Connectors | Connector — Connector elements between points | coming soon (M7) |  |
| Coupled | FSI Interface — The fluid–structure interface | coming soon (M13) |  |

## Solve

| Panel | Capability | Status | Command |
|---|---|---|---|
| Run | Solve — Solve every step in order and show the results | ready | `solve.run` |
| Check | Check Model — Name every missing or bad definition before writing or solving | ready | `model.check` |
| Decks | Abaqus Input — Write an Abaqus input file (.inp) | ready | `deck.write {"format":"abaqus"}` |
| Decks | Nastran Deck — Write Nastran bulk data (SOL 101, 103) | ready | `deck.write {"format":"nastran"}` |
| Decks | LS-DYNA Deck — Write an LS-DYNA keyword file | ready | `deck.write {"format":"lsdyna"}` |
| Decks | Read Deck — Read a deck, keeping keywords FieldCraft does not know | coming soon (M9) |  |
| Implicit | Static — Linear static (or coupled temperature–displacement) solve of a step with FieldCraft's own solver; nonlinear static, buckling and heat transfer in the drop-down | ready | `solve.static`, `solve.nonlinear`, `solve.buckle`, `solve.heat` |
| Implicit | Frequency — Natural frequencies and mode shapes of a frequency step | ready | `solve.frequency` |
| Explicit | Explicit — Transient dynamics of an explicit step by central differences (contact and erosion coming in M10) | ready | `solve.explicit` |
| Explicit | Blast & Impact — EOS-driven fluids with fluid–structure contact (ALE or SPH) | coming soon (M14) |  |
| Fluids | Lattice Boltzmann — Lattice Boltzmann flow of a fluid step on the CPU or a GPU, with wall forces, and the heat it carries | ready | `solve.fluid`, `solve.settings` |
| Fluids | Fluid–Structure — Coupled fluid and structure solve | coming soon (M13) |  |

## Results

| Panel | Capability | Status | Command |
|---|---|---|---|
| Display | Deformed — The deformed shape of a solved step, to scale | ready | `results.deformed`, `results.clear` |
| Display | Contours — Contours of displacement and stress, with a legend | ready | `results.contour` |
| Display | Mode Shapes — Mode shapes and frequencies of a frequency step | ready | `results.modes` |
| Display | Animate — Animate the deformed shape, a mode or an explicit step's frames | ready | `results.animate`, `results.frame` |
| Query | Probe — Displacement, rotation, reaction and stress at nodes, element section forces, a mode shape, or an explicit field frame (or click the model) | ready | `solve.probe`, `solve.history` |
| Query | Reactions — Reaction force and moment sums per set | ready | `results.reactions`, `results.export.csv` |
| Query | Free Body — Free-body forces across a cut | coming soon (M9) |  |
| Plots | XY Plot — Plot a quantity against the load factor, time or mode | ready | `results.xy` |
| Plots | Path Plot — Interpolate a solid nodal field along a reference-coordinate path, with CSV export | ready | `results.path` |
| Plots | Cut — Section cuts through the model | coming soon (M9) |  |
| Plots | Display Group — Show and hide parts of the model | coming soon (M9) |  |
| Plots | Streamlines — A fluid step's velocity, pressure or temperature on slices and obstacles, speed or temperature along streamlines, and signed wall heat per area | ready | `results.streamlines`, `results.slice`, `results.surface` |
| Files | Export VTU — Write a solved step as VTK unstructured grids for ParaView | ready | `results.export.vtu`, `results.export.animation`, `results.render.save`, `results.render.delete` |
| Files | Read Results — Read results in public formats | coming soon (M9) |  |
| Files | Convergence — Mesh convergence studies | coming soon (M9) |  |

## Tools

| Panel | Capability | Status | Command |
|---|---|---|---|
| Model | Model Info — Node and element counts and the model summary | ready | `model.info`, `model.rename` |
| File | Open Model — Open a .fieldcraft model | ready | `file.open`, `file.new` |
| File | Samples — Open a sample model, ready to solve | ready | `file.sample` |
| File | Save — Save the model to its .fieldcraft file | ready | `file.save` |
| File | Save As — Save the model to a new .fieldcraft file | ready | `file.save_as` |
| File | Recover — Open a model autosaved by an app that crashed or closed with unsaved changes | ready | `file.recover`, `file.recovery.list`, `file.recovery.discard` |
| View | Show All — Show every item hidden in the viewport (hide, show and isolate items from the browser) | ready | `view.show_all`, (`view.hide`), (`view.show`), (`view.isolate`), (`view.save`), (`view.restore`), (`view.delete`), (`view.list`) |
| History | Undo — Undo the last command that changed the model | ready | `edit.undo` |
| History | Redo — Redo the command undone last | ready | `edit.redo` |
| Help | Report a Bug — Draft a GitHub bug report with the version, the OS and the last error; you review it, and nothing is sent until you submit it | ready | `diagnostics.bug_report`, `diagnostics.crash_reports`, `diagnostics.crash_report`, `diagnostics.crash_reports.delete`, `diagnostics.crash_test`, `app.update.check` |
