# Examples

Every example is a built-in sample: open it from **File ▸ Samples** (or the start screen), with
`fieldcraft --sample NAME`, or with the `file.sample` command (`{"name": "drop-test"}`), then
press **Solve** (`solve.run`). Each opens fully defined — geometry, mesh, materials, boundary
conditions, loads and steps — in mm, t, s, N and MPa, and opens and solves within half a minute. The geometry is our
own: SolidCraft scripts run when the sample opens (a recipe from [`tests/corpus`](../tests/corpus)
or one written in `crates/engine/src/sample.rs`). Each sample has a test in
`crates/engine/src/sample.rs` or its modules that opens it, checks it clean, solves it and holds its result to
the check given below.

## Drop test (`drop-test`)

![The bracket at 0.1 ms after striking the tilted floor, von Mises stress up to 129 MPa in the bend](media/examples/drop-test.png)

An aluminium L-bracket (60 × 40 mm foot, 50 mm upright, 8 mm thick, R6 fillet) falls 2 m onto a
rigid floor tilted 15°, so the end of its foot strikes first. Explicit dynamics for 1 ms with
gravity: 1 911 swept hex8 elements (one-point, stiffness hourglass control), bilinear plasticity
(yield 240 MPa), the floor a rigid wall. About 8 400 increments, a few seconds.

**Look at:** the frames from the strike (Results ▸ Frame): the compression wave leaves the
struck corner and the upright's inertia bends the fillet, where von Mises peaks near 130 MPa —
below yield, so the drop leaves no permanent set. The energies plot shows the kinetic energy
turning into vibration with a balance under 0.01 %.

**Check:** the mass equals ρ V of the bracket (0.5 %); it meets the floor at √(2 g h) =
6.26 m/s, with ½ m v² = 1 748 mJ of kinetic energy (1 %); energy is conserved (balance < 1 %);
by 1 ms the floor has taken about half its downward momentum.

## Rubber seal (`rubber-seal`)

![The gasket ring squeezed 25 %, bulging outward, von Mises stress peaking at the bonded edges](media/examples/rubber-seal.png)

A rubber gasket ring (20 mm inside radius, 30 mm outside, 10 mm high) squeezed 2.5 mm (25 %)
between two plates bonded to it. Generic rubber, E = 6 MPa, ν = 0.49, Mooney–Rivlin with
C01 = 0.2 MPa; 1 000 hybrid (u–p) hex8 elements in four mapped quarter blocks. Nonlinear static
with large deformation: Newton–Raphson in nine increments, under a second.

**Look at:** the deformed shape at scale 1 — the free inner and outer walls bulge out, and
von Mises peaks (about 7 MPa) where the bulge meets the bonded plates. Results ▸ Increment steps
through the converged increments.

**Check:** the same ring between frictionless plates compresses evenly, and its force matches
the closed form for an incompressible Mooney–Rivlin solid in uniaxial compression,
P = 2 (λ − λ⁻²)(C10 + C01/λ) per original area: 3 430 N against 3 466 N (within 3 %). Bonded,
the bulge makes it stiffer (shape factor 0.5): 6 518 N, held between 1.05 and 2.5 times the
frictionless force. The groove and the plates are prescribed displacements today; with implicit
contact (#161) the ring will be pressed by contact against a groove.

## Flow past a bracket (`flow-past`)

Air (ρ = 1.2 kg/m³, μ = 1.8 × 10⁻⁵ Pa s) at 25 mm/s past the same L-bracket standing on the
floor, the upright's back to the wind. A 240 × 120 × 100 mm box of 2.5 mm cells (184 320),
uniform inlet, outlet at zero pressure, a no-slip floor and slip walls at the sides and top;
D3Q19 lattice Boltzmann (TRT) for 16 s, about 1.7 passes through the 240 mm box, ten frames. Reynolds number
83 on the upright's 50 mm height: laminar, a steady wake. About 10 s.

**Look at:** the step's report (`solve.fluid`, or Results after Solve): the force on every
named face of the bracket and on each side of the box, final and time-averaged, and the mass
history; `solve.probe` with `points` gives the velocity and pressure in any cell, and with
`forces` each face's force history.

**Check:** the drag on the bracket, as a coefficient on the upright's frontal area
(40 × 50 mm), is 2.3 — a bluff plate standing on a wall at low Reynolds number in a box it
blocks by 17 %, held between 1 and 3.5 (a flat plate broadside to the flow is about 1.2 at high
Reynolds number and rises as it falls); the lattice's mass changes by less than 1 %.

![The L-bracket the air flows past, as the sample opens](media/examples/flow-past.png)

The picture shows the bracket as the sample opens. After Solve the viewport draws the flow:
velocity magnitude on a slice through the middle of the box (`results.slice`; pressure with
`field: "pressure"`), the bracket coloured by the flow beside it (`results.surface`) and
streamlines seeded upstream of it, coloured by speed (`results.streamlines`).

![The solved flow: velocity on the mid-slice, a slow wake behind the upright and faster flow over it, with streamlines bending around the bracket](media/examples/flow-past-flow.png)

**Watch it start up:** a fluid step animates its frames (Results ▸ Animate, or `results.animate`),
playing the start-up transient from the uniform 25 mm/s initial field to the settled wake — the
wake forming behind the upright and the drag settling over the 16 s run. A clip is at
[flow-past-transient.mp4](media/examples/flow-past-transient.mp4), regenerated by
[tools/flow-clip](../tools/flow-clip/README.md).

## Heated bracket (`heated-bracket`)

![The steel bracket's temperature, 120 °C at the heated end of the upright falling to 95 °C at the foot](media/examples/heated-bracket.png)

The corpus filleted bracket (`tests/corpus/filleted-bracket`: an 8 mm L-bracket, R6 fillet, a
hole through each leg) in steel, both hole bores clamped. The end of the upright is held at
120 °C by a heater; the rest of its outside loses heat to still air at 20 °C
(h = 10 W/(m² K)); steel conducts with k = 50 W/(m K), expands with α = 12 × 10⁻⁶/K and is
stress free at 20 °C. One coupled temperature–displacement step (steady) on about 23 000 tet10
elements (the count differs a little between platforms): about 5 s to open (meshing) and 13 s to
solve.

**Look at:** Contours ▸ Temperature: the upright stays within a few degrees of the heater and
the foot falls to 95 °C, a steady balance between conduction and the air. Contours ▸ Von Mises
shows the stress of the growth the clamps hold back: tens of MPa through the legs, peaking at the
edges of the clamped bores, where a rigid clamp is a stress singularity (about 1.5 GPa there, a
figure that grows as the mesh is refined: read the stress a little away from the bore).

**Check:** the temperatures lie between the air's and the heater's and fall away from the heater.
For the expansion, the test puts the air at the heater's temperature, so the bracket warms
evenly by ΔT = 100 K, and clamps the foot's hole alone: the bracket then grows freely about it,
and a point moves α ΔT times its distance from the hole's axis. The upright's top moves −0.0432
mm along x (within 2 %); the foot's end, 16 mm from the clamped bore, which cannot grow, moves
0.0228 mm against 0.024 mm (within 6 %).

## Hex-meshed bracket (`bracket-hex`)

![Von Mises stress on the all-hex bracket, O-grids round both holes, peaking at the upright's hole](media/examples/bracket-hex.png)

The same filleted bracket split into 30 cells in SolidCraft
(`tests/corpus/filleted-bracket-blocks`: planar cuts through the legs and diagonal cuts round
the holes and the fillet, each cell a six-sided block), meshed by `mesh.hex.blocks` as one
conforming all-hex mesh: 5 778 hex20 elements of 3 mm, O-grids round the holes, graded through
the fillet, every shared face meshed once. Steel; the upright's end clamped, 1 kN down on the
foot's end. Linear static: about 2 s to open (the cells are rebuilt from the script) and 3 s to
solve.

**Look at:** the mesh (Prepare ▸ Mesh) — mapped blocks meeting face to face with no pyramids or
tets — and the von Mises contours (tick Deformed for the deformed shape): the stress runs round the fillet and
peaks at 295 MPa at the upright's hole.

**Check:** the tip of the foot moves 0.509 mm down, within 1 % of a 48 016-element tet10
solution of the unpartitioned part (0.5093 mm, `crates/testkit/tests/bracket_blocks.rs`), and the
clamp carries the whole 1 kN.

## Mirrored crank arm (`crank-mirror`)

![Von Mises stress on the crank arm from below: the pocket ringed by mapped rows, the half and its mirror image one mesh](media/examples/crank-mirror.png)

The model of `crank-arm` under the name of its mirror workflow (the sample kept so scripts and
links made with it still open): seen from below, where the pocket is. How the half is cut,
meshed and mirrored is told under `crank-arm`.

**Look at:** the pocket from below: straight mapped rows along its sides and floor, and at each
end its corner fillets as three-sided triangle primitives, cut off by the mesher's own blocks;
the two halves, mirror images across the arm's middle.

**Check:** as `crank-arm`.

## Bicycle crank arm (`crank-arm`)

![Von Mises stress on the structured all-hex crank arm, the arm root and the square bore carrying the most](media/examples/crank-arm.png)

A 170 mm aluminium crank built in SolidCraft (`crates/engine/samples/crank-arm.json`): a 16 mm
arm tapering from 38 to 24 mm, root arcs into a Ø46 × 22 spindle boss with a 12.7 mm square bore
(R1.5 corners) and a Ø30 × 20 pedal boss with a Ø14 hole, and an 82 × 14 × 8 pocket in its back
(R5 corners). It is a stack of extrusions along z, so it sweeps all-hex. To keep it sweepable the
hex variant leaves out the pocket's floor fillet and the bosses' 1.5 mm chamfers, and its root
arcs meet the bosses about 20° off tangent (a tangent meeting leaves a zero-angle cusp on the
arm's top face, where any hex is degenerate). 6061-T6 (E = 68 900 MPa, ν = 0.33); the square
bore's faces fixed, 1500 N on the pedal hole in the crank's plane, across the arm (−y). Linear
static: about 5 s to open and 7 s to solve.

The mesh needs one datum. A datum plane on the symmetry plane (`Sym`, y = 0) partitions the
crank (`Half`); `mesh.hex.generate` extrudes the half `CrankArm/a` (hex20 of 3 mm) and cuts its
cross-section into mappable blocks itself, inside the mesher (nothing is added to the model):
the half ring round each bore as four-sided blocks, the pocket's R5 corners as three-sided
triangle primitives, and four-sided blocks along the arm, split where the pocket's fillets meet
it. Interval assignment gives every block's sides matching counts, with no seeds.
`mesh.mirror` reflects the half onto `CrankArm/b`: 8 664 hex20 elements, all valid, scaled
Jacobian min 0.37 (at the corners where the root arcs meet the bosses 20° off tangent), mean
0.91. Before the mesher cut blocks itself the same crank took seven partitions and eleven exact
seeds for 0.37 (and the swept half without them: min 0.009, its middle free-paved).

**Look at:** the mesh (Prepare ▸ Mesh): rows running straight along the arm, the O-grids round
the bores, eight layers through the arm, and the von Mises contours: the bending stress grows
along the arm toward the root; the peak, 117 MPa, sits at a square-bore corner where the rigid
clamp meets the free face, a support singularity.

**Check:** σxx on the arm's sides at x = 29 mm is −62.44 and +62.44 MPa, equal and opposite as
the mirror-image halves are, against the beam formula
M/Z = 1500 · 141 / (16 · 35.6² / 6) = 62.6 MPa (0.25 %); a 41 815-element tet10 solution of the
full part gives 62.3 MPa. The pedal bore moves 0.75 mm down (tet10: 0.70–0.75 mm; Euler–Bernoulli
with a rigid clamp: 0.58 mm, without shear or the boss's compliance), and the spindle carries the
whole 1500 N and 255 N·m.

## Cantilever (`cantilever`)

![Von Mises stress on the clamped steel beam under its end load, highest at the clamped root](media/examples/cantilever.png)

A steel beam of 160 hex8 elements, clamped at one end with a load on the other: a static step
and a frequency step. **Look at:** the bending stress, highest at the clamped root, and the
Modes step's mode shapes.

## Bar impact (`impact`)

![Von Mises stress in the steel beam at the end of its impact on the rigid wall](media/examples/impact.png)

The same beam flying onto a rigid wall: an explicit dynamic step of ten frames. **Look at:** the
stress wave running back along the bar (Results ▸ Frame) and the energy history in the Results
panel.

## Hex and tet parts (`parts`)

![A hex8 block beside a tet4 cube](media/examples/parts.png)

A hex8 block beside a tet4 cube: meshes to look at (element types, nodes, picking), nothing to
solve.

## I-beam cantilever (`ibeam`)

![Von Mises stress on the swept hex8 I-beam, clamped at x = 0 with 1 kN on its free end](media/examples/ibeam.png)

The 200 mm steel I-beam of the README video (`tests/corpus/cantilever-i-beam`), opened as
SolidCraft geometry with persistent face names, swept all-hex (hex8 of 4 mm, incompatible modes),
clamped at `F2:start` and loaded with 1000 N down through a kinematic coupling to `F2:top`.
**Look at:** the bending stress in the flanges, highest at the clamp; the tip moves 0.18 mm down
(beam theory with shear: 0.184 mm).

## Heated flow past a bracket (`heated-flow`)

Air at 20 °C enters at 1 mm/s and passes the original upright bracket, whose faces are
held at 80 °C. The 130 × 80 × 70 mm domain has 5 mm cells (5,824 cells), a pressure outlet,
slip sides and top, and a no-slip floor. There are 180 steps and six saved output intervals,
plus the initial frame. Open `file.sample {"name":"heated-flow"}` and solve on CPU or GPU.
The slow prescribed flow keeps this coarse lattice away from its small-relaxation-time limit.

**Look at:** Results → Temperature, or `results.slice {"field":"temperature"}`, for the
warming air. `results.surface {"field":"wall_heat_flux"}` shows heat received by each
bracket face per area: the hot wall heating the air is negative. In this model's mm–tonne–second
units the wall flux is mW/mm². The initial frame has no recorded wall heat and stays uncoloured.

**Check:** the fluid warms above its inlet temperature while staying below the held wall
temperature (with a small numerical tolerance); recorded hot walls give heat to the fluid.
This coarse, short transient is a quick display demonstration, not a converged convection
benchmark. It has neither buoyancy nor solid heat conduction. The sample is isolated in
`crates/engine/src/sample/heated_flow.rs` and its gallery registration.
