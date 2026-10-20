# Examples

Every example is a built-in sample: open it from **File ▸ Samples** (or the start screen), with
`fieldcraft --sample NAME`, or with the `file.sample` command (`{"name": "drop-test"}`), then
press **Solve** (`solve.run`). Each opens fully defined — geometry, mesh, materials, boundary
conditions, loads and steps — in mm, t, s, N and MPa, and opens and solves within half a minute (the electronics cooling studies at the end take minutes). The geometry is our
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

# Electronics cooling studies

Three larger thermal lattice-Boltzmann studies, samples like the rest but slower to solve: the heat
sink takes about four minutes on a CPU, the enclosure about twenty, and the AIO radiator (5 M cells)
needs a GPU (`solve.fluid` with `"backend": "gpu"`). Their setup is
`crates/engine/src/sample/cooling.rs`, their geometry SolidCraft recipes (the `cooling-*` parts of
[`tests/corpus`](../tests/corpus) and `crates/engine/samples/cooling-aio*.json`); trunk tests solve
the heat sink and the enclosure briefly, heavy tests the AIO's start and the heat sink's and the
enclosure's numbers below.
`tools/cooling/make.py` (see its [README](../tools/cooling/README.md)) opens each sample headless,
solves it, prints the numbers and renders the pictures. They show what the fluid solver can do for
electronics cooling today, and where it stands in for what it cannot yet do:

- **No conjugate heat transfer into solids yet**
  ([#157](https://github.com/tangency-labs/fieldcraft/issues/157)): a high-conductivity heat sink
  or radiator is an **isothermal surface**, and a component a **held-temperature face**. Real
  junction temperatures need conduction through the solid, and a component given its power in W
  needs [#609](https://github.com/tangency-labs/fieldcraft/pull/609).
- **Fans are a uniform inlet velocity** here. The solver now has fan curves and porous vents
  ([#596](https://github.com/tangency-labs/fieldcraft/pull/596)); the studies keep a fixed face
  velocity because the two-relaxation-time collision is resolution-bound near its stability limit
  ([#557](https://github.com/tangency-labs/fieldcraft/issues/557)), which caps them at a gentle
  0.1 m/s, below where a fan curve's operating point would sit. A fan zone in the AIO's shroud is
  the natural next step.
- **No surface-to-surface radiation**: the sealed enclosure is **convection only** (radiation
  carries a real share of the heat in a sealed box).
- The colour legend under each picture is the temperature range it is clamped to.

Background for the planned electronics-cooling release is in
[discussion #584](https://github.com/tangency-labs/fieldcraft/discussions/584).

## Heat sink in forced air (`cooling-heat-sink`)

![Air at 25 °C entering a six-fin heat sink held at 50 °C, warming between the fins, the warm layer lifting off the sink into a rising wake, streamlines coloured by air temperature](media/examples/cooling-heat-sink.png)

A 40 × 40 mm aluminium heat sink with six 18 mm fins, held isothermal at 50 °C (a
high-conductivity sink with fin efficiency near 1), in a closed 120 × 44 × 30 mm duct. Air at
25 °C enters at 100 mm/s (0.1 m/s, a gentle fan), the outlet at zero pressure. **Buoyancy is on**
(Boussinesq, β = 3.4 × 10⁻³ K⁻¹, g along −z): at 0.1 m/s this is mixed convection (Re ≈ 270 on the
fin length, Ri = Gr/Re² ≈ 3), and the 25 K rise keeps βΔT at 0.085, well inside Boussinesq. D3Q19
TRT with the D3Q7 thermal lattice on 158,400 cells of 1 mm, 12,000 steps of 0.25 ms to 3 s. The
wake behind the sink keeps shedding, so the steady stop rule never fires; the heat holds within
about 1 % from 2 s on (1.92–1.95 W out to 20 s).

**Numbers:** the sink sheds **1.93 W**, a sink-to-inlet **resistance of 13.0 K/W**. The duct's
air could carry 3.98 W if it all left at the sink's temperature; it carries **48 %** of that
(**ε = 0.48**), a bulk rise of **12 K** at the outlet. So at this airflow the sink is
airflow-limited: per area it does about as well as the bare chip below (h ≈ 7 W/m²K over its
110 cm² against the chip's 7.9), and its lower resistance comes from its area. Buoyancy lifts the
warm layer off the fins into a rising wake and adds about 4 % (1.85 W without it at 3 s).

**Look at:** the temperature slice between the middle fins and the streamlines coloured by
temperature (`results.slice` / `results.streamlines` with `field: "temperature"`), the hot layer
rising behind the sink. The heat each face gives the air is in the `solve.fluid` report (`heat`,
negative into the fluid).

**Caveat:** the air's minimum is **20.4 °C against the 25 °C inlet** (4.6 K under), at the sharp
fin corners: the D3Q7 anti-bounce-back wall at a convex edge, not a physical undershoot. The
maximum overshoots the sink by 0.5 K at the same corners.

**Regenerate:** `python3 tools/cooling/make.py heat-sink`.

## Sealed enclosure with a PCB (`cooling-enclosure`)

![Side section through a chip on a vertical PCB in a sealed box: the warm plume rising off the chip, spilling over the top of the board and recirculating down the walls](media/examples/cooling-enclosure.png)

A small package held at 40 °C on a vertical PCB inside a sealed 60 × 40 × 60 mm air box whose
walls are held at 25 °C, no fan. Buoyancy (Boussinesq, β = 3.4 × 10⁻³ K⁻¹, g along −z) drives a
plume up off the chip, over the top of the board and down the cool walls. D3Q19 TRT with the D3Q7
thermal lattice and buoyancy on 144,000 cells of 1 mm at a fixed 0.2 ms step, run to **steady
state** with the steady stop rule ([#595](https://github.com/tangency-labs/fieldcraft/pull/595):
flow, temperature and the chip's heat each changing by under 0.1 % per 2 s window): it stops at
18 s, after 90,000 steps, when the walls take back the chip's heat to within 0.01 %.

**Numbers:** the chip sheds **62.7 mW** by natural convection: over its 528 mm² of wetted faces at
15 K that is **h ≈ 7.9 W/m²K**, against 9.5 from the Churchill–Chu correlation for a vertical
plate of the chip's height, and a chip-to-wall resistance of **239 K/W**. The air stays between
the walls' and the chip's temperatures (24.9–39.5 °C). The plume's 72 mm/s peak matches its scale
√(gβΔT·L) ≈ 77 mm/s. Solved on 2 mm cells instead, the steady heat is 1.3 % lower.

**Look at:** the side section through the chip (`results.slice` with `axis: "x"` at 35): the
plume rising off the chip, spilling over the board's top edge, and the recirculation the
streamlines draw down the walls.

**Caveats:** the buoyant **start-up** needs a small fixed time step: the solver's automatic step
over-accelerates the first instant and goes unstable (lattice speed 0.505 at step 38). This is
**convection only**: a sealed enclosure also loses a real share of its heat by **surface-to-surface
radiation**, which the solver does not model yet.

**Regenerate:** `python3 tools/cooling/make.py enclosure`.

## AIO cooler radiator (`cooling-aio`)

![A 240-class all-in-one liquid cooler: a cold plate and a two-tier pump block, two coolant tubes running side by side to the radiator, whose fine air channels sit in a recessed pocket behind a rounded raised border](media/examples/cooling-aio.png)

A 240-class all-in-one liquid cooler, modelled in SolidCraft for this study
([`tools/cooling/aio_recipe.py`](../tools/cooling/aio_recipe.py) writes the `cooling-aio` and
`cooling-aio-half` recipes): a cold plate and a two-tier pump block, two Ø6.5 mm coolant tubes
running as a parallel pair from barbs on the pump to ports on the radiator, and a
230 × 118 × 30 mm radiator whose fin pocket sits recessed behind a filleted, raised border, with
2.4 mm air channels cut through its whole depth.

**Only the radiator's air side is simulated.** The water-cooling half of the cooler (coolant
pumped through the cold plate, the tubes and the radiator's core, and the heat conducted from the
CPU into it) is still to do: the cold plate, the pump and the tubes are drawn for context, and the
radiator stands in for the coolant by being held at an assumed coolant temperature.

**The study:** a fan shroud seals to the radiator's frame, so air at 25 °C entering at 0.1 m/s
across the face has to pass through the channels; the radiator is held at a 45 °C coolant
temperature. The sample is the half of the cooler with the ports (`cooling-aio-half`), with a
symmetry plane on the radiator's middle, which is also what lets the picture look into the
channels. D3Q19 TRT with the D3Q7 thermal lattice on 5.07 M cells of 0.6 mm (10.1 M for the whole
radiator), 20,513 steps to 1.6 s, solved on the GPU (NVIDIA GB10, single precision) in about four
minutes.

![Cutaway through the radiator's middle plane: cool air arriving from the left, each channel heating its air between the fins, warm jets leaving into the exhaust behind, slow recirculation behind the frame's top and bottom](media/examples/cooling-aio-airflow.png)

**Numbers:** 163 L/min of air (5.8 cfm) goes through the radiator, and the radiator gives it
**61 W**: an effectiveness of **0.93**, of the 66 W the air could carry if it left at the
coolant temperature. At this airflow the radiator is nearly saturated: the air leaves almost as
warm as the coolant.

**Checks:** refining the cells from 0.6 to 0.5 mm moves the heat by 0.4 %. That shows the heat
is converged, but says little about the channels' heat transfer coefficient: at NTU ≈ 3 a 10 %
change in h moves the heat by only about 1.6 %. The time step matters more, through the weakly
compressible error of the fast channel flow: the heat falls from 69 to 64 to 61 W as the inlet
lattice speed goes from 0.03 to 0.02 to 0.013 (a density swing of 9.3, 3.9 and 1.6 %), linearly
in its square, and extrapolates to about 59 W (effectiveness 0.90) at zero time step. The study
uses 0.013, within about 3 % of that. The air's minimum is 20.8 °C against the 25 °C inlet, at
the channels' sharp corners (the same wall artefact as the heat sink's).

**Look at:** the cutaway temperature slice (`results.slice` with `axis: "x"` on the cut): cool air
in front, each channel heating its air between the fins, warm jets into the exhaust, and slow
recirculation behind the frame's top and bottom.
[The clip](media/examples/cooling-aio-airflow.mp4), recorded in the desktop app on a longer domain
(130 mm of exhaust behind the radiator, 12.6 M cells to 2.4 s), shows the start-up: the warm front
leaving the channels with a vortex curling at either end of the frame, then the settled exhaust.

**Caveats:** 0.1 m/s is the ceiling of today's collision operator: the channel flow reaches about
470 mm/s, a cell Reynolds number of about 19 at 0.6 mm, where the heat-sink study found the
two-relaxation-time collision going unstable
([#557](https://github.com/tangency-labs/fieldcraft/issues/557)). A pair of 120 mm fans pushes 1–2
m/s through a real radiator, where it no longer saturates and the fin design decides the heat,
so this study cannot rank radiator designs yet. The water loop is to do: the radiator is
isothermal, with no coolant flow and no conduction inside it
([#157](https://github.com/tangency-labs/fieldcraft/issues/157) and a coupled coolant loop), so
the coolant temperature is an input here rather than a result. The heavy tests do not hold the
AIO's numbers: its run is 10¹¹ cell updates, hours on a CPU; a heavy test solves a coarse start.

**Regenerate:** `python3 tools/cooling/make.py aio --backend gpu` (the CPU solver also runs it,
far more slowly at 5 M cells).
