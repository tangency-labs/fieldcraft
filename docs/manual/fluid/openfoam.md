# OpenFOAM cases

`deck.openfoam` writes a fluid step as an OpenFOAM case directory, an export for running the case
in OpenFOAM. FieldCraft writes the case and nothing else: OpenFOAM runs as its own program from its
own installation. The case is written for ESI OpenFOAM v2406 or later (openfoam.com), and every file
says so in its header. FieldCraft's own fluid results are checked against theory, experiments and
published benchmarks (see [verification](../../verification.md)), never against OpenFOAM's.

```json
{"command": "deck.openfoam", "params": {"path": "channel_case", "processes": 8,
  "refine_boxes": [{"min": [40, 0, 0], "max": [190, 8, 8], "spacing": 0.5}]}}
```

Then, in the case directory with OpenFOAM's environment: `sh Allrun`, then open `case.foam` in
ParaView. The script runs on the `processes` ranks the case was written for (its
`decomposeParDict`); write the case again to change them.

Every setting the case depends on is written into it, the mesh-quality limits included: nothing is
taken from the installation's own files (no `#include`), so the same case meshes and runs the same
on every installation that reads it. Each fixed entry names, in the source beside it, the section of
OpenFOAM's public user guide, tutorial guide or documentation site it follows, or how its value is
derived.

## What is written

| Model | Case |
|---|---|
| fluid domain box | `system/blockMeshDict`: one block, cells of `spacing × cell_factor` (1 by default); patches `x_min` … `z_max` |
| obstacle bodies | `constant/triSurface/<body>.stl`: ASCII, one solid per named face (tessellated as `solve.fluid` does), in metres |
| | `system/snappyHexMeshDict`: obstacle surfaces refined back to the spacing, `refine_boxes` to their own; feature edges from `surfaceFeatureExtract` |
| `velocity` face | `U` fixed, `p` zero gradient |
| `pressure` face | `p` fixed at pressure / density (OpenFOAM's incompressible pressure is kinematic), `U` inflow/outflow (`pressureInletOutletVelocity`) |
| `slip` face | `U` slip, on a plain patch |
| `symmetry` face | `symmetryPlane` (box faces are planar) |
| `wall` face, obstacle face | no-slip; a `moving_wall` face or a named obstacle surface condition: `U` fixed at its velocity |
| Newtonian material | `constant/transportProperties`: ν = μ / ρ |
| step turbulence | `constant/turbulenceProperties`: laminar without one; Smagorinsky's model as OpenFOAM's LES Smagorinsky, its `Ck` derived from the step's constant |
| step length and frames | `system/controlDict`: `pimpleFoam`, `endTime`, a result every length / frames, the time step adapted to `max_courant` (1 by default, the most it may be) |
| step averaging | a running mean of U and p from the step's own averaging start (none when the step does not average) |
| `processes` | `system/decomposeParDict` (scotch) and the ranks `Allrun` uses |

Each obstacle face becomes the patch `<body>_<face>`: characters OpenFOAM does not take become `_`,
and `+`/`-` become `p`/`m` (`Plate_F1_px`). The reply maps every patch back to its model face.

`snappyHexMesh` keeps the side of the obstacles that holds `location`, which must lie on no cell
face. The point written is always the centre of a cell of the finest refinement level, and the
lattice's exact inside test decides that it is in the fluid:
- when `location` is given, the centre of the finest cell holding it (the reply shows the point
  written); a point inside an obstacle is refused;
- when it is left out, the first fluid lattice cell whose neighbours are all fluid, flooding in from
  the inlets and outlets (or from the box's faces when none is open), so the point is in the fluid
  the open faces reach, never in a closed pocket.

Lengths are given in the input units; the case is SI, as the model stores it.

## Reading the results back

`results.openfoam {"path": "channel_case"}` reads the results of a case `deck.openfoam` wrote and
OpenFOAM ran as the fluid results of its step, for viewing them with FieldCraft's slices, probes,
lines and measurements. OpenFOAM's results are shown as OpenFOAM computed them; they are never
used to judge FieldCraft's own.

- Each written time is a frame, sampled at the centres of a lattice over the domain box (the
  domain's spacing, or `spacing`). Each fluid lattice cell takes the values of the OpenFOAM cell
  whose centre is nearest; the obstacles' cells are solid.
- Velocity, and the gauge pressure: OpenFOAM's incompressible pressure is kinematic, so it is
  multiplied by the density the case was written with (recorded in `fieldcraft.json` beside the case).
- The cell centres come from `0/C`, which `Allrun` writes. The fields must be ASCII, as the case is
  written.
- A case that was not reconstructed (only `processor0` … `processorN`) is read from its processor
  directories, each rank's centres and fields joined in rank order.
- Fewer frames: `from` and `to` (the written times in this range, input units, inclusive) and
  `every` (every n-th of them, always with the last).

What is refused, rather than read as something else:
- a directory without `fieldcraft.json`, or one whose record lacks a value;
- a case none of whose cell centres lies in the step's domain (another model's case);
- a fluid lattice cell with no OpenFOAM cell near it, which would otherwise show as zero flow. "Near"
  is bounded by the case's largest cell: the background cell's diagonal, plus the distance
  snappyHexMesh may move points onto a surface;
- a field whose count differs from its mesh's cells, binary files, and counts beyond the memory the
  fluid results may take.

## Compressible cases

`"compressible": {"ambient": 0.101325, "temperature": 20, "gas_constant": 287.05, "far_field": 50}`
writes a `rhoCentralFoam` case instead: a laminar ideal gas, for flows the lattice cannot run, such
as a choked nozzle. Every value is required; no gas is assumed.

| Input | What it is |
|---|---|
| `ambient` | the absolute pressure the model's gauge pressures count from |
| `temperature` | the initial temperature, and the total temperature of the reservoirs |
| `gas_constant` | the specific gas constant R/M of the gas (air: about 287 J/(kg·K)) |
| `far_field` | the far-field relaxation length of the non-reflecting outlets (OpenFOAM's `lInf`): the distance over which an outlet's pressure relaxes towards its value |
| `initial_pressure` | boxes of the fluid that start at their own gauge pressure (`setFields`), such as a charged reservoir |
| the material | the specific heat, viscosity and conductivity (isotropic); its density is not used |

| Model | Compressible case |
|---|---|
| `pressure` face above the ambient | a reservoir: total pressure and total temperature |
| other `pressure` faces | a non-reflecting outlet at that pressure |
| `velocity` face | `U` fixed, `T` at `temperature` |
| walls, obstacle faces | no-slip and adiabatic |
| time step | adapts to `max_courant`: 0.2 unless given, at most ½ (see below) |

A turbulence model is refused in a compressible case: it is laminar.

### The time-step limit

The Courant number of a cell is C = Σ over the axes i of (|u_i| + c) Δt / Δx_i, the sum of the
wave speeds crossing each cell face's axis, as OpenFOAM's `CourantNo` sums it. The case adapts
the step to `max_courant`.

The case writes the central-upwind flux of Kurganov, Noelle and Petrova (2001, *SIAM J. Sci.
Comput.* 23), which the solver applies in semi-discrete form.

The default is 0.2: the CFL number at which Greenshields, Weller, Gasparini and Reese (2010, *Int. J.
Numer. Meth. Fluids* 63, §4.1 and §4.2) ran their central-scheme cases, a shock tube and a forward
step, and their figures 2 and 3 use it too.

The ceiling ½ is the stability limit of the central schemes that the same paper states in §4.1
("below the stability limit of 0.5 of the central schemes"). The paper states it and does not derive
it, and this document does not prove it. `max_courant` may not pass it.

### Specific heats and γ

γ = c_p / (c_p − R) is derived from the material's c_p and the gas constant, not given. An ideal gas
has c_v ≥ 3R/2 (equipartition over the three translational quadratic terms; Kittel and Kroemer,
*Thermal Physics*, 2nd ed., ch. 3), so γ ≤ 5/3, with equality for a monatomic gas. A γ above 5/3
means the data are not an ideal gas's. The case is still written, with a note on the case naming γ:
data rounded to four figures can sit just above 5/3 for a monatomic gas (argon, with c_p = 520.3 and
R = 208.13, is 0.0032 % above it).

Reading a compressible case back gives each cell's gauge pressure (the absolute pressure less
`ambient`), its temperature, and its density ρ = p/(R T). A cell whose temperature or pressure is
not positive is refused. The ports' mass flow, mixed temperature and heat then weigh each cell by
its own density.

Results saved by FieldCraft 0.8 and earlier (the `MCRSLT05` results encoding and older) open with their model and
without their results; solve the step again.

## Not written yet

The command fails and names what it cannot write yet:
- heat in the fluid
- porous zones, fans and thin walls
- free surfaces
- body accelerations
- periodic faces
- parabolic inlet profiles
- Vreman's and the Spalart–Allmaras models; any turbulence model in a compressible case
- non-Newtonian viscosity, and a density or viscosity that follows the temperature
- a material without a density or a viscosity (no value is assumed)
- a case directory that already holds files: choose a new or empty one (nothing is overwritten or deleted)

What it writes differently from the lattice step is listed in the reply's `notes`:
- the lattice's time step is not used
- the steady stop rule is not mapped
- the wall models are not mapped (walls are resolved no-slip)
- a convective or characteristic outflow becomes an inflow/outflow pressure face
- Van Driest's damping is written as OpenFOAM's law, min(κ y D, C Δ) with Cdelta the step's constant, not the lattice's C Δx D

## Licence and clean room

The case files are written from OpenFOAM's public documentation (case layout, dictionary syntax,
boundary conditions); no OpenFOAM source or tutorial case is read, copied or linked. OpenFOAM is
run by the user as a separate program.
