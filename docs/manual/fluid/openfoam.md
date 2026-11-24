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

## Not written yet

The command fails and names what it cannot write yet:
- heat in the fluid
- porous zones, fans and thin walls
- free surfaces
- body accelerations
- periodic faces
- parabolic inlet profiles
- Vreman's and the Spalart–Allmaras models
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
