# OpenFOAM cases

`deck.openfoam` writes a fluid step as an OpenFOAM case directory, so the same model can be run
with a finite-volume solver beside the lattice: to check the lattice's answer, or for flow it cannot
run (high Mach numbers, body-fitted wall layers). FieldCraft only writes the case; OpenFOAM runs as
its own program from its own install.

```json
{"command": "deck.openfoam", "params": {"path": "nozzle_case", "processes": 8,
  "refine_boxes": [{"min": [40, 0, 0], "max": [190, 8, 8]}]}}
```

Then, in the case directory with OpenFOAM's environment: `sh Allrun` (`NP=4 sh Allrun` changes the
ranks), and open `case.foam` in ParaView.

## What is written

| Model | Case |
|---|---|
| fluid domain box | `system/blockMeshDict`: one block, cells of `spacing × cell_factor`; patches `x_min` … `z_max` |
| obstacle bodies | `constant/triSurface/<body>.stl`: ASCII, one solid per named face (tessellated as `solve.fluid` does), metres |
| | `system/snappyHexMeshDict`: obstacle surfaces and `refine_boxes` refined back to the spacing; feature edges from `surfaceFeatureExtract` |
| `velocity` face | `U` fixed, `p` zero gradient |
| `pressure` face | `p` fixed at pressure / density (kinematic), `U` inflow/outflow (`pressureInletOutletVelocity`) |
| `slip` face | `symmetryPlane` (box faces are planar) |
| `wall` face, obstacle face | no-slip; a `moving_wall` face or named obstacle surface: `U` fixed at its velocity |
| Newtonian material | `constant/transportProperties`: ν = μ / ρ |
| `turbulence: smagorinsky` | `constant/turbulenceProperties`: LES Smagorinsky with `Ck = (Cs² √Ce)^(2/3)`, `Ce` 1.048, cube-root-volume filter width; `van_driest` → Van Driest damping with that A⁺ |
| step length and frames | `system/controlDict`: `pimpleFoam`, `endTime`, `writeInterval` = length / frames, time step adapted to `max_courant` (0.8), a running mean of U and p from a third of the run |
| `processes` | `system/decomposeParDict` (scotch) and the ranks `Allrun` uses |

Each obstacle face becomes the patch `<body>_<face>` (characters OpenFOAM does not take become `_`,
and `+`/`-` become `p`/`m`: `Nozzle_F1_px`). The reply maps every patch back to its model face.
`snappyHexMesh` keeps the side of the obstacles that holds `location`; when it is left out, a point
just inside a pressure face (else a velocity face, else the box's centre) that lies in no obstacle
is used.

Everything is SI, as the model stores it and OpenFOAM reads it. The lengths of `refine_boxes` and
`location` are given in the input units.

## Not written yet

The command fails and names what it cannot write yet:
- heat in the fluid
- porous zones, fans and thin walls
- free surfaces
- body accelerations
- periodic faces
- parabolic inlet profiles
- Vreman's model
- non-Newtonian viscosity

What it writes differently from the lattice step is listed in the reply's `notes`:
- the lattice's time step is not used
- the steady stop rule is not mapped
- the wall models are not mapped (walls are resolved no-slip)
- a convective outflow becomes an inflow/outflow pressure face

## Licence

The case files are written from OpenFOAM's public user guide (case layout, dictionary syntax,
boundary conditions); no OpenFOAM source is read, copied or linked. OpenFOAM is run as a separate
program by the user.
