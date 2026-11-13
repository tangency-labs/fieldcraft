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

## Compressible cases

`"compressible": true` (air at 1 atm and 20 °C) or `{ambient, temperature, molar_mass,
specific_heat, prandtl}` writes a `rhoCentralFoam` case: a laminar ideal gas (Kurganov's central
scheme), for flows the lattice cannot run, such as a choked nozzle and its underexpanded jet.

| Model | Compressible case |
|---|---|
| `pressure` face above 0 | inlet: total pressure `ambient` + gauge, total temperature `temperature` |
| `pressure` face at or below 0 | non-reflecting outlet (`waveTransmissive`) at `ambient` + gauge |
| `velocity` face | `U` fixed, `T` at `temperature` |
| walls, obstacle faces | no-slip, adiabatic |
| material | viscosity; `specific_heat` from its thermal data unless given (else air's); the density is not used |
| time step | adapts to `max_courant`, 0.2 unless given |

A high-pressure inlet against still ambient air starts as a shock tube at the inlet face, and a
strong pressure jump that crosses snapped wall cells can drive the temperature negative in the
first steps. `initial_pressure: [{min, max, pressure}]` (gauge, input units) starts boxes of the
fluid at their own pressure (`setFields`). Charge the reservoir and nozzle up to the exit plane, so
the jump starts where the wall is a grid-aligned face:

```json
{"command": "deck.openfoam", "params": {"path": "blast_case", "processes": 8,
  "compressible": {"ambient": 0.1, "temperature": 20},
  "initial_pressure": [{"min": [-1, -1, -1], "max": [70, 100, 100], "pressure": 29.9}],
  "refine_boxes": [{"min": [65, 0, 0], "max": [112, 12, 12], "spacing": 0.4},
                   {"min": [55, 0, 0], "max": [100, 4, 4], "spacing": 0.2}]}}
```

## Reading the results back

`results.openfoam {"path": "blast_case"}` reads every written time of a reconstructed case as
a frame of the step's fluid results, sampled at the centres of a lattice over the domain box (the
domain's spacing, or `spacing`) from the nearest OpenFOAM cell; the obstacles' cells are solid. The
slices, probes, lines and measurements of fluid results then show OpenFOAM's flow:
- velocity
- gauge pressure: above the ambient (compressible), or the kinematic pressure times the density
- temperature (compressible)
- density (compressible): each cell's, from the ideal gas ρ = p / (R T); probes show it

Compressible results weigh each port cell by its own density: the inlet/outlet mass flows,
mixed temperatures and sensible heat of `results.fluid.measure` and the results panel are
`Σ ρ |u · n| Δx²`, not the material's reference density (which the lattice's nearly
incompressible flows use).

The reply gives the peak speed and, compressible, the peak Mach number. The cell centres come
from `0/C`, which `Allrun` writes; fields must be ASCII, as the case is written.

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
