# Composite shells and material directions

Composite and oriented shell sections currently use the flat `tri3`/`quad4` kernels.
The curved `quad8` MITC8 path accepts only homogeneous isotropic sections and refuses
layups and material orientations; its local frame is not a ply-coordinate implementation.

Define a material, assign its directions and an ordered layup to a shell set, solve a linear
static, frequency or flat-plate buckling step, then inspect each ply through `solve.probe`.
The CLI, scripts and MCP use these same commands. A first-ply failure index does not change
the stiffness or simulate progressive damage.

## Define and inspect

This example assumes an existing shell element set named `panel`. Explicit units keep the
input independent of the active display system:

```json
[
  {"command":"model.material.create","params":{
    "name":"Ply",
    "elastic":{"type":"engineering_constants","young":["120 GPa","9 GPa","9 GPa"],
      "poisson":[0.25,0.25,0.3],"shear":["5 GPa","4 GPa","3 GPa"]},
    "density":"1600 kg/m^3",
    "ply_strength":{"tension":["1000 MPa","40 MPa"],"compression":["800 MPa","120 MPa"],
      "shear12":"50 MPa","shear23":"60 MPa","tsai_wu_interaction":-0.5}}},
  {"command":"model.section.create","params":{
    "name":"Skin","set":"panel",
    "kind":{"type":"composite_shell","layup":{"plies":[
      {"material":"Ply","thickness":"0.125 mm","angle":"0 deg","points":3},
      {"material":"Ply","thickness":"0.125 mm","angle":"90 deg","points":3},
      {"material":"Ply","thickness":"0.125 mm","angle":"90 deg","points":3},
      {"material":"Ply","thickness":"0.125 mm","angle":"0 deg","points":3}],"offset":0}},
    "orientation":{"frame":{"type":"rectangular","axis1":[1,0,0],"axis2":[0,1,0]},
      "shell_reference_axis":1}}},
  {"command":"model.material.get","params":{"name":"Ply"}},
  {"command":"model.section.get","params":{"name":"Skin"}},
  {"command":"model.section.properties","params":{"name":"Skin"}}
]
```

`model.material.edit` and `model.section.edit` accept the corresponding `get` reply.
A composite section omits the homogeneous `material` field: each ply names its own material.
Materials referenced by a ply cannot be deleted without resolving the reference.

| Parameter | Meaning |
| --- | --- |
| `young`, `poisson`, `shear` | E1/E2/E3, ν12/ν13/ν23, G12/G13/G23 in material axes |
| `plies` | 1–256 plies, from the negative to the positive element normal |
| `thickness` | Positive thickness of that ply; their sum is total thickness h |
| `angle` | Right-hand rotation about the positive normal after projecting the reference axis; radians unless display angle units are changed or specified explicitly |
| `points` | 1–15 extra Gauss samples per ply; bottom/middle/top are always available. Elastic thickness moments are exact even with one sample |
| `offset` | Nodal reference surface relative to the laminate centre, divided by h, between −0.5 and +0.5; ply coordinates have centre z = −offset·h |
| `orientation.frame` | Rectangular `axis1` and independent `axis2`, or cylindrical `origin` and axis direction `axis` |
| `shell_reference_axis` | Axis 1/2/3 to project onto the shell; cylindrical axes are radial/hoop/axial. A direction parallel to the normal is rejected |

Shell directions are evaluated at each facet's node-average centre. Omitted directions use
the element frame, which need not equal global X/Y. Cylindrical **solid** directions are
evaluated at each stiffness integration point. A solid section uses the same `orientation`
with `shell_reference_axis:1`, but no layup. See [solid directions](../../solvers/solid-material-axes.md).

`model.section.properties` returns A/B/D, transverse shear, thickness, centre z and mass
moments I0/I1/I2 in the section reference axes. The reply specifies each unit. Missing ply
density gives `mass_moments:null` and `missing_density`; frequency or gravity solves require
every ply's density. Use consistent mass to retain offset or unequal-density I1 coupling.

## Solve and read layers

Create the usual `static`, `frequency` or `buckle` step with its supports and loads, and use
`solve.static`, `solve.frequency` or `solve.buckle`. For a stored static result:

```json
{"command":"solve.probe","params":{"step":"Load","elements":[1]}}
```

Each shell element contains `axes`, resultants and legacy global bottom/middle/top stresses.
Its `plies[node_slot][ply_slot]` keeps each layer separate, with `angle`, `bottom`, `middle`,
`top`, and `points`. Each point has z, `stress`, `strain`, and optional `failure`. These are
element-nodal extrapolations; neighbouring elements and either side of a ply interface are
not averaged. Full vectors are 11,22,33,12,13,23 in material axes; shear strains are engineering
shears. Stress is tensile positive. z, stress and angle follow display units, while strain
and failure indices remain dimensionless. A frequency mode's amplitude is arbitrary and
does not define a physical failure index.

Stored ply fields are checked against the solver's `max_factor_entries × sizeof(f64)` memory
budget before assembly. A dense stack can exceed that budget even on a small mesh; the error
names the required storage instead of allocating an unbounded result.

The strength arrays contain positive tensile and compressive **magnitudes** along material
1 and 2. Missing strengths give `failure:null`. With strengths, the output gives maximum
stress, signed-product Tsai–Hill, raw Tsai–Wu, and the active 1980 Hashin sign branches and
their maximum. An index of 1 is the threshold. Tsai–Wu may be negative; quadratic indices
are not their square roots or inverse reserve factors. `tsai_wu_interaction=-0.5` is an
assumption, not a measured biaxial interaction. S23 is independent and is used in Hashin
matrix compression. That active Hashin branch is also retained if its polynomial is
negative; `maximum` compares the active fibre and matrix values. Equations and citations:
[laminate kernels](../../solvers/laminates.md).

## Verification and scope

Theory cases cover ABD, pure bending and thin/thick Navier statics and frequencies. The
NAFEMS R0093/1 case checks ply stress; the carbon/epoxy cantilever experiment checks measured
frequencies. [Oracle results](../../oracle.md) distinguish reference tolerances from frozen
regression fixtures. The FEM tests also check known-stress failure values, orientation
objectivity, offsets, unequal densities and thread-count determinism.

MITC3/MITC4 use first-order shear deformation with fixed 5/6 correction. Through-thickness
shear stress is a constitutive approximation, not an equilibrated interlaminar field. It
does not verify delamination. Buckling requires a flat, coplanar preload with no transverse
slope or director tilt; an unsymmetric layup's bending response can violate that condition.
Thermal laminate resultants, finite rotations, progressive damage, continuum-shell layups,
Puck and LaRC remain outside this capability.

The [deck mappings](../../composite-decks.md) state the supported directions, sample counts
and named refusals for Abaqus, Nastran and LS-DYNA. None of those mappings is a certificate
from running a commercial solver.
