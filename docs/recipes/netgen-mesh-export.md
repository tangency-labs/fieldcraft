# Export a mesh to Netgen or NGSolve

Choose **Mesh → Export → Export Netgen Mesh**. The file field starts with the model name
and a `.vol` extension. Enter or paste a path, or use **Browse…** to choose a destination.
FieldCraft adds `.vol` if you omit the extension. Confirm the form to write the mesh.
The exporter writes counterclockwise linear Tri3 and Quad4 meshes in the XY plane at z = 0 and linear Tet4 and Hex8
meshes in 3D.
Coordinates are in meters, regardless of the model's display units. The file carries the
existing connectivity and region names. Define material properties, boundary conditions, and
loads in your analysis script. See the [VOL format description](../formats/vol.md).

- **Region sets** name the region containing each element without assigning a material. When supplied, the sets
  must cover every element. When omitted, the exporter uses section material names, CAD body
  names for unassigned elements, and `default` for untagged elements.
- **Boundary sets** name 2D edges or 3D faces, on the outside or at an interior interface.
  Select a shared edge or face from either side or both sides of the same set. The file stores it
  once with both adjacent domains. In 2D, a name used across multiple domain pairs keeps that
  name on a distinct segment index for each pair. Separate meshes retain their separate nodes.
  Other exterior boundaries and interfaces between regions use their qualified CAD names when
  available or fall back to `default` and `interface`, respectively.
  The exporter joins the CAD body and entity names with a period and replaces whitespace, control
  characters, and `#` with underscores. If different boundaries produce the same VOL name,
  assign distinct boundary sets.

Open the result as a mesh in Netgen or load it from your existing NGSolve script:

```python
from ngsolve import Mesh

mesh = Mesh("team13.vol")
print(mesh.GetMaterials())
print(mesh.GetBoundaries())
```

## TEAM 13 nonlinear magnetostatics

Johannes Maierhofer built the [editable SolidCraft model](team13.solidcraft) for the
FieldCraft workshop from the published TEAM 13 drawing. Open it in SolidCraft to inspect the
center plate, two offset steel channels, eight winding segments, and surrounding air. The design
uses millimeters. The steel is 3.2 mm thick, the plate-to-channel clearance is 0.50 mm, and the
winding is 100 mm high with inner and outer corner radii of 25 and 50 mm. The workshop adds a
500 mm cube of air. [COMPUMAG problem 13, Fig. 1 and Table 1](https://www.compumag.org/wp/wp-content/uploads/2018/06/problem13.pdf)
· [Nakata, Takahashi, and Fujiwara, 1995](https://doi.org/10.1108/eb010140)

The [Rust generator](../../apps/fieldcraft-cli/examples/netgen_team13.rs) builds the
[FieldCraft command recipe](netgen-team13.json) from the planes of the steel and winding.
It maps square rings to the circular winding corners and returns to a Cartesian grid at the
steel flanges. The partition planes and 3,600 single-cell blocks are original workshop
meshing choices. FieldCraft divides each block into six linear tetrahedra, with shared nodes
across material interfaces. The recipe runs independently of the CAD file.

Run the checked-in recipe to export 21,600 tetrahedra to `team13.vol`.

```sh
fieldcraft-cli run docs/recipes/netgen-team13.json
```

To regenerate the recipe from the repository root, run the generator.

```sh
cargo run --quiet -p fieldcraft-cli --example netgen_team13 > docs/recipes/netgen-team13.json
```

The export has ten named regions. The `iron` set combines the plate and both channels,
and the eight winding segments retain their own names. The `outer` boundary encloses the air,
`gamma_core` and `gamma_coil` border steel and winding against air, and `winding_seams` joins
adjacent winding segments. For an already-open mesh, select these sets in **Export Netgen Mesh**.

## Planar MEMS electrostatic actuator

The [MEMS pull-in mesh recipe](netgen-mems-pull-in.json) triangulates the air gap and silicon
beam as separate XY regions. It uses the element sets `air` and `beam` as region names and boundary sets
`ground`, `outer`, `beam_electrode`, `tip`, `beam_top`, and `clamp`. Both sides of
`beam_electrode` retain the same label and their own mesh nodes.

Run the recipe to build and export `pull-in.vol`.

```sh
fieldcraft-cli run docs/recipes/netgen-mems-pull-in.json
```
