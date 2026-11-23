# Netgen VOL export

FieldCraft writes a text mesh file with one mesh dimension. Coordinates are in meters. Point,
domain, and boundary numbers start at one. Domain zero denotes the space outside the mesh.
Every counted section starts with its record count on a separate line.

## Sections

The file starts with `mesh3d`, followed by `dimension` and its value (2 or 3) and `geomtype`
and its value (0). The final marker is `endmesh`. The table shows the columns in each section
FieldCraft emits. `n` is the number of point indices that follow a cell or face record.

| Section | Dimension | Record |
|---|---|---|
| `surfaceelements` | 3 | Surface descriptor, boundary number, inside domain, outside domain, `n`, point indices |
| `surfaceelements` | 2 | Domain, domain, 0, 0, `n`, point indices |
| `volumeelements` | 3 | Domain, `n`, point indices |
| `volumeelements` | 2 | Empty |
| `edgesegmentsgi2` | 2 | Segment index, 0, first point, second point, −1, −1, inside domain, outside domain, first edge index, first parameter, second edge index, second parameter |
| `points` | 2 or 3 | X, Y, Z in meters |
| `materials` | 2 or 3 | Domain, material name |
| `bcnames` | 3 | Boundary number, boundary name |
| `bcnames` | 2 | Segment index, boundary name |

Boyle's workshop labels the surface and volume record columns. Its EIT-specific 2D layout
differs from FieldCraft's `dimension 2` arrangement. The meshio VOL implementation shows the
standard dimension-dependent section handling, one-line `edgesegmentsgi2` fields, and name
table layout. FieldCraft writes −1 for both incident-triangle references in a segment and
fills the two edge indices with the segment index. The endpoint parameters are 0 and 1.
[Boyle, *zedhat Workshop*, slide 14](https://aboyle.ca/pubs/boyle2019d-workshop.pdf)
[meshio VOL reader and writer](https://github.com/nschloe/meshio/blob/main/src/meshio/netgen/_netgen.py)

## Connectivity and orientation

The table numbers each position in FieldCraft's reference element starting at one. The order
shown is the sequence written to the file.

| Element | FieldCraft positions in file order | Geometry |
|---|---|---|
| Tri3 | 1, 2, 3 | Counterclockwise viewed from +Z |
| Quad4 | 1, 2, 3, 4 | Counterclockwise viewed from +Z |
| Tet4 | 1, 3, 2, 4 | The normal of file positions 1–2–3 points away from position 4 |
| Hex8 | 1, 4, 3, 2, 5, 8, 7, 6 | Positions 1–2–3–4 form the outward-oriented base. Positions 5–6–7–8 are the corresponding vertices on the opposite face |

The official planar mesh guide gives counterclockwise triangle and quadrilateral connectivity.
For volume cells, the meshio VOL table gives the conversion between a positive reference cell
and VOL's outward-oriented first face. For a Tet4, the signed volume in this convention is
`det(p3 - p1, p2 - p1, p4 - p1) / 6`. For an affine Hex8, it is
`det(p4 - p1, p2 - p1, p5 - p1)`, where `p1` through `p8` refer to file positions.
FieldCraft requires planar cells to lie at Z = 0 with strictly counterclockwise corners.
The exporter checks the sign of tetrahedral volume and the orientation at each Hex8 corner
with exact predicates. A Hex8 also needs a positive Jacobian over its entire reference cube.
[Netgen manual mesh guide](https://docu.ngsolve.org/latest/netgen_tutorials/manual_mesh_generation.html)
[meshio VOL connectivity table](https://github.com/nschloe/meshio/blob/main/src/meshio/netgen/_netgen.py)

In 3D, the cross-product normal of a surface record points out of its inside domain and into
its outside domain. In 2D, a segment's inside domain is on its left when facing from the first
point toward the second. FieldCraft writes a single record for a shared edge or face and
preserves each side's domain number.
[Netgen working with meshes guide](https://docu.ngsolve.org/latest/netgen_tutorials/working_with_meshes.html)
[Netgen 2D geometry guide](https://docu.ngsolve.org/v6.2.2501/netgen_tutorials/define_2d_geometries.html)

## Domain names and physical materials

The VOL `materials` section stores **domain names**. FieldCraft writes neither physical
properties nor material laws into that section. The export command returns a `domains`
array with each domain's file number, file name, and naming source.

| Command reply source | How cells receive their domain name |
|---|---|
| `element_set` | The optional `materials` argument names element sets that cover the mesh. Each set name becomes a domain name, including when its elements have a section material. |
| `section_material` | With no explicit sets, a cell assigned a section uses that section's material name. |
| `cad_body` | An unassigned cell classified on a CAD body uses the body's name. |
| `unassigned` | An unassigned cell without a body uses `default`. |

For example, an unassigned cell in the CAD body `Bracket` appears as the domain `Bracket`.
Assigning a section with material `Aluminum` makes its automatic domain `Aluminum`.
Passing an element set named `Mount` explicitly gives its elements the domain `Mount`.
If an unassigned CAD body and a section material share a name, supply distinct element sets
to name their domains. Define the properties for those domains in the receiving analysis.

## Name lookup

The names `default` and `interface` are reserved for unclassified boundary facets.
Choose a distinct name for an explicit boundary set.

For 3D faces, use column 2 of `surfaceelements` to look up `bcnames`. Column 1 is a surface
descriptor. FieldCraft assigns a distinct descriptor to each boundary number and domain pair.
Several descriptors can share one boundary number, as when a named TEAM 13 interface touches
several regions.

For 2D edges, use column 1 of `edgesegmentsgi2` to look up `bcnames`. FieldCraft assigns a
distinct segment index to each name and domain pair, then repeats the name for each index.
The meshio reader uses these same cell-index columns for its codimension-one data and reads
each `bcnames` entry as an index followed by a name.
[meshio VOL cell indexing and name lookup](https://github.com/nschloe/meshio/blob/main/src/meshio/netgen/_netgen.py)
