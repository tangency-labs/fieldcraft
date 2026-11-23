# Face-tagged meshes (3MF with B-rep face tags)

Besides opening a SolidCraft design or a CAD file, geometry can enter FieldCraft as a triangulated boundary
surface written by another program, with each triangle tagged by the B-rep face it lies on.
FieldCraft reads it into a `SurfaceMesh` (`fieldcraft_geometry::tagged_mesh`, command
`geometry.import.mesh`) without involving the kernel.

An imported surface is triangles without a B-rep: it is shown, used for reference and as a fluid
obstacle (its face tags name the faces whose wall forces and heat are reported), and
`mesh.surface.check` checks it, but it is never volume-meshed. `mesh.tet.generate`,
`mesh.hex.generate`, `mesh.hex.sweep`, `mesh.hex.mapped`, `mesh.hex.octree`, `mesh.hex.project` and
`mesh.cells` need a design body and refuse a surface, tags or not, naming its file and the way out:
open the CAD file it came from (#1049). `geometry.stitch` closes its gaps under a tolerance, so that it bounds a
solid as an obstacle needs.

## Other triangle files

STL, OBJ, PLY and glTF files import through the same command, read by SolidCraft's mesh readers
(`fieldcraft_geometry::triangles`); which extensions these are comes from SolidCraft's format table
(`fieldcraft_geometry::formats`). Each mesh in the file becomes one surface with one face.

- **Units.** STL, OBJ and PLY state no unit: they are read in the command's `unit` (millimetres by
  default), as is a 3MF whose model states none. The unit is kept with the model (format 4), so
  opening it again reads the file alike. A 3MF that states its unit, and glTF (metres by its
  specification), are read in their own, and a `unit` for them is refused. The reply gives the
  unit read in and the bounding box's size in the display units, which is what shows a file read
  in the wrong unit.
- **STL.** A binary STL has no signature, so its size is its check: exactly 84 + 50·n bytes for
  the n triangles its header gives. Anything else that is not ASCII STL (it starts with `solid`
  and has facets) is refused, as is an ASCII facet without exactly three vertices, a file cut off
  inside a facet, and any coordinate SolidCraft's reader would leave out (not finite, or not below
  10⁹): a facet left out would leave the surface open where the file is damaged.
- **Left out.** SolidCraft's readers leave bad OBJ lines, PLY faces and glTF parts or placements
  past their caps out with only a warning; FieldCraft refuses such a file, naming what was left
  out. Warnings that leave no geometry out (a glTF's textures, skins and animations) are passed on
  in the reply.
- **Limits.** SolidCraft's readers cap the file size and the vertex and triangle counts, the
  result is held to the limits below, and the call runs under SolidCraft's panic guard, so a
  damaged file is an error.

## The format is SolidCraft's

SolidCraft defines and writes the format (**File → Export 3MF with Face Tags**, `file.export` with
`format: "3mf-tagged"`). The specification is
[`docs/formats/3mf-face-tags.md`](https://github.com/tangency-labs/solidcraft/blob/main/docs/formats/3mf-face-tags.md)
in the SolidCraft repository. In short:

- The file is an ordinary 3MF Core 1.3 package.
- Each `<triangle>` carries `ft:face` (its 0-based B-rep face index), and each `<object>` carries
  `ft:faces` (its face count).
- The namespace is `https://github.com/bherbruck/solidcraft/formats/3mf-face-tags/2026/10`
  (version 1). It is not listed in `requiredextensions`.
- Coordinates are written in shortest round-trip `f64` form, so they read back exactly.
- Within an object, either every triangle is tagged or none is.

FieldCraft follows that specification and does not extend it. Because the format is 3MF,
untagged 3MF from any other program goes through the same reader, with one face per object.

## What FieldCraft's reader does

- **Surfaces.** Every build item becomes one `TaggedSurface`: a `SurfaceMesh` with its name and
  face count. Lengths are converted to millimetres from `unit`. The surface is named after its
  object, or `Mesh<n>` when the object is unnamed.
- **Transforms.** An item's `transform` is applied, and a mirroring transform reverses the
  winding. A singular transform (one that flattens the object) is refused. Components are not read
  yet; the reader returns an error that says so.
- **Tags.**
  - Tags are matched by namespace URI, not by prefix.
  - An untagged object becomes one face (tag 0), and an unknown face-tag namespace reads as
    untagged.
  - An object that mixes tagged and untagged triangles is rejected, as is any tag not below
    `ft:faces`.
- **Limits.** All input is treated as untrusted:
  - The file is at most 512 MB and must be a regular file, read through the virtual file system.
  - Each part is at most 1 GB uncompressed, in at most 10 000 zip entries.
  - At most 20 million vertices and 20 million triangles, counted over every placement.
  - At most 10 000 build items, and XML nesting at most 256 deep.
  - Coordinates must be finite and within ±10⁹ mm (`fieldcraft_model::MAX_COORD`), both as read
    and once placed.
  - Every vertex index must exist.
- **Errors.** Any violation is an error that names the problem. Nothing panics.
