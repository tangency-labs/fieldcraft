# Face-tagged meshes (3MF with B-rep face tags)

Besides opening a `.solidcraft` design, geometry can enter FieldCraft as a triangulated boundary
surface written by another program, with each triangle tagged by the B-rep face it lies on.
FieldCraft reads it into a `SurfaceMesh` (`fieldcraft_geometry::tagged_mesh`, command
`geometry.import.mesh`) without involving the kernel.

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
