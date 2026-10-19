# The `.fieldcraft` file (format `fieldcraft/2`)

A `.fieldcraft` file holds one model: its definition, its mesh, the geometry it is built on and
the commands that built it. It is written and read by `fieldcraft-io` (`write`, `read`; bytes in,
bytes out, so the browser build uses the same code) and by the commands `file.save`,
`file.save_as` and `file.open`.

## Package

The file is a ZIP package with exactly four parts, in this order, all stored uncompressed (so any
unzip tool reads them: `unzip -p bracket.fieldcraft model.json`):

| Part | Content |
|---|---|
| `fieldcraft.json` | Manifest: `format` (`"fieldcraft/2"`, [Versions](#versions-and-migration)), `generator` (`"FieldCraft 0.7.0"`), `geometry`, `view` (left out when empty), and `parts`: the byte size and SHA-256 of each other part. |
| `model.json` | The model definition: unit system, sets, materials, sections, reference points, couplings, boundary conditions, loads, steps, datum planes, axes and points, and the partitions that cut the geometry (in order; applied again when the design is opened). Pretty-printed JSON with a fixed key order (maps sorted by name), so a change shows as a line diff. Its form is the model's own JSON form, which commands take as parameters. |
| `mesh.bin` | Nodes, elements and geometry tags, binary (below). |
| `journal.json` | The commands that built the model, `[{"command": id, "params": {...}}, ...]`, one per line. It is a script `fieldcraft-cli run` replays. |

The same model gives the same bytes on every platform: entry times are fixed (1980-01-01),
permissions are fixed, and every map is written in sorted order. To diff two versions with git,
add `*.fieldcraft diff=fieldcraft` to `.gitattributes` and
`git config diff.fieldcraft.textconv "unzip -p"`; the definition and the journal show line by
line, the mesh part as binary.

## Geometry

```json
"geometry": {
  "design": {"path": "plate.solidcraft", "sha256": "…"},
  "meshes": [{"path": "/parts/tet.3mf", "sha256": "…"}]
}
```

`design` is the SolidCraft design (`geometry.open`), `meshes` the face-tagged 3MF files
(`geometry.import.mesh`) in import order. Each records the SHA-256 of the bytes that were opened.
Paths inside the model file's folder are stored relative to it (with `/`), others absolute.
Opening the model opens each file again: a missing, unreadable or changed file (its hash differs)
is reported as a warning and the model still opens; sets bind to geometry names, so a remesh
brings the mesh in line with a changed design. A geometry file that could not be opened keeps its
reference, so saving again does not lose it.

`view` is view state, not part of the definition: the items the viewport hides (`view.hide`,
`view.isolate`), as `{"hidden": [{"kind": "datum", "name": "Sym"}, …]}` with `kind` one of
`body`, `surface`, `datum`, `boundary_condition`, `load`, `reference_point`, `coupling`.

## `mesh.bin`

Little-endian:

```
"MCMESH01"                                    magic (the mesh encoding's version is 01)
u64 nodes, u64 elements, u64 entities, u64 tagged elements, u64 tagged nodes
nodes      × { u32 label, f64 x, f64 y, f64 z }             labels strictly ascending
elements   × { u32 label, u8 kind, u32 node × n(kind) }     labels strictly ascending
entities   × { u8 kind (0 body, 1 face, 2 edge, 3 vertex),
               u16 length, body name, [u16 length, entity name] unless a body }   UTF-8
tagged elements × { u32 element, u32 region (0xFFFFFFFF none),
                    u8 n, n × { u8 face, u32 tag }, u8 n, n × { u8 edge, u32 tag } }
tagged nodes × { u32 node, u32 tag }
```

Element kind codes: 1 rod2, 2 beam2, 3 tri3, 4 tri6, 5 quad4, 6 quad8, 7 tet4, 8 tet10, 9 wedge6,
10 hex8, 11 hex20, 12 rod3, 13 pyramid5. Node order is the model's (Abaqus/VTK). A tag indexes the entity list, which
lists each geometry entity once. Coordinates keep every bit (signed zero included).

## Reading a file safely

Every file is treated as hostile:

- The file is at most 8 GiB; the manifest, definition and journal parts at most 256 MiB each. No
  allocation is larger than the input.
- Only the four parts, each once. The ZIP CRC and the manifest's size and SHA-256 are checked for
  every part.
- The `format` is read first, before anything else in the file: a later `fieldcraft/N` is refused
  as a version error naming the program that wrote it, any value that is not `fieldcraft/<N>` as
  not a FieldCraft file. Unknown fields are errors (after migration).
- `mesh.bin`: counts are capped (50 million nodes and elements, one million entities) and checked
  against the bytes present before anything is read; labels must ascend; unknown kinds and
  trailing bytes are refused.
- The model is rebuilt through the same checked methods commands use: names, references,
  finite coordinates within ±1e9, element connectivity, tag dimensions. A file cannot hold what
  commands would refuse.

## Versions and migration

`format` is `fieldcraft/<N>`, N an integer that grows by one with every change to what a file
holds or how it is encoded. A part with an encoding of its own carries that encoding's version
too: `mesh.bin` starts with `MCMESH01`; a new mesh encoding is a new magic, brought in by a new
format.

| Format | Written by | Change |
|---|---|---|
| `fieldcraft/1` | FieldCraft 0.4.0 to 0.7.0 | the first format |
| `fieldcraft/2` | FieldCraft after 0.7.0 | versioning and migration; the parts hold the same as format 1 |

FieldCraft reads every format up to its own. A file of a later one is refused, before anything
else in it is read, with a message naming what wrote it: "this model was made by FieldCraft
0.9.0 (format 3); FieldCraft 0.8.0 reads formats up to 2. Update FieldCraft to open it."

An older file is migrated as it is read (`fieldcraft-io`, `migrate.rs`). The package is checked
as it is on disk (the four parts, their sizes and SHA-256); then the JSON parts are parsed as
plain JSON, `mesh.bin` kept as bytes, and the steps `N → N+1` run in order up to the current
format; only then are the parts read into the model, with every check a current file gets.

- A step is a pure function of the parts: the same file always gives the same model.
- **A step is never edited once a release has written its format.** Files that went through it
  are out there and later steps build on its output; a mistake is fixed by a new step and a new
  format.
- Every released format keeps fixtures: real files written by every release that wrote it, one
  folder per release (`crates/io/tests/formats/v<N>/<version>/`, never edited; format 1 has
  0.4.0, 0.5.0, 0.6.0 and 0.7.0), each with what that release's own reader gave for it. The tests
  open every one and check that the current reader gives that record, changed only as the
  migration steps say (compared as raw JSON), that each step does only what it says on the real
  files, and that saving the migrated model and reading it again gives the same model and the
  same bytes.

## Saving

`file.save` writes the bytes to a temporary file beside the target, flushes it, keeps the old
file as `<file>.bak`, renames the new file over the target and flushes the folder: a crash at any
moment leaves the old file or the new one. A path that exists and is not a regular file is
refused.

Opening a file never writes anything, whatever its format. A model opened from an older format
is marked as migrated (`migrated_from` in `file.open`'s reply; the document tab shows
"(from format N)"), and so is work recovered from an autosave of one. A save replaces a file
only when there is none or it is positively a model of the current format; anything else (an
older model it was migrated from, a later one, a damaged one or one whose manifest cannot be
read, or any other file a Save As picks) is first copied beside it, flushed to disk:
`<stem>.v<N>.fieldcraft` (`bracket.v1.fieldcraft`) for a model of format N, `<stem>.backup.<ext>`
for anything else. The reply reports it as `backup`. The copy never replaces a file: when the
name is taken, the next free one is used (`bracket.v1-2.fieldcraft`, …), and when no copy can be
made the save is refused and the original stays as it was. Later saves replace `<file>.bak` as
always, but not the copy. The browser build saves by download, which replaces nothing, so it keeps no copy.

## Autosave and recovery

The desktop app autosaves every 30 seconds while there are unsaved changes, into the recovery
folder: `FIELDCRAFT_RECOVERY_DIR` when set; else, in runs driven by automation (`--control`,
`FIELDCRAFT_CONTROL_PORT`) and on CI (`CI` set), a temporary folder
(`fieldcraft-automation/recovery` in the system's temporary folder, beside that run's recent files
and preferences), so they never add to the user's recoveries; else `FieldCraft/recovery` in the
per-user data folder. An entry is `<id>.fieldcraft`, `<id>.json` (the model's own path, time,
counts) and `<id>.lock`, held while the app runs. Saving (or undoing back to the saved state)
drops the autosave. An autosave is written only into the recovery folder, never over the
model's own file. After a crash, the next launch reports the entries; `file.recovery.list`
lists them newest first, `file.recover` opens one as unsaved work under its original path, and
`file.recovery.discard` deletes one or all of them (File › Recover… shows each with its name,
age and counts, with Recover, Discard and Discard all).

Entries prune themselves: each launch deletes left-behind entries older than 14 days, and those
whose model's own file was saved after the autosave (the file's modification time is later than
the entry's). Entries of running apps are never touched.
