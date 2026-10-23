# The `.fieldcraft` file (format `fieldcraft/3`)

A `.fieldcraft` file holds one model: its definition, its mesh, the geometry it is built on, the
commands that built it and, unless a save leaves them out, its solved results. It is written and read by `fieldcraft-io` (`write`, `read`; bytes in,
bytes out, so the browser build uses the same code) and by the commands `file.save`,
`file.save_as` and `file.open`.

## Package

The file is a ZIP package of four parts, then the optional ones, in this order, all stored
uncompressed (so any unzip tool reads them: `unzip -p bracket.fieldcraft model.json`):

| Part | Content |
|---|---|
| `fieldcraft.json` | Manifest: `format` (`"fieldcraft/3"`, [Versions](#versions-and-migration)), `generator` (`"FieldCraft 0.7.0"`), `geometry`, `view` (left out when empty), `run` (left out when empty), and `parts`: the byte size and SHA-256 of each other part. |
| `model.json` | The model definition, every quantity in SI (docs/units.md): its display units (`display_units`, left out when the model has none of its own), sets, materials, sections, reference points, couplings, boundary conditions, loads, steps, datum planes, axes and points, and the partitions that cut the geometry (in order; applied again when the design is opened). Pretty-printed JSON with a fixed key order (maps sorted by name), so a change shows as a line diff. Its form is the model's own JSON form, which commands take as parameters. |
| `mesh.bin` | Nodes, elements and geometry tags, binary (below). |
| `journal.json` | The commands that built the model, `{"units": system, "commands": [{"command": id, "params": {...}}, ...]}`, one command per line, as they were typed: `units` names the display system their bare numbers are in (a number typed in another system carries its unit, `{"value", "unit"}`); a migrated journal also names `model_units`, the system its model-level numbers were in. It is a script `fieldcraft-cli run` replays (a plain list of commands is a script in the session's units). |
| `geometry/<file>` | Optional: a copy of a geometry file the model is built on, when the file was saved with its geometry embedded (below). |
| `results.bin` | Optional: every solved step's latest result, with the hash of the definition it was solved from (below). `file.save` writes it unless given `results: false`. |

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

A file can also carry a copy of each geometry file (`file.save`/`file.save_as` with
`"geometry": "embed"`, File › Save With Geometry Embedded): each is a part
`geometry/<file name>` and its reference names it, `{"path": "plate.solidcraft", "sha256": "…",
"embedded": "geometry/plate.solidcraft"}`. The reference stays the source: the copy is opened only
when the file it names cannot be read (missing, unreadable or outside the folders an automation
client may read), with a warning saying so; a file that changed is still opened from its path, as
without a copy. A model opened from such a file is saved with its geometry embedded again unless a
save says `"geometry": "refer"`. References are the default: a packed file holds the design and
meshes themselves, so whoever gets it gets the geometry too (`file.save`'s reply lists what a file
holds in `contents`). A content shared by several references is one part; the copies are checked
like every part, must be the content their reference's SHA-256 names, and each part must be the
copy of a reference (no other part is accepted); a part is at most 512 MiB.

`view` is view state, not part of the definition: the items the viewport hides (`view.hide`,
`view.isolate`), as `{"hidden": [{"kind": "datum", "name": "Sym"}, …]}` with `kind` one of
`body`, `surface`, `datum`, `boundary_condition`, `load`, `reference_point`, `coupling`.

`view` also keeps the saved result views (`view.save`) and the one shown when the model was saved,
which opening it shows again (once its step has a result, from `results.bin` or a solve):

```json
"view": {
  "views": [{"name": "Cutaway", "units": "mm-t-s-N-MPa",
             "camera": {"yaw": -70, "pitch": 18, "fit": true},
             "results": {"step": "Flow", "frame": 10,
                         "contour": {"field": "temperature", "bands": 12, "range": {"min": 25, "max": 45}},
                         "slice": {"axis": "x", "position": 14.9},
                         "surface": {"show": false},
                         "streamlines": {"show": true, "count": 64, "max_steps": 4000}}}],
  "current": "Cutaway"
}
```

Each part is the parameters of the command that sets it: `camera` is `ui.camera`'s (and the MCP
`screenshot` tool's), `deformed`, `contour`, `slice`, `surface` and `streamlines` are those
`results.*` commands' (without `step`), `mode` and `frame` the number `results.modes` and
`results.frame` take. `units` names the system their bare numbers are in, as a journal's header
does: they are read in it whatever the model is shown in. Showing a view (`view.restore`, or
opening the file) runs those commands in turn and aims the app's viewport (or the next MCP
picture) at the camera; what no longer applies (a field the step has not) comes back as a
warning. At most 1000 views, names 1 to 64 characters; `current` names one of them. The view
state is outside the definition: it is not in `model.json` and not in a result's definition hash.

`view` keeps the renders too (`results.render.save`): how a saved view becomes a clip.

```json
"renders": [{"name": "Airflow clip", "view": "Cutaway", "step": "Flow", "frames": "all",
             "fps": 12, "hold": {"each": 2, "last": 16}, "size": [1280, 720],
             "caption": "AIO radiator · fan air through the channels", "legend": true,
             "time_label": true, "format": "png-sequence"}]
```

`view` names a saved view; `step` the step whose frames play (the view's when left out); `frames`
is `"all"` or `[from, to]` (from 0, inclusive); each picture shows for `hold.each` frames of the
clip and the last for `hold.last`, at `fps` (1 to 120); `size` is the picture in pixels (16 to
4096 a side), with a caption bar under it when there is a `caption`, a `legend` or a
`time_label`. Left out, a field takes these defaults: all frames, 12 fps, holds of 1, 1280 × 720,
legend and time label on. `results.export.animation {"render", "path"}` draws it offscreen, picture
by picture, with the host's GPU pipeline (the CPU rasterizer without one), and writes the pictures
beside `path` as `<stem>_00000.png`, `<stem>_00001.png`, … (each repeated for its hold); the only
`format` is `png-sequence`: FieldCraft bundles no video encoder, and the reply's `ffmpeg` is the
command line that encodes the sequence as `path` (an `.mp4` or `.gif`). At most 1000 renders, each
starting from a saved view (a view a render uses is not deleted).

`run` keeps each step's preferred run settings (`solve.settings`), a preference outside the
definition (not in `model.json`, not in a result's definition hash):

```json
"run": {"Flow": {"backend": "gpu", "fallback": "cpu"}}
```

`backend` is `gpu` or `cpu`, `fallback` (`cpu`, or left out to refuse) what a GPU preference
does on a machine with no suitable GPU: `solve.fluid` and `solve.run` without a `backend` of their
own follow the preference, falling back with a `notice` in the reply (and the Solve dialog's
note). A machine overrides every preference with `FIELDCRAFT_SOLVE_BACKEND=cpu|gpu`, a call with
its own `backend`. Only fluid steps choose a backend today; settings of steps the model no longer
has are not saved.

## `mesh.bin`

Little-endian:

```
"MCMESH02"                                    magic (the mesh encoding's version is 02: coordinates in metres)
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

## `results.bin`

The results part keeps what each solved step's solver reported, so a model opens showing its
results without solving again (`results.contour`, `results.frame`, `results.slice`,
`results.streamlines`, `solve.probe` and the rest work at once). Little-endian:

```
"MCRSLT01"                     magic (the results encoding's version is 01)
u32 steps
steps × { u64 length, step name (UTF-8),
          32 bytes: SHA-256 of the definition it was solved from,
          u8 kind (0 static, 1 frequency, 2 heat transfer, 3 explicit dynamic,
                   4 legacy fluid, 5 fluid with port properties),
          u64 length, the result }                                     names ascending
```

The result is the solver's own, every number in SI and every bit kept: a static step's
displacements, rotations, reactions, nodal, element and integration point stresses, shell and
beam forces, reference points, contacts, and a nonlinear step's increments with a frame per
converged increment; a frequency or buckling step's modes; a heat transfer step's temperature
frames, fluxes and heat flows; an explicit step's field frames, histories, energies, contacts and
eroded elements; a fluid step's lattice (origin, spacing and cells), its cell fields per frame
(velocity, pressure, temperature), wall forces, heat and fan and porous zone histories, run
statistics and steady-state monitors. It is written in a compact binary encoding
(`fieldcraft_io::codec`): fixed-width numbers, a `u64` count before a list or text, one byte before
an optional value and before each field of a structure, a `u32` before a variant. The encoding
follows the order of the solvers' result types, so it has a version of its own in the magic: a
change to those types is a new version. A file whose results are of another version opens without
them, with a warning to solve again: results are what solving gives again, the definition is not.
Kind 5 wraps the unchanged kind-4 fluid payload followed by optional solve-time port properties
(density, optional specific heat). Kind 4 remains readable and has no port properties: mass and
sensible-heat transport are unavailable, rather than inferred from an edited material. Adding a
new kind leaves existing payloads intact; an older reader refuses the unknown results kind and
opens the definition with its usual warning. The result's in-memory `port_properties` field is
skipped by the legacy payload codec and is supplied only by this envelope.

The definition hash is the SHA-256 of what solving the step reads (the mesh, the definitions the
step uses and the step itself, as their JSON) and of the hashes of the geometry files the model is
built on. On opening, a result whose hash is the model's is current; any other is kept and shown
as stale, with a warning (`file.open`'s `results` lists each step with `stale`), and an edit that
brings the model back to its definition makes it current again. Editing the model marks results
stale as it does in a session. A solve is an unsaved change, since the file keeps results;
autosaves, crash reports and bug reports keep none.

## Reading a file safely

Every file is treated as hostile:

- The file is at most 8 GiB; the manifest, definition and journal parts at most 256 MiB each. No
  allocation is larger than the input.
- Only the parts the manifest lists, each once: the four, and `results.bin` and `geometry/<file>`
  parts when they are listed (format 3; older formats list exactly the four). The ZIP CRC and the manifest's size and
  SHA-256 are checked for every part.
- The `format` is read first, before anything else in the file: a later `fieldcraft/N` is refused
  as a version error naming the program that wrote it, any value that is not `fieldcraft/<N>` as
  not a FieldCraft file. Unknown fields are errors (after migration).
- `mesh.bin`: counts are capped (50 million nodes and elements, one million entities) and checked
  against the bytes present before anything is read; labels must ascend; unknown kinds and
  trailing bytes are refused.
- `results.bin`: at most 10 000 steps; every length and count is checked against the bytes left
  before anything is read for it (each element takes at least one byte, so nothing is allocated
  beyond a small multiple of the part); unknown kinds, malformed values and trailing bytes are
  refused; every result is checked before it is used: its step name, ascending node and element
  labels, per node, element and cell arrays the lengths of what they follow, a fluid grid of 1 to
  2³² − 1 cells with a finite origin and a positive spacing. A file is at most 8 GiB in all; a save
  that would write more is refused, naming `results: false`.
- The model is rebuilt through the same checked methods commands use: names, references,
  finite coordinates within ±1e9, element connectivity, tag dimensions. A file cannot hold what
  commands would refuse.

## Versions and migration

`format` is `fieldcraft/<N>`, N an integer that grows by one with every change to what a file
holds or how it is encoded. A part with an encoding of its own carries that encoding's version
too: `mesh.bin` starts with `MCMESH02`; a new mesh encoding is a new magic, brought in by a new
format.

| Format | Written by | Change |
|---|---|---|
| `fieldcraft/1` | FieldCraft 0.4.0 to 0.7.0 | the first format |
| `fieldcraft/2` | FieldCraft after 0.7.0 | versioning and migration; the parts hold the same as format 1 |
| `fieldcraft/3` | FieldCraft after 0.7.0 | SI storage (#579): `model.json` in SI with `display_units` in place of `units`, `mesh.bin` `MCMESH02` (coordinates in metres), the journal with a header naming its units. Step 2 → 3 scales by the file's declared `units` (`unspecified` as mm-t-s-N-MPa) with its own frozen field map and MCMESH01 codec; the journal's numbers are kept as typed. Before its first release it also gained optional parts and fields, each absent in a file without it (so every format 3 file written before them reads as it did): `results.bin` (#644), saved views (#645) and renders (#646) in `view`, embedded geometry parts (#647), and `run` in the manifest (#648). |

Format 3 was not yet written by a release when its optional parts were added, and every addition
is optional (absent means what a file without it always meant), so they are part of format 3
rather than a format 4 with a step that would only change the tag: no released FieldCraft reads
format 3, and every format 3 file opens unchanged.

FieldCraft reads every format up to its own. A file of a later one is refused, before anything
else in it is read, with a message naming what wrote it: "this model was made by FieldCraft
0.9.0 (format 4); FieldCraft 0.8.0 reads formats up to 3. Update FieldCraft to open it."

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

Nothing else FieldCraft writes where a user points it (decks, STL, VTU and CSV exports, XY
tables, pictures from the control channel, MCP and `fieldcraft-cli snapshot`) replaces a
FieldCraft model or a SolidCraft design: a target that is one, by its extension or by its
content, is refused with "choose another file name" (`fieldcraft_io::fs::write_output`).

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
