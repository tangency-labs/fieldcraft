# FieldCraft and AI agents (MCP)

FieldCraft speaks the [Model Context Protocol](https://modelcontextprotocol.io) over stdio
(JSON-RPC 2.0, one message per line; protocol revisions 2025-06-18, 2025-03-26 and 2024-11-05).
An agent builds and solves models the way a person does: every tool runs engine commands, the
same ones the UI, scripts, the command line and the control channel run, and never opens a
dialog or a file picker. The server is `crates/mcp` (`fieldcraft-mcp`); `fieldcraft-cli mcp`
runs it.

## Start it

```sh
# Headless (the default): an in-process session, no window; screenshots are drawn offscreen
# with the viewport's GPU pipeline (hardware or software GPU), or on the CPU without one.
fieldcraft-cli mcp
fieldcraft-cli mcp --headless

# Bridged to the running app, so you can watch the agent work:
fieldcraft --control 7070 &
fieldcraft-cli mcp --connect 127.0.0.1:7070
```

With `--connect`, an app that cannot be reached when the server starts gives a headless session
instead (a note goes to stderr). stdout carries only protocol messages; errors go to stderr.
Bridged screenshots are written to a temporary file by the app and read back by the server, so
both must run on the same computer.

### Claude Code

Use the CLI from your installation:

| Platform | CLI path |
|---|---|
| macOS, app installed in Applications | `/Applications/FieldCraft.app/Contents/MacOS/fieldcraft-cli` |
| Linux, unpacked tar.gz | `<unpacked-folder>/bin/fieldcraft-cli` |
| Linux, installed `.deb` or `.rpm` | `/usr/bin/fieldcraft-cli` |
| Windows, unpacked zip | `<unpacked-folder>\fieldcraft-cli.exe` |
| Built from source | `<repository>/target/release/fieldcraft-cli` (`fieldcraft-cli.exe` on Windows) |

The separate macOS CLI zip also contains `fieldcraft-cli` in its unpacked folder. Use an absolute
path in your MCP configuration; quote it if it contains spaces. For the installed macOS app:

```sh
claude mcp add fieldcraft -- /Applications/FieldCraft.app/Contents/MacOS/fieldcraft-cli mcp
```

```sh
claude mcp add fieldcraft -- /path/to/fieldcraft-cli mcp
# bridged to a running app:
claude mcp add fieldcraft -- /path/to/fieldcraft-cli mcp --connect 127.0.0.1:7070
# kept to one project folder:
claude mcp add fieldcraft -- /path/to/fieldcraft-cli mcp --automation-read-root /work/bracket --automation-write-root /work/bracket/out
```

When building from source, run `cargo build --release -p fieldcraft-cli` first. The same setup as a
project `.mcp.json`, or in another client's
configuration (Claude Desktop's `claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "fieldcraft": { "command": "/path/to/fieldcraft-cli", "args": ["mcp"] }
  }
}
```

## Confining file access

```sh
fieldcraft-cli mcp --automation-read-root /work/bracket --automation-write-root /work/bracket/out
```

The two flags keep the agent to the folders you give it:

- Every file a command reads must lie inside the **read root**: `open`, `file.open`,
  `geometry.open`, `geometry.import.mesh`, the design and meshes a `.fieldcraft` file opens
  again, `file.recovery.list` and `file.recover` with a `dir`.
- Every file a command writes must lie inside the **write root**: `save`, `file.save` and
  `file.save_as` (also a save to the model's own file), `write_deck` and `deck.write`,
  `screenshot {path}`, `file.recovery.discard` with a `dir`.
- Read and write access are separate. Give the same folder twice for both; a root left out
  grants none of its access (with only `--automation-read-root`, nothing can be written). Without
  either flag nothing is confined.
- A path is made absolute against the server's working folder, then its links are followed: a
  file (or folder) link that leads out of a root is outside it. For a file that does not exist
  yet, its nearest existing folder's links are followed and the rest must be plain names, so
  `..` below a folder that does not exist is refused, as is a link that leads nowhere (writing
  through it would create its target). The result is compared with the root by whole folder
  names: `/work2/x` is not inside `/work`. On Windows device names (`CON`, `NUL`, `COM1`…) are
  refused.
- A refused path is an ordinary tool error the agent can read and correct, for example
  `automation path rejected: outside the write root /work/bracket/out: /etc/x.inp` or
  `automation filesystem access is not granted: write authority is absent: /work/bracket/out/x.inp`.
  A design or mesh a model refers to outside the read root is not opened; the model opens
  without it and says so in `warnings`.
- The app's own recovery folder (autosaves, `file.recovery.*` without `dir`) is not confined:
  no agent names a path there.
- The roots must be existing folders; a missing one stops the server with an error. Each flag
  is given once (`--flag <dir>` or `--flag=<dir>`).

The flags confine the headless server and imply `--headless`; `--connect` with them is an error,
since the server would otherwise drive an app that is not confined.

A path is checked, then opened: another program that swaps a link in between can still win that
race, so the roots keep an agent to the files it was given rather than guarding against other
software on the computer.

## Bundled skill prompts

The server advertises `prompts: {"listChanged": false}`. A client can let the user choose
`fieldcraft-meshing` or `fieldcraft-modelling` from `prompts/list`, then retrieve it with
`prompts/get`. These names are stable and take no arguments. The complete small catalogue
is returned in one page; there is no continuation cursor.

```json
{"jsonrpc":"2.0","id":1,"method":"prompts/list"}
{"jsonrpc":"2.0","id":2,"method":"prompts/get","params":{"name":"fieldcraft-modelling"}}
```

The first returned user message contains the exact corresponding `skills/<name>/SKILL.md`
Markdown, including its frontmatter, embedded at compile time with `include_str!`. There is
no separately maintained prompt copy. Retrieving a prompt does not run commands or change
the model. Prompts are user-selected instructions, not automatic model actions.

The modelling prompt's remaining messages are three embedded resources. Their URI suffixes
map to its relative `references/` links:

| Relative reference | Fixed MCP resource URI |
|---|---|
| `references/checks.md` | `fieldcraft://skills/fieldcraft-modelling/references/checks.md` |
| `references/axial-rod.json` | `fieldcraft://skills/fieldcraft-modelling/references/axial-rod.json` |
| `references/steady-conductor.json` | `fieldcraft://skills/fieldcraft-modelling/references/steady-conductor.json` |

They also appear in `resources/list` and are readable with `resources/read`. These resources
contain the exact bundled reference/recipe text; recipes are not executed on retrieval.
Other repository documentation links still refer to the matching checkout. Both headless
and connected-app servers return the same compiled skill/reference contents, even without
a source checkout at runtime. Names and URIs are a fixed allowlist, with no arbitrary path
lookup. Unknown prompt names, nonempty or malformed arguments, and continuation cursors
return JSON-RPC invalid-params errors; unknown resource URIs retain resource-not-found errors.
The existing model/command resources are unchanged.

The wire format follows the [MCP prompt specification](https://modelcontextprotocol.io/specification/2025-06-18/server/prompts).

## Tools

| Tool | Arguments | What it does |
|---|---|---|
| `list_commands` | `filter?` | Every command with its label, parameter documentation and `schema` (its parameters as JSON Schema: types, units as `x-dimension`, defaults, required fields, one-of variants) and `requires` (what the session must hold first: `design`, `geometry`, `mesh`, `material`, `step`, `result`, …); `filter` matches id and label. |
| `execute` | `command`, `params?` | Run any command (`list_commands`). A failure leaves the model unchanged and names what to fix; a command whose `requires` are not met fails before doing anything, naming each missing piece and the command that provides it. |
| `batch` | `commands: [{command, params}]` | Run steps in order, stop at the first error; the reply has every result and the failing step (`failed.index`, `failed.command`, `failed.error`). Steps before it stay applied (`undo` takes them back). |
| `inspect_model` | — | `model.info` (counts, units, definitions), `model.activity` (the step table) and `geometry.bodies` (bodies with measures, persistent face names, which sets name as `{body, face}`, and persistent edge names `<face>|<face>` with lengths, which sets name as `{body, edge}`). |
| `check_model` | — | `model.check`: errors and warnings, each naming its item, and `progress`: the setup checklist with the next stage and its command. |
| `open` | `path`, `model?` | By extension: `.fieldcraft` (`file.open`), `.solidcraft` (`geometry.open`), `.3mf` (`geometry.import.mesh`). For a design, `model` is `keep` (default: the mesh and definitions stay, as a reopen after an edit wants) or `new` (a new model with the design); the reply's `related` is false when the design has none of the bodies the model refers to, and `stale` lists the mesh bodies and sets that no longer resolve. |
| `save` | `path?` | `file.save`: to `path` or the model's own file. |
| `write_deck` | `format` (abaqus, nastran, lsdyna), `path`, `steps?` | `deck.write`. |
| `mesh` | `method`, `params?` | Runs the mesher (`tet`: `mesh.tet.generate`; `hex`, `hex_sweep`, `hex_mapped`, `hex_blocks`: the structured hex meshers; `hex_extrude`: `mesh.hex.extrude`, layers from a section up to a target surface; `cells`: `mesh.cells`, a partitioned part cell by cell; `mirror`: `mesh.mirror`, a meshed half reflected across a plane; `block`: `mesh.block`; `planar`: `mesh.planar`) and replies with its result, `model.info` and `mesh.quality`. |
| `solve` | `analysis` (static, nonlinear, frequency, buckle, explicit, fluid, heat), `step?`, `options?` | `solve.<analysis>` with the step and the solver's caps and settings. |
| `show_results` | `step?`, `field?`, `bands?`, `range?`, `deformed?`, `scale?`, `mode?`, `frame?`, `time?` | The `results.*` commands its arguments ask for (`deformed: false` draws the undeformed mesh; `scale` is a factor, 0 for the undeformed shape so a sequence from 0 swings the deformation out, or `"auto"`; `range: {min?, max?}` clamps the legend, values above `max` drawn in the over colour); the reply has the contour range (smallest and largest value and their nodes: `field: "displacement"` gives the largest displacement), modes or frames. |
| `probe` | `step?`, `nodes?`, `at?`, `elements?`, `points?`, `forces?`, `mode?`, `frame?` | `solve.probe`: stored results at nodes, at the nodes nearest to points `at` (structural steps), elements (explicit) or points (fluid). |
| `reactions` | `step?`, `sets?`, `about?` | `results.reactions`: reaction force and moment per node set and in total. |
| `screenshot` | `view?`, `camera?`, `section?`, `source?` (window, model), `width?`, `height?`, `path?` | A PNG image content block and what it shows. Headless: an offscreen render with the viewport's GPU pipeline (`FIELDCRAFT_GPU=auto|hardware|software|cpu`; the reply's `renderer` names the tier and adapter, or why the CPU drew it) framing the model, the geometry, the analysis setup's symbols (boundary conditions, loads, reference points, couplings; `sampled_symbols` lists definitions drawn with a sample) and the results display (deformed shape, contour bands). Bridged: the whole app window drawn offscreen (`ui.capture`, the default) or a model render with the app's camera (`ui.render`). `view` turns the camera to a standard view first (framing everything drawn); `camera` aims it: `{yaw, pitch, distance?, target?, fit?, fov?}`, in degrees (yaw 0 looks along +Y, −45 is the iso view's; pitch is the elevation, clamped to ±90) and model units. `target` is `[x, y, z]` or `{"center_of": "mass" \| "bbox", "bodies"?: [names]}` (design bodies drawn; an unknown name is an error), by default the centre of mass of the model drawn (element faces, else the bodies), so a sequence of yaws orbits the object and keeps it at the picture's centre. The distance is `distance`; with `fit: true` the one at which the targeted bodies (or the whole model) fill the frame from this direction; without either, the one at which the whole model stays in view from every direction, so a sweep keeps its size. `fov` is the vertical field of view in degrees (0, the default, orthographic; up to 150). Give `view` or `camera`; the camera is kept for the next picture, and the reply's `camera` is the one drawn with (target, distance, angles in radians). Bridged, `camera` moves the app's camera the same way (`ui.camera`; the window is otherwise untouched). `section` cuts the mesh as `fieldcraft-cli snapshot --section` does: `x`, `y`, `z` (at the middle) or `AXIS=VALUE`, hiding the elements on the camera's side (headless only). `fieldcraft-cli snapshot --camera JSON` takes the same camera. |
| `select` | `filter?`, `box_mode?`, `mode?`, `clear?`, `geometry?`, `nodes?`, `elements?`, `faces?`, `at?`, `box?` | Bridged only: changes the app's viewport selection (`ui.select`; control-protocol.md, Selection) and returns it. Headless it is an error: selection belongs to the window. |
| `selection` | — | Bridged only: the viewport selection (`ui.selection`), with `params`, the selection as a set's members for `model.set.create`. |
| `undo`, `redo` | — | `edit.undo`, `edit.redo`. |

Resources: `fieldcraft://model` (`model.info`) and `fieldcraft://commands` (the command catalog).

Errors: an unknown tool, or arguments that do not match the tool's schema (a missing or unknown
argument, a wrong type, a value outside an enum), are JSON-RPC `-32602` errors; an engine failure
is a tool result with `isError: true` and the engine's message.

## Partition, mesh half, mirror

The commands that cut a part into cells and build its mesh from a half, all through `execute`
and `mesh` (each one's parameters: `list_commands`):

1. `geometry.datum {name, by}`: a datum plane — a principal plane moved along its normal
   (`{"type": "offset", "plane": "zx", "offset": 0}`), three points, or a point and a normal.
   `geometry.datum.axis` makes axes.
2. `geometry.partition {name, datum, bodies?}`: cuts the open design's bodies by the plane into
   cells `<body>/a` (negative side) and `<body>/b`. The two faces of the cut are
   `partition:<name>:a` and `partition:<name>:b`; every other face keeps the name of the face it
   is a piece of. The partition is part of the model: saving and opening it, or opening the
   design again after an edit, cuts the design again. `inspect_model` lists the cells.
3. `mesh` with method `hex` on an extruded cell needs no further cuts: its cross-section is cut
   into mappable blocks inside the mesher (O-grids round holes, triangle primitives at fillets,
   four-sided blocks), so a symmetric part takes one datum. Or `mesh` with method `hex`, `hex_sweep`, `hex_mapped` or `tet` on one cell (`{"body":
   "Part/b"}`), or method `cells` for all cells with a method each (`{"methods": {"Part/a":
   "tet"}, "default": "auto"}`): cells whose faces match on a partition are merged, others are
   tied (the reply lists each interface).
4. `mesh` with method `mirror` (`{"plane": "<datum>", "bodies": ["Part/b"]}`): the meshed half
   reflected across the plane, right-handed, joined on it, its elements and faces tagged as the
   other half's, so sets by face name (`{body, face}` on either half) reach both.
5. Or keep the half: `model.boundary_condition.symmetry {name, datum}` holds it on the plane
   (zero normal displacement), with half the load.

`file.sample {"name": "crank-mirror"}` opens the crank arm built this way: one datum (`Sym`), its
partition, `hex` on the half and `mirror` (docs/examples.md).

## Edge seeds

Seed an edge when its element count matters more than the global size: through a thin
section's thickness (at least 2–3 elements for bending), around a fillet or hole where the
stress peaks (a count or a small size, biased toward the stress raiser), or to give the shared
edges of partition cells matching counts.

1. `inspect_model` lists each body's edges by persistent name (`<face>|<face>`, `~k` when several
   separate the same faces) with their lengths; the names follow the faces through edits and
   partitions.
2. `execute` `mesh.seed.edge {edges: [{body, edge}], count | size, bias?: {ratio, toward:
   start|end|both}, constraint?: exact|min|free}`. `exact` keeps the count, `min` keeps at least
   it, `free` is a goal. The reply gives the count each edge gets.
3. Mesh: the hex meshers (`hex`, `hex_sweep`, `hex_mapped`, `hex_blocks`, `cells`) keep exact counts
   through their interval assignment and refuse, naming the seeds, when two exact counts cannot
   both hold (opposite sides of a mapped face, the layers between two levels); the tet mesher
   divides the edge at the seed's nodes and grades its size field from them (its reply's
   `seeds` gives the mesh edges along each). A seed on an edge two partition cells share holds
   in both, so their meshes still conform. A biased seed's spacing carries to the unseeded edges
   tied to it (the parallel edges of a mapped block or cell, the opposite side of a mapped face,
   the other side edges of a sweep), so one biased edge grades the whole layer stack.
4. `mesh.seed.list` lists the seeds and the edges they reach (`stale` when the design has no such
   edge); `mesh.seed.clear {edges?}` removes them. Seeds are saved with the model.

## Size controls

Ask for a local element size where the answer is decided (a fillet, a hole's wall, a contact
patch) instead of making the whole part fine.

1. `execute` `mesh.size.set {name, size, faces: [{body, face}] | bodies: [...] | box: {min, max} |
   sphere: {center, radius}}` adds the control (or replaces the one of that name). Faces are
   named as `inspect_model` lists them.
2. `mesh.size.growth {rate}` sets how fast the size grows away from a control (the largest
   change of element size per unit length; 0.3 by default, 0.05 to 2).
3. Mesh: every mesher honours the controls. The size at a point is the smallest of the mesh
   command's `size` and, over the controls, a control's size plus the growth rate times the
   distance to it. The tet mesher refines to that field; the structured hex meshers, the block
   projection and mesh plans divide their edges by it (a structured mesh carries a count across
   its block, so the refinement spans the row or layer); the octree refines its cells to it (a
   size that needs more levels than `max_depth` is refused before meshing). A control on whole
   bodies is their size in place of the command's.
4. `mesh.size.list` lists the controls, what each covers in the open design (`stale` when the
   design no longer has a face it names) and the growth rate; `mesh.size.clear {names?}` removes
   them. Controls are saved with the model.

## Meshing from the feature graph

A design opened from `.solidcraft` is rebuilt from its timeline, so how each body was made is
known. `execute` `geometry.features {"body"?}` lists, per body, the features that built it in
order: type, operation, evaluated values, placement and direction, the sketch profiles (lines,
arcs, circles), pattern and mirror transforms, targets and tools, and the persistent faces and
edges each one made. `mesh.plan {"body"}` turns that into a saved, editable hex-meshing plan
(templates per feature, O-grid and butterfly rings round circles, the partition cuts their frames
and imprints need, the copies a pattern or mirror made, and the features not mapped yet with why);
`mesh.plan.apply {"plan"}` cuts and meshes the body (and its copies) with the structured meshers
and reports `all_hex`, `conforming`, the smallest scaled Jacobian and what still needs the
projection mesher (`pending`).

## Relaxing a few bad elements

When `mesh.quality` shows a handful of spikes in an otherwise good mesh:

1. `execute` `mesh.node.relax {"around_element": <the scaled Jacobian's worst>}` (or `{"center":
   [x, y, z], "radius": r}`, or `nodes` / `elements` / `set`; `iterations`, `strength`, `falloff`).
   Nodes there relax (smart Laplacian, then optimization of the worst element), each move kept
   only when the worst element around the node does not get worse and nothing inverts. Nodes on a
   B-rep face slide in its UV parameters and stay on it, nodes on a B-rep edge slide along it,
   vertex nodes stay; nodes on a partition cut face move freely (both cells share them). The
   reply gives the nodes moved, the largest displacement, and the region's smallest scaled
   Jacobian, largest edge ratio and worst elements before and after. One undo step.
   For a whole region at once, `mesh.node.relax {"mode": "strong", "bodies": ["Base"], "below": 0.7}`
   morphs the mesh round every element below 0.7: their nodes move to better places and the
   nodes within `reach` (default three element lengths) follow by harmonic interpolation; edge
   and vertex nodes are held, face nodes stay on their faces, nothing outside the reach moves,
   and the smallest scaled Jacobian never drops. The reply lists each round's bad elements,
   clusters and region size (`rounds`). It runs as a job (`job.status`, `job.cancel`).
2. The relaxation is kept as a field of the model (`field`, default `Relax1`, …; saved with it).
   After a remesh, `mesh.node.relax.apply {"field": "Relax1"}` puts the relaxed shape on the new
   mesh (each node takes the displacement interpolated in the old cell it lies in, on its face,
   then one guarded smoothing sweep); with `"reapply": true` on the relaxation every remesh does
   it by itself (the mesher's reply lists `relaxations_applied`).

## Limits

The server treats its input as hostile. A message line may be at most 16 MiB (the control
channel's limit); a longer one is answered with `request too large` and skipped, and the server
keeps serving. A JSON-RPC batch holds at most 100 messages, a `batch` call at most 10 000
commands; JSON nested deeper than 128 levels is a parse error. Picture sides are clamped to
16–4096 pixels. Every command keeps its own caps (element counts, solver degrees of freedom and
factor size, iterations), and an engine failure or panic leaves the model as it was.

## A typical session

```text
→ {"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"me","version":"1"}}}
→ {"jsonrpc":"2.0","method":"notifications/initialized"}
→ {"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"open","arguments":{"path":"/work/steps.solidcraft"}}}
→ {"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"inspect_model","arguments":{}}}
← … "bodies": [{"name": "Steps", "faces": [{"name": "F1:-x"}, {"name": "F1:+x"}, …]}] …
→ {"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"mesh","arguments":{"method":"tet","params":{"body":"Steps","tolerance":0.5}}}}
→ {"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"batch","arguments":{"commands":[
     {"command":"model.units","params":{"system":"mm-t-s-N-MPa"}},
     {"command":"model.set.create","params":{"name":"all","geometry":[{"body":"Steps"}]}},
     {"command":"model.set.create","params":{"name":"root","geometry":[{"body":"Steps","face":"F1:-x"}]}},
     {"command":"model.set.create","params":{"name":"tip","geometry":[{"body":"Steps","face":"F1:+x"}]}},
     {"command":"model.material.create","params":{"name":"Steel","elastic":{"type":"isotropic","young":210000,"poisson":0.3},"density":7.85e-9}},
     {"command":"model.section.create","params":{"name":"solid","set":"all","material":"Steel","kind":{"type":"solid"}}},
     {"command":"model.boundary_condition.create","params":{"name":"clamp","type":"fixed","target":{"set":"root"}}},
     {"command":"model.load.create","params":{"name":"push","type":"force","target":{"set":"tip"},"vector":[0,0,-1000]}},
     {"command":"model.step.create","params":{"name":"Load","procedure":{"type":"static"},"boundary_conditions":["clamp"],"loads":["push"]}}]}}}
→ {"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"solve","arguments":{"analysis":"static","step":"Load"}}}
→ {"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"show_results","arguments":{"step":"Load","field":"displacement"}}}
← … "contour": {"field": "displacement", "range": {"max": {"node": 25, "value": 0.00485}, …}} …
→ {"jsonrpc":"2.0","id":8,"method":"tools/call","params":{"name":"screenshot","arguments":{"view":"iso"}}}
← {"jsonrpc":"2.0","id":8,"result":{"content":[{"type":"image","mimeType":"image/png","data":"iVBORw0…"},{"type":"text","text":"{\"source\":\"model\",…}"}]}}
```

`crates/mcp/tests/workflow.rs` runs this session on a part of the meshing corpus
(`tests/corpus/stepped-block`) and checks the results only through MCP; `crates/mcp/src/tests.rs`
sends malformed messages, huge parameters, paths outside the roots and unknown tools.
