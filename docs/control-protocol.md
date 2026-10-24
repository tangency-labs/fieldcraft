# Control protocol

`fieldcraft --control PORT` (or `FIELDCRAFT_CONTROL_PORT=PORT`) listens on `127.0.0.1:PORT` for
JSON lines. Each request is one line `{"id": …, "method": "…", "params": {…}}`; each reply is one
line `{"id": …, "ok": true, "result": …}` or `{"id": …, "ok": false, "error": "…"}`. Requests are
handled on the UI thread between frames; pointer and keyboard input is injected as real egui
events, one pointer event per frame. The protocol has the same shape as SolidCraft's.

```text
fieldcraft --sample --control 7070          # or --sample cantilever | impact
{"id":1,"method":"engine.execute","params":{"command":"model.info","params":{}}}
→ {"id":1,"ok":true,"result":{"elements":{"hex8":64,"tet4":162},"nodes":199}}
{"id":2,"method":"ui.render","params":{"path":"/tmp/model.png","width":900,"height":600}}
```

The channel is untrusted input. A request line may be at most 16 MiB; a longer one is answered
with `{"ok": false, "error": "request too large"}` and the connection is closed, as are
connections that send eight lines that are not JSON. At most 16 connections are served at once,
at most 64 queued requests run per frame so a flood cannot starve drawing, and a request the UI
thread has not answered within 60 seconds gets `{"ok": false, "error": "timeout"}` (a solve, a
mesher or an open is a job, answered when it ends however long it runs: see [Jobs](#jobs)). Numbers are
checked and clamped (drag steps, wheel lines, window size); `ui.render` refuses an image size out of range.

## Engine

| Method | Params | Result |
|---|---|---|
| `engine.execute` | `command`: command id, `params`: JSON object, `wait?` (default true) | the command's result; a solve, a mesher or an open runs as a [job](#jobs) |
| `job.status` | — | `job`: the running job (`title`, `command`, `elapsed` seconds, `progress` `{step: {name, index, count}, stage, done, total}`, `fraction` 0–1 when bounded, `line` as the window shows it), or null |
| `job.wait` | — | the running job's result when it ends (null `job` when none runs) |
| `job.cancel` | — | cancels the job (`cancelled`: what it was): its solver or mesher stops at its next increment, time step, supernode or refinement pass and what it made is dropped; the model stays as it was and every request waiting for it gets `{"ok": false, "cancelled": true}` |
| `engine.script` | `commands`: `[{command, params}, …]` | stops at the first error; `journal`: commands that changed the model so far |
| `engine.commands` | — | every command: `id`, `label`, `params` (the parameter documentation), `schema` (its parameters as JSON Schema) and `requires` (what the session must hold before it runs: `design`, `geometry`, `nodes`, `mesh`, `analysable`, `material`, `step`, `result`, …; `engine.execute` of a command missing one fails without doing anything, naming each missing piece and its toolbar entry) |

Commands never open forms or file pickers when run this way. `fieldcraft-cli commands` prints the
same registry. `file.sample {name}` opens a sample (`parts`, `cantilever`, `impact`) and
`solve.run {}` solves every step; the Results workspace comes up when a result appears (a static step with its contours on the undeformed shape: `results.deformed {}` shows the deformed shape). The
ribbon follows the commands run here, as it does for scripts, MCP, the command search and the
browser: a command that changes the model or the results display brings the tab where its
toolbar entry lives (Mesh for `mesh.*`, Properties, Steps, Loads & BCs… for `model.*`, Solve,
Results), as a click there would. Queries (`model.check`, `mesh.quality`…) and the Tools tab's
commands (files, undo, `model.info`) leave the ribbon where it is.

## UI

| Method | Params | Notes |
|---|---|---|
| `ui.inspect` | — | view settings (`ui`), camera, viewport rect (points), renderer (`gpu`/`cpu`) and adapter, open command form, `file_picker` (whether a native chooser is pending), status bar, node and element counts, the results display (`results`: what the `results.*` commands reply: view, `stale`, contour range, modes), the probed node (`probe`), `solve` (the title bar's Solve: `enabled`, and `unavailable` as `ui.toolbar`'s entries give it, listing everything `solve.run` misses), `busy` (the background job: a file or sample opening, a solve, a mesher run from its form by a click on OK), `job` (that job as `job.status` gives it), `prompt` (`{"unsaved": path}`, `{"design": path}` (a design of another part opened over a meshed or defined model), `"recover"`, `"samples"`), `recent` (recent files), `updates` (docs/updates.md: `automatic`, `checked_at`, `skip`, `available` (a newer version the status bar offers), `busy`, `staged`, `unsupported`; null when the run does not check), `selected` (the browser's selected definition `{kind, name}`, or seed `{kind: "seed", edge}`), `browser` (each top-level browser node's setup state: `[{node, state, tip, add}]`, `state` `done`, `next`, `blocked` (⚠: its command cannot run yet; `tip` says what it needs and where), `todo`, `optional` or null for nodes that are no setup stage, `add` its ghost rows `{stage, label, command, status}`; a ghost row's region is `layout`'s `browser/add/<stage>`), `nav` (the left-drag mode), `layout` (each region of the window drawn last frame as `[x, y, w, h]` in points: `app_bar`, `ribbon`, `workspace_switcher`, `document_tab`, `app_bar:<action>`, `browser`, `viewport`, `nav_bar`, `axis_view`, `selection_bar`, `legend` (while results are shown), `context`, `command_line`, `timeline`, `toolbar_dropdown`, `palette`), version |
| `ui.set` | any view setting: `style` (`shaded_edges`, `shaded`, `wireframe`), `show_nodes`, `show_geometry`, `show_mesh`, `show_loads` (boundary condition and load symbols), `show_interactions` (reference points, couplings, contact and tie surfaces, rigid walls), `symbol_size` (the setup symbols' size on screen in points, 6–48, at any zoom), `result_symbols` (the shown step's setup symbols over a results plot; off by default), `result_extremes` (markers at a contour's largest and smallest value), `perspective`, `dark`, `show_tree` (the browser), `show_timeline`, `nav_preset` (the mouse navigation: `fieldcraft`, `fusion`, `solidworks`, `inventor`, `alias`, `tinkercad`, `openscad`; an unknown id acts as `fieldcraft`) and `reverse_zoom` (both remembered per user), `workspace` (`prepare`, `model`, `analysis`, `results`), `tab` (the toolbar tab: `geometry`, `mesh`, `properties`, `assembly`, `steps`, `loads`, `interactions`, `solve`, `results`, `tools`; a tab of another workspace brings that workspace), `palette` (the command search open), `menu` (the open menu: `file`, `help`, `workspace`, `view`, `more`, or a panel's drop-down as `<tab>/<panel>` such as `mesh/Tet`; `null` closes it) | an unknown setting or value is an error; returns the settings |
| `ui.view` | `view`: `front`, `back`, `left`, `right`, `top`, `bottom`, `iso`, `home` (iso and fit), `fit` | every view frames the model as drawn (the deformed shape included; the setup's symbols keep their size on screen and are not framed) in the part of the viewport the overlays (`layout`'s `legend`, `selection_bar`, `axis_view`, `nav_bar`) leave free; `fit` keeps the direction of view. Opening another model, sample or design frames it in the home view (the same design opened again after an edit keeps the camera), and results that come up with part of the model out of sight are framed. Returns the camera |
| `ui.camera` | `yaw`, `pitch` (degrees), `distance?` (in the display units' length, like `target`'s point), `target?`, `fit?`, `fov?` (degrees) | turns the camera to any orientation about its target (with `target`, `[x, y, z]` or `{"center_of": "mass" \| "bbox", "bodies"?}` as MCP `screenshot`'s camera (docs/mcp.md), it moves to look at that target from `distance`, from the distance that fits the targeted bodies with `fit`, or from the one that keeps the whole model in view from every side; `fov` sets the perspective, 0 orthographic): `yaw` about Z (0 looks along +Y, the front view; −45 is the iso view's, −90 the right view's), `pitch` the elevation above the XY plane; `distance` moves the eye (the zoom of an orthographic view). Angles must be finite and the distance above 0: the yaw is wrapped, the pitch clamped to ±90 and the distance to 1 µm–10⁹ m; anything else is an error. A sequence of calls orbits smoothly. Returns the camera (`ui.inspect`'s: angles in radians, `distance` and `target` in the display units) |
| `ui.toolbar` | — | the ribbon: `workspace`, `workspaces`, `tab`, `tabs` (the workspace's), `menu` (the open menu), `panels` (each panel label's `menu` key and `rect`), and every catalogue entry (docs/features.md) with `id`, `label`, `tab`, `panel`, `command` (null while coming soon), `more` (its related commands), `coming_soon` (the milestone, or null), `enabled` (false while coming soon or while its command misses what it needs), `unavailable` (null, or what its command misses: `reason`, the tooltip's lines such as "Needs a material — Properties › Material", and `missing`, each with `requirement`, `what`, `provider` (the command that provides it) and `place` (its toolbar entry)), `fit` — `full` (a button), `collapsed` (its panel is one button), `more` (its panel is in the » menu) or `menu` (in its panel's drop-down) — and `rect` (screen points, last frame, on the shown tab; for `menu` entries the row while the drop-down is open) |
| `ui.start` | `command`: a command or toolbar entry id | like a click on the toolbar, a drop-down row or a search result: a command without parameters runs, one with parameters opens its form in the context panel; an entry that is coming soon is an error naming its milestone. Never opens a native file picker: file entries (open, save as, write deck) open their command's form |
| `ui.select` | `kind` (`set`, `material`, `section`, `step`, `boundary_condition`, `load`, `reference_point`, `coupling`, `contact`, `tie`, `fluid_domain`), `name`; nothing to deselect | selects a definition as a browser click does; its properties (with Edit and Delete) show in the context panel |
| `ui.select` (viewport) | any of `filter`, `box_mode`, `mode`, `clear`, `at`, `box`, `geometry`, `nodes`, `elements`, `faces` (and no `kind`) | changes the viewport's selection (see [Selection](#selection)); returns `ui.selection`'s reply |
| `ui.selection` | — | the viewport's selection: `filter`, `box_mode` (box selection on), `count`, `capped`, `summary`, `items`, `params` (the selection as a set's members), `hover` (what a click would pick under the pointer) |
| `ui.prompt` | `answer`: for `unsaved` `save`, `discard` or `cancel`; for `design` `new` (a new model with the design), `keep` (the design over the current model) or `cancel`; for `replace_mesh` `replace` (runs the form, replacing the mesh) or `cancel` (the form stays open); for `delete` (the item menu's Delete… or the Delete key: `{delete: question, kind, name, users, allowed}`) `delete` (with a datum's partitions when it has any; refused when `allowed` is false) or `cancel`; for `recover` `recover` (with `id` from `file.recovery.list`) or `close`; for `samples` `open` (with `name`) or `close`; for `crashed` (the start after a crash: `{crashed: [reports]}`) `report` (opens Report a Bug for the newest) or `close`; for `crash_reports` (Help › Crash Reports…, or `ui.start diagnostics.crash_reports`: `{crash_reports: [reports]}`) `report` (with `id`: opens Report a Bug for it), `delete` (with `id`), `delete_all` or `close`; for `report_bug` (`{report_bug: draft}`) `set` (any of `title`, `observed`, `share_session`, `include_model`; replies with the draft) or `close` — the browser opens only from a click (docs/crash-reports.md); for `consent` (the first run's question about anonymous usage) `share` or `decline`; for `settings` (Help › Settings…) `set` (any of `usage_sharing`, `usage_automation`; or `update_check`, `update_skip` (a version, `""` for none); or `units`, the display units for models without their own: a preset such as `"m-kg-s-N-Pa"` or `{"preset", "overrides"}`, docs/units.md) or `close` (docs/privacy.md, docs/updates.md); for `update` (a downloaded and verified update: `{update: state}`) `restart` (saves or asks about unsaved changes first, then closes the window to install) or `later` (installs when the app quits) | answers the window's question (`ui.inspect`'s `prompt`) |
| `ui.settings` | — | opens Help › Settings… (answered with `ui.prompt`); replies with the privacy state as `ui.inspect`'s `privacy` gives it: `available` (this run may share usage), `usage_sharing` (`null` until asked), `usage_automation`, `recording`, `queued` |
| `ui.item` | `kind` (`body`, `surface`, or a definition kind as `ui.select` takes it, `datum`, `partition`), `name`, `action`: `edit`, `rename`, `suppress`/`unsuppress` (a datum or partition), `hide`, `show`, `isolate`, `show_all`, `delete`, or `activate`/`deactivate` with `step` | an item menu's entry, as a right click on the item and a click on the entry: runs the same command (`rename` opens the `model.rename` form filled in; `delete` opens the delete prompt, answered with `ui.prompt`) |
| `ui.dialog` | `params?` (JSON object, fills the form), `use_selection?`, `run?`, `submit?`, `cancel?` | drives the open command form, in that order; `use_selection` is the form's Target row (see [Selection](#selection)); `run` returns the command's result and closes the form on success; `submit` is the form's OK button: a mesher that would replace a body's mesh first asks (`prompt` `replace_mesh`, with the question and what else goes: sets listing the old labels, solved results), and meshers run in the background (`busy`); a form opened by another form's **Create …** button (a reference field with nothing to pick) goes back to that form on `run` (with the new item picked) or `cancel` |
| `ui.click` | `x`, `y` (screen points), `button?`: `left`/`right`/`middle`, `shift?`, `ctrl?`, `alt?`, `double?` | real pointer input; in the viewport a click selects what is under it (Shift adds, Ctrl/Cmd adds or removes, a click on nothing clears the selection), and a click on a results plot also probes the nearest node of the face under it (`solve.probe`; `ui.inspect` shows the values); with the Relax Nodes form open (`mesh.node.relax`) a click on the mesh places its sphere of influence instead (Alt + wheel sizes it) |
| `ui.move` | `x`, `y` | |
| `ui.drag` | `x0`, `y0`, `x1`, `y1`, `button?`, `shift?`, `ctrl?`, `steps?` (1–200) | press, move, release: in the viewport the drag does what the navigation preset (`ui.set {nav_preset}`) binds; by default a left drag orbits about the model point under the press, or the centre of what is drawn from empty space (with box selection on, it selects the visible nodes, element faces or elements in the box; Ctrl + left from empty space still orbits), a right or middle drag (or Shift + left) pans |
| `ui.scroll` | `x`, `y`, `delta?` (wheel lines, −100–100; positive zooms in) | zooms toward the model point under `x`, `y` (the point at the target's depth over empty space), by the same factor per line in orthographic and perspective views; `ui.set {reverse_zoom}` turns it the other way |
| `ui.key` | `key` (egui key name, `Enter`, `Escape`…), `cmd?`, `shift?` | `F` over the viewport frames the model; `Escape` clears the selection (and then cancels an open form) |
| `ui.text` | `text` | typed text (into the focused field) |
| `ui.drop` | `paths`: 1–16 file paths | drops the files on the window: the next frame opens them as a drop from the file manager does, and as a file picked in the Open dialog does (a `.fieldcraft` model, a `.solidcraft` design, a `.3mf` tagged mesh; `ui.inspect`'s `busy` while they open). The native picker cannot be driven over the channel; this is the same route |
| `ui.worldToScreen` | `points`: `[[x, y, z]…]` (at most 1000, in the display units) | screen points in the viewport drawn last frame (null when not on screen) |
| `ui.resize` | `width`, `height` | window inner size in points (a maximized window is restored first) |
| `ui.fullscreen` | `on?` (default true) | into or out of full screen; on macOS the title bar's left group moves into the traffic lights' room while there are none |
| `ui.screenshot` | `path?` | PNG of the window, taken a few frames later; waits for a presented frame; when none comes within 8 seconds (a hidden window, a locked screen) it answers with `ui.capture`'s offscreen picture (`"offscreen": true`). Without `path`: `fieldcraft-shot.png` in the temporary folder |
| `ui.capture` | `path?`, `width?`, `height?` (points; default the window's), `scale?` (pixels per point, 0.5–4; default the window's) | PNG of the whole window (title bar, ribbon, open menus and search, panels, forms and the viewport) drawn offscreen, without fade-in animations: the UI pass runs in a separate egui context and is rasterised on the CPU, so it works with the window hidden, minimised or the screen locked. At most 8192 pixels a side and 16 megapixels. Without `path`: `fieldcraft-shot.png` in the temporary folder |
| `ui.render` | `path`, `width?` (default 1280), `height?` (default 800): whole pixels, each 1–8192 and at most 16 megapixels in all; `scale?` (pixels per point, 0.5–4, default the window's) | CPU render of the model with the current camera, view settings and selection highlight (no window or GPU needed). Symbols, lines and nodes are sized in points, so a render at the window's device resolution matches the window. A size out of range, or a scene over 20 million primitives, is an error |
| `app.quit` | — | closes the window |

## Selection

The viewport's selection belongs to the window, not the model: there is no engine command for
it, and nothing about it is saved or journalled. What is selected is kept in the model's
persistent terms, so it reads as a set definition: geometry by name (`{"body", "face"}`, or
`edge`, `vertex`, or the body alone; a face picked on the mesh is named through the mesh's
geometry tags, so it stays selected through a remesh) and the mesh by label (`nodes`,
`elements`, element `faces` as `{"element", "face"}`). Labels the model no longer has drop out
after an edit.

`ui.select` on the viewport takes, in this order:

- `filter`: what a pick selects: `faces`, `edges`, `vertices`, `bodies` (geometry), `nodes`,
  `element_faces`, `elements` (the mesh). The selection bar at the top of the viewport shows it.
- `box_mode`: `true` or `false` turns box selection on or off (a left drag draws the box).
- `mode`: how the items change the selection: `replace` (default), `add`, `remove`, `toggle`.
- `clear: true`: empty the selection first.
- `geometry`, `nodes`, `elements`, `faces`: items by name and label, as `model.set.create`
  takes them (at most 100 000 in each list; labels from 1).
- `at`: `[x, y]` in window points: pick there as a click does, with the filter (nothing under
  the point selects nothing; with `replace` it clears the selection).
- `box`: `[x0, y0, x1, y1]` in window points: the visible nodes, element faces or elements
  inside the box (other filters refuse a box). Faces in front hide what is behind them.

A selection holds at most 100 000 items; past that it is cut and `capped` is true. Picking
runs on the CPU against the triangles the viewport draws (on the GPU too), with a watertight
ray–triangle test: the nearest face wins, nodes, vertices and edges are the nearest on screen
within a few pixels that are not hidden, and ties go to the lower index, so the same pick
gives the same answer every time.

Command forms that apply to part of the model have a Target row with Use selection
(`ui.dialog {"use_selection": true}`): `model.set.create` and `model.set.edit` take the selection
as members; `model.boundary_condition.*`, `model.load.*`, `model.section.*` and
`model.coupling.*` get a set made from it (`model.set.create`, named after the definition:
`Clamp-set`, then `Clamp-set-2`…) in `target: {"set"}` or `set`. A right click in the viewport
offers Create set from selection, which opens `model.set.create`'s form filled in.

```text
{"id":1,"method":"ui.view","params":{"view":"right"}}
{"id":2,"method":"ui.select","params":{"filter":"element_faces","box":[300,120,1100,700]}}
→ {"id":2,"ok":true,"result":{"count":8,"params":{"faces":[{"element":20,"face":4},…]},…}}
{"id":3,"method":"ui.start","params":{"command":"model.set.create"}}
{"id":4,"method":"ui.dialog","params":{"params":{"name":"end"},"use_selection":true,"run":true}}
```

## Files

User clicks (title bar File menu, toolbar entries, Cmd/Ctrl+O, S, Shift+S), drops on the window and
files named on the command line go through native file pickers and open in the background
(`ui.inspect`'s `busy`); commands sent here take their paths as parameters and run at once
(`file.open`, `file.save_as`, `geometry.open`, `geometry.import.mesh`, `deck.write`).

## Jobs

`engine.execute` of a solve (`solve.run`, `solve.static`, `solve.nonlinear`, `solve.heat`,
`solve.frequency`, `solve.buckle`, `solve.explicit`, `solve.fluid`) or a mesher (`mesh.tet.generate`,
`mesh.hex.*`, `mesh.node.relax`, `mesh.plan.apply`, `mesh.cells`) or a partition (`geometry.partition`)
or an open (`file.open`, `file.recover`, `file.sample`, `geometry.open`, `geometry.import.mesh`)
runs it off the window thread, as the Solve button, a mesher form's OK and the Samples window do: the
window keeps drawing its progress (step, stage and count, elapsed time) with Cancel. The reply is
the command's result when it ends, however long that takes; with `"wait": false` it is
`{"job": n}` at once, and `job.status`, `job.wait` and `job.cancel` follow it. Cancel (the button
or `job.cancel`) interrupts the running command: every solver and mesher checks for it at its
natural points (increments and Newton iterations, time steps, fluid steps, Lanczos blocks, the
factorization's supernodes, the tet mesher's refinement and improvement passes, the hex meshers'
attempts) and stops with `cancelled`, leaving the model unchanged. One job runs at a
time; another long command meanwhile is refused, naming the job. The MCP server's `solve` and
`mesh` tools wait for the reply.

## Heated flow results

A solved fluid step with temperature in every saved frame offers `temperature` in the
Results field list and through `results.contour`, `results.slice`, `results.surface` and
`results.streamlines` (`{"field":"temperature"}`). Slices and obstacle surfaces sample
fluid-cell temperature; streamline colours sample temperature along the velocity traces.
The legend and existing `solve.probe` use the model's temperature scale, without conversion.
An isothermal step hides the field and rejects an explicit request with an explanation.
The frame slider and `results.frame` select the field's time, including the initial frame.

`results.surface {"field":"wall_heat_flux","frame":10}` shows the signed heat received
by each obstacle face divided by that face's complete tessellated geometric area. Positive
values heat the wall; a hot wall heating the fluid is negative. This is a face average,
not a resolved local distribution. The legend uses the model's consistent power/area unit:
for mm–tonne–second it is mW/mm² (not W/mm²); other systems follow their own units.

Wall heat uses the recorded face-group histories, linearly interpolated at the shown frame
only between recorded times. The initial frame precedes the history and remains uncoloured;
no heat is invented or extrapolated. Recorded adiabatic faces are zero. Unrecorded/imported
faces remain in their own colour. The command reply includes `wall_heat.faces` with body,
face, area, power and flux (null when unavailable), plus the sign and sampling convention.
This surface-only field suppresses the slice and leaves streamlines plain; selecting a cell
field restores the saved slice settings. `results.contour` also selects the wall field and
controls its bands/range. The desktop and headless picture use the same mapping and palette.
