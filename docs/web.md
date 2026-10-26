# FieldCraft in the browser

`apps/fieldcraft-web` runs the desktop app's interface (`fieldcraft_ui_egui::ModelApp`) in a web
page through eframe's web runner, drawing with WebGPU where the browser has it and WebGL2
otherwise. The engine, the meshers and the solvers are the same code, compiled to
`wasm32-unknown-unknown`.

Live: <https://tangency-labs.github.io/fieldcraft/app/>, deployed with the documentation site from each release (a push to `main`,
`.github/workflows/pages.yml`).

## Run it locally

```sh
rustup target add wasm32-unknown-unknown
cargo install --locked trunk          # or a release binary from github.com/trunk-rs/trunk
cd apps/fieldcraft-web
trunk serve                           # http://127.0.0.1:8790, rebuilt on every change
trunk serve --release                 # as deployed (optimised, wasm-opt)
trunk build --release                 # the site in dist/web, as Pages gets it
```

Trunk fetches the matching `wasm-bindgen` and `wasm-opt` itself.

Address flags: `?sample` (or `?sample=cantilever`, `?sample=impact`) opens a sample; `?webgl`
draws with WebGL2 even where WebGPU is available.

CI keeps the crate compiling for the browser: `cargo xtask wasm` (part of `cargo xtask ci`)
checks `fieldcraft-web` with the library crates, and `pages.yml` builds the whole site with
Trunk on every merge to `develop` that touches code, without deploying it.

## What works

Everything that runs on a model: opening samples, opening SolidCraft designs and 3MF tagged
meshes, every mesher, materials, sections, steps, loads and boundary conditions, the implicit,
explicit and lattice Boltzmann solvers, the results display (contours, deformed shape, probe,
reaction sums), undo and redo, the command line and command search, and writing Abaqus, Nastran
and LS-DYNA decks. The solvers run on one thread in the browser; `fieldcraft-parallel` gives the
same numbers for any thread count (docs/adr/0001), so results are the desktop's.

The browser has no file system the app can reach, so files live in the page
(`fieldcraft_io::fs` keeps them in memory) and the host services move them in and out:

| Desktop | Browser |
|---|---|
| File › Open, a native picker | The browser's file picker; pick a model together with the design and meshes it refers to |
| Drag and drop | Drop files on the page; drop a model together with its geometry files |
| Save, Save As | The `.fieldcraft` file is downloaded (each save downloads it again) |
| Write Deck | The deck is downloaded |
| Autosave to the recovery folder, File › Recover… | Autosaves are kept in the browser's local storage; a page closed or crashed with unsaved changes offers them at the next visit (File › Recover…) |
| Leaving with unsaved changes asks first | The browser's "leave site?" prompt |
| Recent files | This visit's files only |
| Help › Check for Updates installs a new release (docs/updates.md) | A new deployment shows **Reload to update** in the status bar |

A `.fieldcraft` model refers to its design and meshes by path. Paths in the model's own folder
are stored relative to it, so picking or dropping the model with those files opens them again; a
model whose geometry is elsewhere opens without it and says so in the status bar.

## What doesn't

- **Long jobs block the page while they run.** The browser build has no worker threads: a mesher
  or a solve started from its form shows its progress first, then runs on the page's thread, and
  the page stops responding until it finishes (seconds for the samples, longer for large parts).
  Cancel does nothing once it runs. Worker threads need `SharedArrayBuffer`, which needs
  cross-origin isolation headers that GitHub Pages cannot send.
- **No control channel or MCP bridge**: a page cannot listen on a socket. Agents use the desktop
  app or `fieldcraft-cli mcp`.
- **Autosave space is small**: local storage holds a few megabytes per site, so the autosave of a
  large mesh may not fit; the status bar then says so, and saving (a download) keeps the work.
- **Two tabs share one store**: a tab offers the autosave of another FieldCraft tab that is still
  open as recoverable (recovering it opens a copy).
- **A lost GPU is not recovered**: if the browser takes the graphics device away, reload the page
  (unsaved work is in the autosave).
- **Download size**: about 26 MB of WebAssembly (10 MB compressed), loaded once per release.

## Requirements

A current Chrome, Edge, Firefox or Safari with WebGPU or WebGL2 and graphics acceleration turned
on. Without either the page says so instead of starting.
