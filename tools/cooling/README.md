# Electronics cooling examples

`make.py` solves the electronics-cooling samples in [docs/examples.md](../../docs/examples.md) and
regenerates their numbers and pictures in `docs/media/examples/` (and their Start-card thumbnails in
`crates/ui-egui/thumbnails/`). The cases are engine samples, starting with
`cooling-heat-sink`: their setup is `crates/engine/src/sample/cooling.rs` and their geometry the
`tests/corpus/cooling-*` SolidCraft recipes. `make.py` opens a sample with `file.sample`, solves it,
prints the headline numbers and renders a still with a colour legend.

```sh
python3 tools/cooling/make.py heat-sink                # build, solve, render into docs/media/examples
python3 tools/cooling/make.py heat-sink --no-build     # with fieldcraft-cli already built (release)
python3 tools/cooling/make.py heat-sink --no-render    # numbers only
```

Everything runs headless: the model is driven through `fieldcraft-cli mcp --headless` (the
`execute` tool for every command, the `screenshot` tool for the pictures), so no window or display
is needed. `FIELDCRAFT_GPU` picks the renderer and, with `--backend gpu`, the solver's GPU
(`software` for a CPU-only machine). To solve on another machine, point `--mcp` at a server there,
e.g. `--mcp "ssh gpu-box fieldcraft-cli mcp --headless"`. Needs Python 3 (standard library
only) and `ffmpeg` for the legend. Honours `CARGO_TARGET_DIR`.

What the solver does today, and the stand-ins these cases use for what it does not (conjugate heat
transfer into solids #157, a collision operator stable at fan speeds #557), are described in
docs/examples.md under *Electronics cooling studies*. The geometry is our
own SolidCraft recipes; nothing comes from a commercial product.
