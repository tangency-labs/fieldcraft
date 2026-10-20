# Electronics cooling examples

`make.py` solves the electronics-cooling samples in [docs/examples.md](../../docs/examples.md) and
regenerates their numbers and pictures in `docs/media/examples/` (and, with `--thumbnail`, their
Start-card thumbnails in `crates/ui-egui/thumbnails/`). The cases are engine samples,
`cooling-heat-sink`, `cooling-enclosure` and `cooling-aio`: their setup is
`crates/engine/src/sample/cooling.rs` and their geometry SolidCraft recipes
(`tests/corpus/cooling-heat-sink`, `tests/corpus/cooling-enclosure` and
`crates/engine/samples/cooling-aio*.json`). `make.py` opens a sample with `file.sample`, solves it,
prints the headline numbers and renders a still with a colour legend (and, for the AIO, a short MP4
and a still of the whole cooler).

```sh
python3 tools/cooling/make.py heat-sink                # build, solve, render into docs/media/examples
python3 tools/cooling/make.py enclosure --no-build     # with fieldcraft-cli already built (release)
python3 tools/cooling/make.py aio --backend gpu        # the radiator study on the GPU, with its clip
python3 tools/cooling/make.py heat-sink --no-render    # numbers only
```

Everything runs headless: the model is driven through `fieldcraft-cli mcp --headless` (the
`execute` tool for every command, the `screenshot` tool for the pictures), so no window or display
is needed. `FIELDCRAFT_GPU` picks the renderer and, with `--backend gpu`, the solver's GPU
(`software` for a CPU-only machine). To solve on another machine, point `--mcp` at a server there,
e.g. `--mcp "ssh gpu-box fieldcraft-cli mcp --headless"`; for the AIO's whole-cooler still give
`--design` a compiled `crates/engine/samples/cooling-aio.json` design it can read (the `compile_recipe` example
of `fieldcraft-geometry` writes one). Needs Python 3 (standard library only) and `ffmpeg` for the
legend and the MP4. Honours `CARGO_TARGET_DIR`.

The AIO geometry comes from `aio_recipe.py`, which writes `crates/engine/samples/cooling-aio.json`
(the whole cooler) and `cooling-aio-half.json` (the radiator cut at its middle plane, the half the
study solves). They stay out of `tests/corpus`, whose trunk test plans every part: these two would
add minutes to it. It builds the two coolant tubes as a true normal offset of one centreline and
refuses to write a recipe whose tubes would touch or fold. Rerun it after changing the cooler.

What the solver does today, and the stand-ins these cases use for what it does not (conjugate heat
transfer into solids #157, surface-to-surface radiation, a collision operator stable at fan speeds
#557), are described in docs/examples.md under *Electronics cooling studies*. The geometry is our
own SolidCraft recipes; nothing comes from a commercial product.
