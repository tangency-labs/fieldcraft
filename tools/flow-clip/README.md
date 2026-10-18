# Flow-past clip

`make.py` regenerates the gallery clip and still of the `flow-past` sample (docs/examples.md):
`docs/media/examples/flow-past-transient.mp4` (the start-up transient of air flowing past the
L-bracket, from the uniform 25 mm/s initial field to the settled wake) and
`docs/media/examples/flow-past-flow.png` (a clean iso render of the settled flow).

```sh
python3 tools/flow-clip/make.py                   # build, run, encode into docs/media/examples
python3 tools/flow-clip/make.py --no-build        # with the release binaries already built
python3 tools/flow-clip/make.py --frames /tmp/f   # also keep the rendered pictures
python3 tools/flow-clip/make.py --gif             # also write a GIF (not committed; larger)
```

Needs Python 3 (standard library only), a display for the FieldCraft window (the pictures are
drawn offscreen with `ui.render`, so a locked or virtual display is fine — on a headless machine
run it under `xvfb-run`), and `ffmpeg` with libx264 on the `PATH` (or `FFMPEG`). It honours
`CARGO_TARGET_DIR`; set `FIELDCRAFT_GPU=software` for a software renderer.

What it does:

1. Builds `fieldcraft` and `fieldcraft-cli` (release) and starts the app with `--control PORT`.
2. Opens the `flow-past` sample and, so the transient plays smoothly, edits its fluid step to
   store `FLOW_CLIP_FRAMES` frames (48 by default) over the 16 s run (`model.step.edit`), then
   solves the lattice Boltzmann step on the CPU (`solve.fluid`).
3. Fixes the velocity legend clamp (0 to 45 mm/s) so the colours mean the same thing from frame
   to frame, sets a light theme and an iso view, and renders the model-only viewport at each
   field frame (`results.frame` then `ui.render`): velocity on the mid-slice, the bracket
   coloured by the flow beside it, and streamlines seeded upstream, coloured by speed.
4. Captions each frame with the elapsed flow time in a band below, drawn by the demo-video tool's
   standard-library TrueType caption drawer (`tools/demo-video/caption.py`, Inter SemiBold).
5. Encodes the frames at 12 fps (the settled wake held at the end) into the MP4 (H.264) and writes
   the settled still. `--gif` also writes a GIF (800 px, its own palette); it is not committed
   because it is several times larger than the MP4.

The clip shows the wake forming behind the upright and the drag settling; it does not shed, which
needs a finer grid or a collision operator stable near τ⁺ = 0.5 (see issue #557). The geometry is
the original upright-bracket SolidCraft recipe (the `flow-past` sample); nothing comes from a
commercial product.
