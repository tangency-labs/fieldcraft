# Demo video

The inlet/outlet feature clip for 0.7.5 is documented in
[docs/port-transport-video.md](../../docs/port-transport-video.md). Its original geometry and
normal engine commands are `port-transport-bracket.json` and `port-transport.json`; the clip
records the complete app window, including the new Results table and save/reopen.

`make.py` regenerates the README's demo, `docs/media/ibeam-analysis.mp4` and
`docs/media/ibeam-analysis.gif`: an agent analyses an I-beam cantilever over MCP while the
FieldCraft desktop app shows each step.

```sh
python3 tools/demo-video/make.py                  # build, run, encode into docs/media
python3 tools/demo-video/make.py --no-build       # with the release binaries already built
python3 tools/demo-video/make.py --frames /tmp/f  # also keep the captured pictures
cd tools/demo-video && python3 -m unittest        # test the caption drawing
```

Needs Python 3 (standard library only), a display for the FieldCraft window (the pictures are drawn
offscreen, so a locked screen is fine) and `ffmpeg` with libx264 on the `PATH` (or `FFMPEG`).
It honours `CARGO_TARGET_DIR`.

What it does:

1. Builds `fieldcraft`, `fieldcraft-cli` and the testkit's `corpus` tool (release).
2. Writes the corpus cantilever I-beam (`tests/corpus/cantilever-i-beam`: 40 x 40 mm, 3 mm
   flanges, a 2 mm web upright, 200 mm long along x) as a `.solidcraft` design to a temporary
   folder with `corpus --part cantilever-i-beam --design`.
3. Starts FieldCraft (`fieldcraft`) with `--control PORT` and `fieldcraft-cli mcp --connect 127.0.0.1:PORT`, and
   speaks MCP (JSON-RPC over stdio) to it. Every step is an MCP tool call: `open` the design,
   `mesh` (`hex_sweep`, 4 mm), `batch` (steel, a solid section, the root clamp, a kinematic
   coupling of the free end at x = 200 mm to a reference point loaded with 1000 N down, −z), `check_model`, `solve`
   (`static`), and `show_results` for the von Mises contours on the undeformed shape
   (`scale: 0`) and then the deformed shape, its `scale` swung (eased) from 0 to 150, back to 0
   and out to 150 again. Every picture is the MCP `screenshot` tool with `source: window` (the
   whole window drawn offscreen, `ui.capture`): the setup in the `iso` view, the results through
   its `camera` argument (the iso orientation drawn back to 420 mm so the contour legend leaves
   the beam clear), and to close, a slow orbit of 60 pictures, one `camera` angle each, from the
   iso view half way round the beam. Each stage holds for 2 to 3 seconds.
4. Captions each stage with the MCP tool it calls and what it does (`open · tests/corpus/…`,
   `mesh · hex sweep 4 mm`, `solve · static`, …) in a band under the window. `caption.py`
   draws the bars in the standard library: a minimal TrueType reader (`cmap` format 4, `loca`,
   simple and composite `glyf` outlines, `hmtx` advances), a scanline rasteriser (non-zero
   winding, anti-aliased) and a PNG writer, with the interface's own font
   (`crates/ui-egui/fonts/Inter-SemiBold.ttf`, SIL OFL 1.1). ffmpeg pads each window picture
   with the band and lays the bar over it with its `overlay` filter.
5. Encodes the pictures at 10 frames per second, each held for its stage, into the MP4 (H.264,
   1280 x 848) and the GIF (960 px, its own palette), and prints the tip deflection and the peak
   von Mises stress.

The theme is FieldCraft's default (dark).

The script checks the solution against beam theory and stops when it is off. Tip deflection:
P L³ / 3 E I = 0.1429 mm in bending plus P L / G A_web = 0.0335 mm in shear, 0.1763 mm; the
solution gives 0.1764 mm (within 3 %). Root bending stress: P L c / I = 45.0 MPa
(I = 88 871 mm⁴); the solution's stress_xx on the top flange, fitted along the beam between one
depth from either end and extrapolated to the root, gives 44.8 MPa (within 2 %). At the root
itself the clamp, which stops the section's Poisson contraction, adds a local peak (62 MPa
stress_xx, 56 MPa von Mises) that beam theory does not have.
