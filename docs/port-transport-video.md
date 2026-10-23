# Inlet/outlet transport: in-app feature demonstration

[Final solve report](media/features/port-transport-report.json)

The 42-second recording and its still belong in attachments on
[PR #704](https://github.com/tangency-labs/fieldcraft/pull/704), outside git.
The reproduction inputs and report below describe the Results panel without requiring a video download.

This clip demonstrates port diagnostics, saved-frame playback, save/reopen without solving,
and display-unit conversion. It is a **single air-side transient demonstration**: air at 20 °C
flows past a bracket held at 80 °C. There is no coolant circuit or solved solid conduction. The
45 s state is not claimed steady or validated exchanger performance; the nonzero flow mismatch
and net advective heat are displayed. This is part of #701's per-PR evidence, not the completed
two-fluid release tutorial.

## Inputs and reproduction

The implementation is commit `fdc11c9` (`results: inspect inlet and outlet transport from saved
fluid frames`), on top of upstream develop `7d9c027`. The desktop binary still identifies itself
as 0.7.0; 0.7.5 is the target release, not a version already published.

1. Build the desktop app/CLI and the existing recipe compiler:

   ```sh
   cargo build -p fieldcraft -p fieldcraft-cli
   cargo run -p fieldcraft-geometry --features scripts --example compile_recipe -- tools/demo-video/port-transport-bracket.json port-transport-bracket.solidcraft
   ```

2. Replay `tools/demo-video/port-transport.json` through the app's control channel or MCP `batch`.
   `fieldcraft-cli run tools/demo-video/port-transport.json` reproduces the model and solution
   headlessly. The geometry path is relative to the working folder; give `geometry.open` its
   actual path if running elsewhere. The recipe is the sample's original `UPRIGHT_BRACKET`,
   also supplied separately so the design can be embedded in the saved project.
3. In the app, choose a light theme and an iso view, turn off the Max/Min markers, and expand
   **Results → Inlets and outlets**. The recorded inner window is 1440 × 1000 points. The
   temperature section is normal to y; its clamp is 20–80 °C, with coloured obstacles and
   streamlines. Show frames 0, 2, 4, 6, 8, 10 and 12. Each interval in this subset is 7.5 s of
   physical time and about 1.3 s of video: playback is accelerated and labelled.
4. Save with `file.save_as {path, geometry: "embed"}`, save a view with `view.save`, and
   `file.save`. Open another document with `file.new`, then `file.open` the saved model.
   Show `results.frame {step: "HeatedFlow", frame: 12}` and expand the port table again.
   The recorded `ports` report was checked equal before/after reopening, without a solve.
5. Run `model.units {system: "m-kg-s-N-Pa"}` to show kelvin, then return to
   `mm-t-s-N-MPa` with `mass_flow: "mg/s"` and `pressure: "Pa"` overrides. The stored fields
   and solve-time properties stay unchanged.

Every picture is the complete FieldCraft desktop window from `ui.screenshot`, captured through
the normal localhost control channel. Caption bands use `tools/demo-video/caption.py` with the
existing Inter SemiBold font. ffmpeg encodes H.264, 1440 × 1048, 10 fps, yuv420p, CRF 22 with
faststart. The still is a video frame. There are no synthetic solver fields or commercial assets.

## Run settings and scope

Windows, AMD Ryzen 9 7950X; CPU flow solver with four threads; viewport on NVIDIA GeForce
RTX 4090 through Vulkan. The development-profile LBM crate is optimised at level 3. This is a
small demonstration, not a throughput measurement.

The original sample domain is −30…100, −20…60, 0…70 mm at 5 mm spacing: 26 × 16 × 14 = 5824
cells before subtracting the bracket. It uses the sample's air properties, boundaries and
held-temperature obstacle faces. The supplied engine recipe changes only display units and the
step: Δt = 0.025 s, end time 45 s, 1800 steps, 12 output intervals (13 saved frames), CPU.
There is no steady stop rule or capacity acceleration. The final report includes no solver
warnings. Exact properties, relaxation times, Mach number, frame times and directional port
quantities are in the accompanying report.

The displayed final inlet/outlet mixed temperatures are 20.00/47.58 °C and the reference-density
mass rates are 6.720/6.358 mg/s. Advective sensible-heat rates are 135.1/304.0 mW relative to
0 °C. They are instantaneously sampled boundary-layer quantities: do not interpret their
difference as a complete energy balance or a converged heat duty. See
[the diagnostic conventions](port-transport.md) and
[the still-open coupled-exchanger contract](solvers/heat-exchanger.md).
