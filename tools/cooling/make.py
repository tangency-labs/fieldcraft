#!/usr/bin/env python3
"""Solve the electronics-cooling samples and regenerate their gallery numbers and pictures.

The cases are engine samples (`crates/engine/src/sample/cooling.rs`, geometry in
`tests/corpus/cooling-*` and `crates/engine/samples/cooling-aio*.json`): this script opens one with `file.sample`, solves it, prints the
headline numbers (heat moved, temperatures, flow, effectiveness) and renders a still with a colour
legend (and, for the AIO, an MP4) into `docs/media/examples/`.

Everything runs headless: the model is driven through `fieldcraft-cli mcp --headless` (the `execute`
tool for every command, the `screenshot` tool for the pictures), so no window or display is needed.
`FIELDCRAFT_GPU` picks the renderer (and, with `--backend gpu`, the solver's GPU).

Python standard library only; ffmpeg draws the legend under the stills and encodes the MP4.
"""

import argparse
import base64
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "demo-video"))
import caption  # noqa: E402  (the demo-video caption drawer, standard library only)

# Copies of values set in Rust, for the derived numbers and the legend: keep them in step with
# crates/engine/src/sample/cooling.rs (the cases' settings) and crates/engine/src/view.rs (colours).
RHO_CP = 1.2e-12 * 1.005e9  # mJ/(mm^3 K): air's volumetric heat capacity
WIDTH, HEIGHT = 1280, 720
CAPTION_FONT = ROOT / "crates" / "ui-egui" / "fonts" / "Inter-SemiBold.ttf"
CAPTION_HEIGHT = 44
FPS = 12
# The caption bar under a picture.
BAR_UNDER = f"[0:v]pad={WIDTH}:{HEIGHT + CAPTION_HEIGHT}:0:0[w];[w][1:v]overlay=0:{HEIGHT}"
# The headless renderer's contour colours and band count (crates/engine/src/view.rs, results.rs).
CONTOUR = [(0x2B, 0x4F, 0xD8), (0x2F, 0xB8, 0xE0), (0x3C, 0xC8, 0x64), (0xF0, 0xD2, 0x3C), (0xE0, 0x40, 0x2E)]
BANDS = 12
POWER = {"W": 1.0, "mW": 1e-3, "kW": 1e3}
# The AIO's half shroud at its 0.6 mm cells, mm: keep in step with cooling::aio_shroud.
AIO_SPACING = 0.6
AIO_SHROUD = ([-100.0, 70.0, 4.8], [15.2, 150.4, 123.0])
# The graphics card (crates/engine/src/sample/gpu_card.rs, tools/cooling/gpu_card_recipe.py): the
# intake's speed (mm/s) and radius (mm), the powers (W), and where the sections go (mm).
CARD_INTAKE_SPEED = 2500.0
CARD_INTAKE_R = 30.0
CARD_POWER = {"Die": 275.0, "Memory": 20.0, "VRM": 25.0}
CARD_INTAKE_C = 38.0
CARD_MID_Y = 55.5            # the die's, the base's and the middle heat pipe's centre line
CARD_FIN_Z = (9.0, 36.0)     # the fins' height over the base


def ffmpeg(*args):
    exe = os.environ.get("FFMPEG", shutil.which("ffmpeg") or "ffmpeg")
    subprocess.run([exe, "-hide_banner", "-loglevel", "error", "-y", *map(str, args)], check=True)


def log(msg):
    print(f"cooling: {msg}", file=sys.stderr, flush=True)


def target_dir():
    return Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target"))


def build():
    cargo = os.environ.get("CARGO", "cargo")
    subprocess.run([cargo, "build", "--release", "-p", "fieldcraft-cli"], cwd=ROOT, check=True)


def compile_recipe(recipe, out):
    cargo = os.environ.get("CARGO", "cargo")
    subprocess.run(
        [cargo, "run", "-q", "--release", "--example", "compile_recipe", "-p", "fieldcraft-geometry",
         "--features", "scripts", "--", str(recipe), str(out)],
        cwd=ROOT, check=True,
    )
    if not Path(out).exists():
        raise SystemExit(f"compile_recipe wrote no design for {recipe}")


class Mcp:
    """A headless FieldCraft session over MCP (JSON-RPC on stdio)."""

    def __init__(self, command):
        self.proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        self.next_id = 0
        self.request("initialize", {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "cooling", "version": "1"}})
        self.proc.stdin.write(json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}) + "\n")
        self.proc.stdin.flush()

    def request(self, method, params):
        self.next_id += 1
        self.proc.stdin.write(json.dumps({"jsonrpc": "2.0", "id": self.next_id, "method": method, "params": params}) + "\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline()
        if not line:
            raise SystemExit(f"the MCP server closed its output ({method})")
        reply = json.loads(line)
        if "error" in reply:
            raise SystemExit(f"{method}: {reply['error']}")
        return reply["result"]

    def tool(self, name, arguments):
        result = self.request("tools/call", {"name": name, "arguments": arguments})
        content = result.get("content", [])
        if result.get("isError"):
            raise SystemExit(f"{name} {json.dumps(arguments)[:160]}: {' '.join(c.get('text', '') for c in content)}")
        return content

    def execute(self, command, params=None):
        content = self.tool("execute", {"command": command, "params": params or {}})
        return json.loads(content[0]["text"]) if content else {}

    def screenshot(self, path, camera):
        content = self.tool("screenshot", {"width": WIDTH, "height": HEIGHT, "camera": camera})
        png = next((c["data"] for c in content if c.get("type") == "image"), None)
        if png is None:
            raise SystemExit("screenshot: no image in the reply")
        Path(path).write_bytes(base64.b64decode(png))

    def close(self):
        self.proc.stdin.close()
        self.proc.wait(timeout=60)


def wall_heat_w(report, body):
    """Heat the body's walls give the air, in W (the report's power is positive into the wall, in
    the unit its `units` map names)."""
    unit = report.get("units", {}).get("power")
    if unit not in POWER:
        raise SystemExit(f"unexpected heat unit {unit!r}")
    return -sum(h["final"] for h in report.get("heat", []) if h.get("body") == body) * POWER[unit]


def common(report):
    return {
        "air_temperature_C": report.get("temperature"), "max_speed_mm_s": round(report.get("max_speed", 0.0), 1),
        "cells": report.get("cells"), "steps": report.get("steps"),
        "stop": {k: (report.get("stop") or {}).get(k) for k in ("reason", "time")} if report.get("stop") else None,
        "density_change": report.get("mass", {}).get("relative_change"), "warnings": report.get("warnings"),
    }


# ---------------------------------------------------------------------------- the cases
# Each case's settings as crates/engine/src/sample/cooling.rs sets them (speed mm/s, inlet or wall
# and hot temperatures C, the flow's cross-section mm^2), for the numbers derived from the report.

def exchanger(q, speed, area, cold, hot):
    """What the air could carry if it all left at the hot temperature, and the share it carries."""
    limit = RHO_CP * speed * area * (hot - cold) / 1000.0
    return {"airflow_limit_W": round(limit, 2), "effectiveness": round(q / limit, 3),
            "bulk_outlet_rise_K": round((hot - cold) * q / limit, 2)}


def heat_sink(report, cold, hot):
    """Mixed convection past a six-fin sink in a 44 x 30 mm duct at 0.1 m/s."""
    q = wall_heat_w(report, "HeatSink")
    return {"heat_removed_W": round(q, 3), "thermal_resistance_K_per_W": round((hot - cold) / q, 2),
            **exchanger(q, 100.0, 44.0 * 30.0, cold, hot), **common(report)}


def enclosure(report, cold, hot):
    """Natural convection: a chip on a PCB in a sealed box with held walls."""
    wetted = 528.0  # the chip's faces but the one on the board, mm^2
    q = wall_heat_w(report, "Chip")
    walls = -wall_heat_w(report, None)  # the box's walls take it back
    return {"heat_shed_W": round(q, 4), "chip_to_walls_K_per_W": round((hot - cold) / q, 1),
            "h_W_per_m2K": round(q / (wetted * 1e-6) / (hot - cold), 2), "walls_take_W": round(walls, 4), **common(report)}


def aio(report, cold, hot):
    """The half radiator: the numbers for the whole one (twice the half's)."""
    lo, hi = AIO_SHROUD
    area = (hi[0] - lo[0]) * (hi[2] - lo[2])
    half_q = wall_heat_w(report, "Radiator")
    half = exchanger(half_q, 100.0, area, cold, hot)
    return {"air_flow_l_per_min": round(2 * 100.0 * area * 60 / 1e6, 1), "heat_to_air_W": round(2 * half_q, 2),
            "airflow_limit_W": round(2 * half["airflow_limit_W"], 2), "effectiveness": half["effectiveness"], **common(report)}


def gpu_card(report, cold, hot):
    """The card: every solid's temperatures and the heat it gives the air, the airflow the
    intake sets, the die-to-intake resistance, and a hand check of the air's rise."""
    unit = report.get("units", {}).get("power")
    if unit not in POWER:
        raise SystemExit(f"unexpected heat unit {unit!r}")
    solids = {}
    for sr in report.get("solids", []):
        t = sr.get("temperature", {})
        solids[sr["body"]] = {"max_C": round(t.get("max", float("nan")), 1), "mean_C": round(t.get("mean", float("nan")), 1),
                              "heat_to_air_W": round((sr.get("heat_to_fluid", {}).get("final") or 0.0) * POWER[unit], 1)}
    power = sum(CARD_POWER.values())
    flow = CARD_INTAKE_SPEED * 3.141592653589793 * CARD_INTAKE_R ** 2  # mm^3/s
    die = solids.get("Die", {}).get("max_C", float("nan"))
    return {"solids": solids, "power_W": power, "airflow_l_per_s": round(flow / 1e6, 2), "airflow_cfm": round(flow / 1e6 * 2.11888, 1),
            "die_to_intake_K_per_W": round((die - CARD_INTAKE_C) / power, 3),
            "hand_check_air_rise_K": round(power * 1000.0 / (RHO_CP * flow), 1), **common(report)}


CASES = {
    "heat-sink": {
        "sample": "cooling-heat-sink", "numbers": heat_sink, "out": "cooling-heat-sink", "clip": False,
        "temperatures": [25.0, 50.0],
        "render": {"slice": {"axis": "y", "position": 20}, "streamlines": {"count": 24},
                   "camera": {"yaw": -40.0, "pitch": 22.0, "target": {"center_of": "bbox"}, "fit": True},
                   "caption": "heat sink · forced air 0.1 m/s, buoyancy on"},
    },
    "enclosure": {
        "sample": "cooling-enclosure", "numbers": enclosure, "out": "cooling-enclosure", "clip": False,
        "temperatures": [25.0, 40.0],
        "render": {"slice": {"axis": "x", "position": 35},
                   "streamlines": {"count": 14, "max_steps": 1400, "line": {"from": [35, 24, 30], "to": [35, 38, 30]}},
                   "camera": {"yaw": -90.0, "pitch": 6.0, "target": {"center_of": "bbox"}, "fit": True},
                   "caption": "sealed enclosure · natural convection"},
    },
    "aio": {
        # no clip: docs/media/examples/cooling-aio-airflow.mp4 is now recorded in the desktop app on a
        # longer domain (see its ATTRIBUTION.md row); a headless clip here would overwrite it.
        "sample": "cooling-aio", "numbers": aio, "out": "cooling-aio-airflow", "clip": False, "still": "cooling-aio",
        "temperatures": [25.0, 45.0],
        # the last layer of cells before the symmetry plane
        "render": {"surface": False, "slice": {"axis": "x", "position": round(AIO_SHROUD[1][0] - AIO_SPACING / 2, 4)},
                   "camera": {"yaw": -70.0, "pitch": 18.0, "fit": True},
                   "caption": "AIO radiator · fan air through the channels (cutaway)"},
    },
    "gpu-card": {
        # The steady field only (one frame): the clips move the section or the camera through it.
        "sample": "cooling-gpu-card", "numbers": gpu_card, "out": "cooling-gpu-card", "clip": False,
        "temperatures": [CARD_INTAKE_C, 95.0], "hide": ["Shroud"],
        "render": {"slice": {"axis": "z", "position": 22.0}, "streamlines": {"count": 40},
                   "camera": {"yaw": -38.0, "pitch": 30.0, "target": {"center_of": "bbox"}, "fit": True},
                   "caption": "graphics card · 320 W, blower air at 38 °C"},
        "sweeps": [
            {"out": "cooling-gpu-card-fins", "caption": "air through the fins · horizontal section",
             "slice": {"axis": "z"}, "positions": [CARD_FIN_Z[0] + (CARD_FIN_Z[1] - CARD_FIN_Z[0]) * k / 47 for k in range(48)],
             "streamlines": {"count": 40}, "camera": [{"yaw": -38.0, "pitch": 30.0, "target": {"center_of": "bbox"}, "fit": True}]},
            {"out": "cooling-gpu-card-section", "caption": "die, base, heat pipe and fins · vertical section",
             "slice": {"axis": "y", "position": CARD_MID_Y}, "positions": None,
             "camera": [{"yaw": -60.0 + 40.0 * k / 47, "pitch": 12.0, "target": {"center_of": "bbox"}, "fit": True} for k in range(48)]},
        ],
    },
}


# ---------------------------------------------------------------------------- render + encode

def legend_bar(font, title, lo, hi, unit, detail):
    """The caption bar: the title, the contour's bands from `lo` to `hi`, and a detail (the time)."""
    canvas = caption.Canvas(WIDTH, CAPTION_HEIGHT, caption.BAR)
    size, baseline = 18, CAPTION_HEIGHT / 2 + 0.73 * 18 / 2
    polygons, _ = caption.text_polygons(font, title, size, 16, baseline)
    canvas.paint(polygons, caption.TOOL)
    sw, sh = 22, 16
    left_label = f"{lo:g} {unit}"
    x0 = WIDTH - 16 - (caption.text_width(font, detail, size) + 28 if detail else 0) - caption.text_width(font, f"{hi:g} {unit}", size) - 8 - BANDS * sw
    x = x0 - 8 - caption.text_width(font, left_label, size)
    polygons, _ = caption.text_polygons(font, left_label, size, x, baseline)
    canvas.paint(polygons, caption.TEXT)
    top = (CAPTION_HEIGHT - sh) / 2
    for k in range(BANDS):
        t = (k + 0.5) / BANDS * (len(CONTOUR) - 1)
        i = min(int(t), len(CONTOUR) - 2)
        f = t - i
        color = tuple(round(a + (b - a) * f) for a, b in zip(CONTOUR[i], CONTOUR[i + 1]))
        xa = x0 + k * sw
        canvas.paint([[(xa, top), (xa + sw, top), (xa + sw, top + sh), (xa, top + sh)]], color)
    polygons, x = caption.text_polygons(font, f"{hi:g} {unit}", size, x0 + BANDS * sw + 8, baseline)
    canvas.paint(polygons, caption.TEXT)
    if detail:
        polygons, _ = caption.text_polygons(font, detail, size, x + 28, baseline)
        canvas.paint(polygons, caption.TEXT)
    return canvas.png()


def with_bar(shot, bar, out):
    """`shot` with `bar` under it, through ffmpeg."""
    ffmpeg("-i", shot, "-i", bar, "-filter_complex", BAR_UNDER, "-frames:v", "1", out)


def apply_view(mcp, render):
    mcp.execute("results.contour", {"field": render["field"], "range": {"min": render["clamp"][0], "max": render["clamp"][1]}})
    mcp.execute("results.surface", {"show": render.get("surface", True)})
    if render.get("slice"):
        mcp.execute("results.slice", {"field": render["field"], **render["slice"]})
    sl = render.get("streamlines")
    if sl:
        mcp.execute("results.streamlines", dict({"show": True}, **sl))


def thumbnail(shot, out):
    """The Start card's thumbnail: the middle of `shot` at 320 x 213."""
    ffmpeg("-i", shot, "-vf", f"crop={round(HEIGHT * 320 / 213)}:{HEIGHT},scale=320:213:flags=lanczos", "-frames:v", "1", out)


def render_still(mcp, render, font, work, out_png, thumb=None):
    apply_view(mcp, render)
    shot = work / "still.png"
    mcp.screenshot(shot, render["camera"])
    bar = work / "still-bar.png"
    bar.write_bytes(legend_bar(font, render["caption"], *render["clamp"], "°C", ""))
    with_bar(shot, bar, out_png)
    log(f"wrote still {out_png}")
    if thumb:
        thumbnail(shot, thumb)
        log(f"wrote thumbnail {thumb}")


def render_clip(mcp, report, render, font, frames_dir, mp4):
    times = report.get("frames", [])
    if len(times) <= 1:
        return None
    frames_dir.mkdir(parents=True, exist_ok=True)
    shots = []
    for k, t in enumerate(times):
        mcp.execute("results.frame", {"frame": k})
        shot = frames_dir / f"f{k:04d}.png"
        mcp.screenshot(shot, render["camera"])
        bar = frames_dir / f"c{k:04d}.png"
        bar.write_bytes(legend_bar(font, render["caption"], *render["clamp"], "°C", f"t = {t:4.2f} s"))
        shots.append((shot, bar))
    seq = frames_dir / "seq"
    (seq / "caps").mkdir(parents=True, exist_ok=True)
    n = 0
    for k, (shot, bar) in enumerate(shots):
        for _ in range(16 if k == len(shots) - 1 else 2):  # 6 fps of frames, the settled flow held at the end
            os.link(shot, seq / f"{n:05d}.png")
            os.link(bar, seq / "caps" / f"{n:05d}.png")
            n += 1
    mp4.parent.mkdir(parents=True, exist_ok=True)
    ffmpeg("-framerate", FPS, "-i", seq / "%05d.png", "-framerate", FPS, "-i", seq / "caps" / "%05d.png",
           "-filter_complex", BAR_UNDER, "-c:v", "libx264", "-preset", "veryslow", "-crf", "26", "-pix_fmt", "yuv420p",
           "-movflags", "+faststart", mp4)
    log(f"wrote clip {mp4} ({n} frames)")
    return mp4


def render_sweep(mcp, render, sweep, font, frames_dir, mp4):
    """A clip of the steady field: the section moving through `positions`, or the camera along
    its path, one picture per step, held at the end."""
    frames_dir.mkdir(parents=True, exist_ok=True)
    mcp.execute("results.streamlines", {"show": False})
    if sweep.get("streamlines"):
        mcp.execute("results.streamlines", dict({"show": True}, **sweep["streamlines"]))
    positions = sweep["positions"] or [sweep["slice"]["position"]]
    n = max(len(positions), len(sweep["camera"]))
    bar = frames_dir / "bar.png"
    bar.write_bytes(legend_bar(font, sweep["caption"], *render["clamp"], "°C", ""))
    seq = frames_dir / "seq"
    (seq / "caps").mkdir(parents=True, exist_ok=True)
    count = 0
    for k in range(n):
        mcp.execute("results.slice", {"field": render["field"], **sweep["slice"], "position": positions[min(k, len(positions) - 1)]})
        shot = frames_dir / f"f{k:04d}.png"
        mcp.screenshot(shot, sweep["camera"][min(k, len(sweep["camera"]) - 1)])
        for _ in range(FPS if k == n - 1 else 1):
            os.link(shot, seq / f"{count:05d}.png")
            os.link(bar, seq / "caps" / f"{count:05d}.png")
            count += 1
    mp4.parent.mkdir(parents=True, exist_ok=True)
    ffmpeg("-framerate", FPS, "-i", seq / "%05d.png", "-framerate", FPS, "-i", seq / "caps" / "%05d.png",
           "-filter_complex", BAR_UNDER, "-c:v", "libx264", "-preset", "veryslow", "-crf", "26", "-pix_fmt", "yuv420p",
           "-movflags", "+faststart", mp4)
    log(f"wrote clip {mp4} ({count} frames)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=sorted(CASES), help="which example to run")
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "media" / "examples")
    ap.add_argument("--no-build", action="store_true", help="use the release binaries already built")
    ap.add_argument("--no-clip", action="store_true", help="still only, skip the MP4")
    ap.add_argument("--no-render", action="store_true", help="numbers only")
    ap.add_argument("--backend", choices=["cpu", "gpu"], default="cpu", help="the fluid solver's backend (default cpu)")
    ap.add_argument("--thumbnail", action="store_true", help="also rewrite the sample's Start-card thumbnail "
                    "(crates/ui-egui/thumbnails/), from the same picture as the still")
    ap.add_argument("--frames-dir", type=Path, help="keep the clip's pictures in this folder")
    ap.add_argument("--mcp", help="the MCP server command (default: this build's `fieldcraft-cli mcp --headless`), "
                    "e.g. one run on a GPU machine through ssh or `kubectl exec -i`")
    ap.add_argument("--design", help="AIO: the whole cooler's compiled design for its still (a path the MCP server can "
                    "read) instead of compiling crates/engine/samples/cooling-aio.json")
    args = ap.parse_args()

    case = CASES[args.case]
    if not args.no_build and not args.mcp:
        build()
    work = Path(tempfile.mkdtemp(prefix="fieldcraft-cooling-"))
    args.out.mkdir(parents=True, exist_ok=True)
    font = caption.Font.load(CAPTION_FONT)
    mcp = Mcp(shlex.split(args.mcp) if args.mcp else [str(target_dir() / "release" / "fieldcraft-cli"), "mcp", "--headless"])
    try:
        if case.get("still") and not args.no_render:
            design = args.design
            if design is None:
                design = work / "cooler.solidcraft"
                compile_recipe(ROOT / "crates" / "engine" / "samples" / "cooling-aio.json", design)
            mcp.execute("geometry.open", {"path": str(design), "model": "new"})
            mcp.screenshot(args.out / f"{case['still']}.png", {"yaw": -52.0, "pitch": 24.0, "fit": True})
            log(f"wrote the still {case['still']}.png")
        mcp.execute("file.sample", {"name": case["sample"]})
        log(f"solving {case['sample']} ({args.backend}) ...")
        report = mcp.execute("solve.fluid", {"backend": args.backend})
        print(json.dumps(case["numbers"](report, *case["temperatures"]), indent=2), flush=True)
        if not args.no_render:
            render = {"field": "temperature", "clamp": case["temperatures"], **case["render"]}
            if case.get("hide"):
                mcp.execute("view.hide", {"items": [{"kind": "body", "name": b} for b in case["hide"]]})
            thumb = ROOT / "crates" / "ui-egui" / "thumbnails" / f"{case['sample']}.png" if args.thumbnail else None
            render_still(mcp, render, font, work, args.out / f"{case['out']}.png", thumb)
            if case["clip"] and not args.no_clip:
                render_clip(mcp, report, render, font, args.frames_dir or work / "frames", args.out / f"{case['out']}.mp4")
            for sweep in case.get("sweeps", []) if not args.no_clip else []:
                render_sweep(mcp, render, sweep, font, (args.frames_dir or work / "frames") / sweep["out"], args.out / f"{sweep['out']}.mp4")
    finally:
        mcp.close()
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
