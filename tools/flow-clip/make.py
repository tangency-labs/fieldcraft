#!/usr/bin/env python3
"""Regenerate the flow-past gallery clip and still image: the start-up transient of air flowing
past the L-bracket, from the uniform 25 mm/s initial field to the settled wake.

The FieldCraft desktop app opens the `flow-past` sample, solves the lattice Boltzmann step with a
dense frame output, and draws the flow (velocity on a mid-slice, the bracket coloured by the flow
beside it, streamlines seeded upstream). Each frame of the step is rendered through the control
channel (`ui.render`, the model only, no window furniture), with a fixed legend clamp so the
colours mean the same thing from frame to frame. ffmpeg encodes the frames as `docs/media/examples/
flow-past-transient.mp4`, captioned with the elapsed flow time, and a clean iso render of the
settled flow is written as `docs/media/examples/flow-past-flow.png`. `--gif` also writes a GIF (not
committed — it is several times larger than the MP4; use it for hosting the clip elsewhere).

Python standard library only; ffmpeg encodes the frames. caption.py (the demo-video tool's
standard-library TrueType caption drawer) draws the time bars.
"""

import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "demo-video"))
import caption  # noqa: E402  (the demo-video caption drawer, standard library only)

# The flow-past sample: air at 25 mm/s past the upright L-bracket (Reynolds number 83 on the
# 50 mm upright). The sample solves 16 s (about 1.7 passes through the 240 mm box at 25 mm/s); the
# clip plays that start-up transient. FRAMES field frames are stored (the sample's 10 is too coarse
# for a smooth clip), so the clip holds FRAMES + 1 pictures from the uniform initial field (frame 0,
# 25 mm/s everywhere) to the settled flow.
SAMPLE = "flow-past"
STEP = "Flow"
DOMAIN = "Wind"
END_TIME = float(os.environ.get("FLOW_CLIP_END_TIME", 16.0))
FRAMES = int(os.environ.get("FLOW_CLIP_FRAMES", 48))

# The render: a clean iso view of the domain, each frame the model only (no window chrome). The
# legend clamp is fixed over the whole run so the wake's colours are comparable frame to frame;
# 45 mm/s sits just above the settled peak (about 42 mm/s), 0 for the still fluid in the wake.
WIDTH, HEIGHT = 1280, 720
CLAMP = (0.0, 45.0)
CAMERA = {"yaw": -40.0, "pitch": 24.0, "target": {"center_of": "bbox"}, "fit": True}

# The video: frames per second, and how long the settled flow is held at the end.
FPS = 12
HOLD_END = 18  # frames (1.5 s) on the settled wake

# The caption bar, a band under each picture (so it never hides the flow).
CAPTION_FONT = ROOT / "crates" / "ui-egui" / "fonts" / "Inter-SemiBold.ttf"
CAPTION_HEIGHT = 44


def log(msg):
    print(f"flow-clip: {msg}", file=sys.stderr, flush=True)


def target_dir():
    return Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target"))


def build():
    cargo = os.environ.get("CARGO", "cargo")
    subprocess.run([cargo, "build", "--release", "-p", "fieldcraft", "-p", "fieldcraft-cli"], cwd=ROOT, check=True)


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class Control:
    """A client of the FieldCraft control channel (one JSON line per request and reply)."""

    def __init__(self, port, timeout=180.0):
        end = time.monotonic() + timeout
        while True:
            try:
                self.sock = socket.create_connection(("127.0.0.1", port), timeout=5.0)
                break
            except OSError:
                if time.monotonic() > end:
                    raise SystemExit(f"the control channel on port {port} did not open within {timeout:.0f} s")
                time.sleep(0.5)
        self.sock.settimeout(600.0)
        self.file = self.sock.makefile("rwb")
        self.next_id = 0

    def send(self, method, params):
        self.next_id += 1
        self.file.write((json.dumps({"id": self.next_id, "method": method, "params": params}) + "\n").encode())
        self.file.flush()
        reply = json.loads(self.file.readline())
        if not reply.get("ok"):
            raise SystemExit(f"{method} {params}: {reply.get('error')}")
        return reply.get("result")

    def execute(self, command, params=None):
        return self.send("engine.execute", {"command": command, "params": params or {}})

    def close(self):
        try:
            self.send("app.quit", {})
        except (OSError, SystemExit):
            pass
        self.file.close()
        self.sock.close()


def solve(ctl):
    """Open the sample, make the step store FRAMES frames, and solve the flow."""
    ctl.execute("file.sample", {"name": SAMPLE})
    # The sample's step stores 10 frames; the clip wants a smoother transient.
    ctl.execute("model.step.edit", {"name": STEP, "procedure": {"type": "fluid", "domain": DOMAIN, "end_time": END_TIME, "frames": FRAMES}})
    log("solving the flow (lattice Boltzmann) …")
    report = ctl.execute("solve.fluid", {"backend": "cpu"})
    forces = report.get("forces") or report.get("faces")
    log(f"solved: {json.dumps(report.get('stats', report.get('summary', {})))[:200]}")
    return forces


def render_frames(ctl, frames_dir):
    """Render the model-only viewport at each field frame, with a fixed legend clamp and the iso
    view, and draw its time caption. Return the list of (picture, caption) pairs."""
    # The flow display the sample opens with: velocity on the mid-slice, the obstacle coloured and
    # streamlines upstream. Fix the legend clamp so the colours are comparable across frames.
    # A light background and shaded part, to match the gallery's still pictures.
    ctl.send("ui.set", {"dark": False, "style": "shaded_edges", "show_tree": False})
    ctl.execute("results.contour", {"field": "velocity", "range": {"min": CLAMP[0], "max": CLAMP[1]}})
    ctl.execute("results.surface", {"show": True})
    ctl.execute("results.streamlines", {"show": True})
    ctl.send("ui.camera", CAMERA)
    font = caption.Font.load(CAPTION_FONT)
    shots = []
    for k in range(FRAMES + 1):
        ctl.execute("results.frame", {"frame": k})
        shot = frames_dir / f"frame{k:04d}.png"
        ctl.send("ui.render", {"path": str(shot), "width": WIDTH, "height": HEIGHT})
        t = END_TIME * k / FRAMES
        bar = frames_dir / f"cap{k:04d}.png"
        bar.write_bytes(caption.caption_bar(font, (WIDTH, CAPTION_HEIGHT), "flow past a bracket", f"t = {t:4.1f} s   ·   Re 83"))
        shots.append((shot, bar))
    return shots


def sequence(shots, folder):
    """Lay the clip out as numbered frames and captions (the last picture held), returning the two
    ffmpeg input patterns and the frame count."""
    shutil.rmtree(folder, ignore_errors=True)
    (folder / "caps").mkdir(parents=True)
    n = 0
    held = [(shot, bar, 1) for shot, bar in shots]
    held[-1] = (held[-1][0], held[-1][1], HOLD_END)
    for shot, bar, count in held:
        for _ in range(count):
            os.link(shot, folder / f"{n:05d}.png")
            os.link(bar, folder / "caps" / f"{n:05d}.png")
            n += 1
    return str(folder / "%05d.png"), str(folder / "caps" / "%05d.png"), n


def encode(shots, frames_dir, mp4, gif=None):
    """Encode the frames into the MP4 (committed to the repo) and, when `gif` is given, a GIF too
    (not committed — it is several times larger; use it for hosting the clip elsewhere)."""
    ffmpeg = os.environ.get("FFMPEG", shutil.which("ffmpeg") or "ffmpeg")
    pattern, caps, count = sequence(shots, frames_dir / "sequence")
    log(f"encoding {count} frames ({count / FPS:.1f} s) into {mp4.name}" + (f" and {gif.name}" if gif else ""))
    src = ["-framerate", str(FPS), "-i", pattern, "-framerate", str(FPS), "-i", caps]
    # Each picture padded with a band below, the caption bar laid over the band.
    captioned = f"[0:v]pad={WIDTH}:{HEIGHT + CAPTION_HEIGHT}:0:0[w];[w][1:v]overlay=0:{HEIGHT}"
    common = [ffmpeg, "-hide_banner", "-loglevel", "error", "-y"]
    mp4.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [*common, *src, "-filter_complex", captioned, "-c:v", "libx264", "-preset", "veryslow", "-crf", "26",
         "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(mp4)],
        check=True,
    )
    if gif:
        subprocess.run(
            [*common, *src, "-filter_complex",
             f"{captioned},scale=800:-1:flags=lanczos,split[a][b];[a]palettegen=stats_mode=diff[p];"
             "[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle", str(gif)],
            check=True,
        )


def still(ctl, path):
    """A clean iso render of the settled flow for the gallery still."""
    ctl.execute("results.frame", {"frame": FRAMES})
    ctl.send("ui.render", {"path": str(path), "width": WIDTH, "height": HEIGHT})
    log(f"wrote the still {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "media" / "examples", help="folder for the clip and still")
    ap.add_argument("--frames", type=Path, help="keep the rendered pictures in this folder")
    ap.add_argument("--gif", action="store_true", help="also write the GIF (not committed; for hosting the clip elsewhere)")
    ap.add_argument("--no-build", action="store_true", help="use the binaries already built")
    args = ap.parse_args()

    if not args.no_build:
        build()
    bin_dir = target_dir() / "release"
    work = Path(tempfile.mkdtemp(prefix="fieldcraft-flow-clip-"))
    frames_dir = args.frames.resolve() if args.frames else work / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    port = free_port()
    app_log = open(work / "app.log", "w")
    env = {**os.environ, "FIELDCRAFT_GPU": os.environ.get("FIELDCRAFT_GPU", "software")}
    app = subprocess.Popen([str(bin_dir / "fieldcraft"), "--control", str(port)], stdout=app_log, stderr=subprocess.STDOUT, env=env)
    ctl = None
    try:
        ctl = Control(port)
        forces = solve(ctl)
        if forces:
            log(f"face forces: {json.dumps(forces)[:400]}")
        shots = render_frames(ctl, frames_dir)
        still(ctl, args.out / "flow-past-flow.png")
        gif = args.out / "flow-past-transient.gif" if args.gif else None
        encode(shots, frames_dir, args.out / "flow-past-transient.mp4", gif)
    finally:
        if ctl:
            ctl.close()
        app.terminate()
        try:
            app.wait(timeout=10)
        except subprocess.TimeoutExpired:
            app.kill()
        app_log.close()
        if not args.frames:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
