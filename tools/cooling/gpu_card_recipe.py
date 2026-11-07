#!/usr/bin/env python3
"""Generate the graphics-card SolidCraft recipe (crates/engine/samples/cooling-gpu-card.json).

A generic two-slot, blower-cooled high-end graphics card, about 267 x 111 x 40 mm (our own layout:
no vendor's design, logo or branding). The card lies with its board in the x-y plane, the cooler
on top (+z) and the backplate below:

- the board (FR-4, 1.6 mm) and an aluminium backplate under it;
- the GPU package and its die (26 x 26 mm, 6.8 cm^2) at the board's middle, eight memory chips
  around it and a column of six VRM blocks towards the bracket;
- a copper base (vapour chamber) on the die, three heat pipes running from it towards the bracket,
  and two aluminium fin stacks: one on the base over the GPU, one on the heat pipes over the VRMs;
- the shroud: a 1.5 mm shell over the cooler with a round intake over the blower;
- the radial blower at the far end: a hub and backward-curved blades (visual: the solver drives the
  air with a fan zone, it has no rotating geometry);
- the I/O bracket at x = 0 with its exhaust slots.

Air enters through the intake, the blower pushes it along the fins towards the bracket (-x), and it
leaves through the slots. Every length is in millimetres.

Python standard library only. Run from the repository root:
    python3 tools/cooling/gpu_card_recipe.py   # writes crates/engine/samples/cooling-gpu-card.json
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "crates" / "engine" / "samples" / "cooling-gpu-card.json"

LENGTH, HEIGHT = 267.0, 111.0      # the board, along x and y
BOARD = 1.6                        # board thickness
TOP = 38.0                         # the shroud's outer top, z (two slots with the backplate)
SHELL = 1.5                        # shroud wall
GPU = (150.0, 55.5)                # package centre on the board
PACKAGE, DIE = 45.0, 26.0          # sides
FAN = (228.0, 55.5)                # blower axis
FAN_R, HUB_R = 34.0, 13.0          # impeller tip and hub radii
INTAKE_R = 30.0                    # the round opening in the shroud over the blower
BLADES = 28
FIN_PITCH, FIN_T = 3.0, 1.0        # fin pitch and thickness (y)


def r(v):
    return round(v, 3)


def box(corner, size, name=None, op=None, targets=None):
    p = {"corner": [r(c) for c in corner], "length": r(size[0]), "width": r(size[1]), "height": r(size[2])}
    if name:
        p["body_name"] = name
    if op:
        p["operation"] = op
        p["targets"] = targets
    return {"command": "solid.box", "params": p}


def centred(cx, cy, z0, sx, sy, h, name=None, op=None, targets=None):
    return box([cx - sx / 2, cy - sy / 2, z0], [sx, sy, h], name, op, targets)


def join(cmd, target):
    """`cmd` (a solid.box or solid.cylinder) joined into body `target`."""
    cmd["params"].pop("body_name", None)
    cmd["params"]["operation"] = "join"
    cmd["params"]["targets"] = [target]
    return cmd


def cylinder(base, axis, radius, height, name):
    return {"command": "solid.cylinder", "params": {"base": [r(v) for v in base], "axis": axis, "radius": r(radius), "height": r(height), "body_name": name}}


def recipe():
    """The bodies in the order that gives them their colours. The viewer's palette cycles
    blue-grey, orange, green, pink, lilac, yellow, teal, tan by the order of the bodies drawn;
    the pictures hide the shroud (body 0) to show the inside, so every other body takes the
    colour one place before its index: the board green, the base and heat pipes copper-orange."""
    c = []
    gx, gy = GPU
    pz = BOARD
    die_top = pz + 2.2
    base_top = die_top + 4.0
    fin_top = TOP - SHELL - 0.5
    # 0 Shroud (blue-grey): a shell with filleted top edges, open at the bracket end and the
    # board, a recessed ring round the round intake and three shallow grooves along the top.
    # (from the bracket's inner face, x 0, as a shroud is fixed to its bracket: a gap between them
    # would be a slit open to the intake's face of the air box)
    c.append(box([0.0, 0.5, BOARD], [LENGTH - 0.5, HEIGHT - 1.0, TOP - BOARD], "Shroud"))
    c.append({"command": "solid.fillet", "params": {"edges": [[133.5, 0.5, TOP], [133.5, HEIGHT - 0.5, TOP], [LENGTH - 0.5, 55.5, TOP]], "radius": 4}})
    c.append(box([2.0, 2.0, BOARD - 1.0], [LENGTH - 4.0, HEIGHT - 4.0, TOP - BOARD - SHELL + 1.0], op="cut", targets=["Shroud"]))
    # (inset half a millimetre from the inner cut's side walls: coplanar faces would leave the
    # shell's surface open)
    c.append(box([-1.0, 2.5, BOARD + 2.0], [4.0, HEIGHT - 5.0, TOP - BOARD - SHELL - 4.0], op="cut", targets=["Shroud"]))
    c.append({"command": "solid.cylinder", "params": {"base": [FAN[0], FAN[1], TOP - 0.6], "axis": "Z", "radius": INTAKE_R + 5.0,
                                                     "height": 1.0, "operation": "cut", "targets": ["Shroud"]}})
    c.append({"command": "solid.cylinder", "params": {"base": [FAN[0], FAN[1], TOP - SHELL - 1.0], "axis": "Z", "radius": INTAKE_R,
                                                     "height": SHELL + 2.0, "operation": "cut", "targets": ["Shroud"]}})
    for gy_ in (18.0, 55.5, 93.0):
        c.append(centred(95.0, gy_, TOP - 0.6, 170.0, 2.5, 1.0, op="cut", targets=["Shroud"]))
    # 1 Backplate (blue-grey), in contact with the board (thermal pads).
    c.append(box([4, 2, -2.0], [LENGTH - 6, HEIGHT - 4, 2.0], "Backplate"))
    # 2 Base (orange, copper): the vapour chamber on the die.
    c.append(centred(gx, gy, die_top, 72.0, 62.0, 4.0, "Base"))
    # 3 Board (green).
    c.append(box([0, 0, 0], [LENGTH, HEIGHT, BOARD], "Board"))
    # 4 VRM (pink): six blocks in a column towards the bracket, one body.
    for k in range(6):
        cmd = centred(78.0, 14.0 + k * 16.6, BOARD, 9.0, 11.0, 5.0, "VRM")
        c.append(cmd if k == 0 else join(cmd, "VRM"))
    # 5 Package (lilac).
    c.append(centred(gx, gy, pz, PACKAGE, PACKAGE, 1.4, "Package"))
    # 6 Memory (yellow): eight chips round the package, one body.
    memory = [(gx - 38, gy - 33), (gx - 38, gy), (gx - 38, gy + 33), (gx + 38, gy - 33), (gx + 38, gy), (gx + 38, gy + 33),
              (gx, gy - 40), (gx, gy + 40)]
    for k, (mx, my) in enumerate(memory):
        cmd = centred(mx, my, BOARD, 14.0, 12.0, 1.2, "Memory")
        c.append(cmd if k == 0 else join(cmd, "Memory"))
    # 7 Bracket: the I/O bracket with its exhaust slots.
    c.append(box([-2.0, -6.0, -4.0], [2.0, HEIGHT + 12.0, TOP + 4.0], "Bracket"))
    for k in range(7):
        c.append(box([-3.0, 8.0 + k * 13.6, 6.0], [4.0, 9.0, TOP - 10.0], op="cut", targets=["Bracket"]))
    # 8 Die (tan).
    c.append(centred(gx, gy, pz + 1.4, DIE, DIE, 0.8, "Die"))
    # 9 Fins (teal, aluminium): one stack on the base over the GPU, one on the heat pipes over the VRMs.
    first = True
    for x0, x1, y0, y1, z0 in ((gx - 36.0, gx + 36.0, 22.0, 89.0, base_top), (42.0, 104.0, 8.0, 103.0, 8.5)):
        root = 1.2
        cmd = box([x0, y0, z0], [x1 - x0, y1 - y0, root], "Fins")
        c.append(cmd if first else join(cmd, "Fins"))
        first = False
        y = y0
        while y + FIN_T <= y1 + 1e-9:
            c.append(box([x0, y, z0 + root], [x1 - x0, FIN_T, fin_top - z0 - root], op="join", targets=["Fins"]))
            y += FIN_PITCH
    # 10 HeatPipes (orange, copper): three pipes from the base to the bracket end, one body.
    pipe_r = 3.0
    for k, py in enumerate((gy - 15.0, gy, gy + 15.0)):
        cmd = cylinder([gx + 20.0, py, base_top + pipe_r], [-1, 0, 0], pipe_r, gx + 20.0 - 40.0, "HeatPipes")
        c.append(cmd if k == 0 else join(cmd, "HeatPipes"))
    # 11 Impeller (visual): backward-swept blades, one sketched and patterned round the axis,
    # joined to the hub.
    fx, fy = FAN
    z0, h = BOARD + 1.0, TOP - SHELL - BOARD - 3.0
    sweep = math.radians(28.0)  # the tip trails the root by this angle

    def at(radius, angle, side):
        a = angle + side * 0.6 / radius  # half a 1.2 mm blade thickness, as an angle at that radius
        return [r(fx + radius * math.cos(a)), r(fy + radius * math.sin(a))]

    pts = [at(HUB_R - 1.0, 0.0, -1), at(FAN_R, -sweep, -1), at(FAN_R, -sweep, 1), at(HUB_R - 1.0, 0.0, 1)]
    c.append({"command": "solid.cylinder", "params": {"base": [fx, fy, r(z0)], "axis": "Z", "radius": HUB_R, "height": r(h), "body_name": "Impeller"}})
    c.append({"command": "sketch.create", "params": {"plane": "XY", "offset": r(z0), "name": "BladeSketch"}})
    ids = ["b1", "b2", "b3", "b4"]
    for k in range(4):
        c.append({"command": "sketch.line", "params": {"points": [pts[k], pts[(k + 1) % 4]], "ids": [ids[k]]}})
    c.append({"command": "sketch.finish", "params": {}})
    c.append({"command": "solid.extrude", "params": {"sketch": "BladeSketch", "profiles": [{"loops": [{"curves": ids, "outer": True}], "sketch": "BladeSketch"}],
                                                    "distance": r(h), "operation": "join", "targets": ["Impeller"], "name": "Blade"}})
    c.append({"command": "solid.pattern.circular", "params": {"features": ["Blade"], "axis": {"origin": [fx, fy, 0], "dir": [0, 0, 1]},
                                                             "count": BLADES, "angle": 360, "name": "Blades"}})
    return c


def recipe_json(commands):
    """The recipe, one command per line (the form of crates/engine/samples/cooling-aio.json)."""
    lines = ",\n".join("  " + json.dumps(cmd) for cmd in commands)
    return '{"commands": [\n' + lines + "\n]}\n"


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cmds = recipe()
    OUT.write_text(recipe_json(cmds))
    print(f"gpu_card_recipe: wrote {OUT.relative_to(ROOT)} ({len(cmds)} commands)", file=sys.stderr)
