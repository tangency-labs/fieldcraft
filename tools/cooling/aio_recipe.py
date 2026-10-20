#!/usr/bin/env python3
"""Generate the AIO cooler SolidCraft recipes (crates/engine/samples/cooling-aio.json and cooling-aio-half.json).

A 240-class all-in-one liquid cooler: a cold plate and a two-tier pump block, a perforated
radiator (a filleted core with a recessed fin pocket behind a raised border, and horizontal air
channels cut all the way through its 30 mm thickness, so a fan blows straight through), and two
coolant tubes running from barbs on the back of the pump to ports on the radiator's bottom border.

The tubes are a close parallel pair: one smooth centreline (straight 10 mm out of the pump, a cubic
Bezier S with vertical tangents at both ends, straight 10 mm into the radiator), offset along its
normal by half the centre spacing on either side. A true normal offset keeps the gap constant along
the whole run (a sideways translate does not: where the curve runs along x the copies close up).
The script checks the clearance and the bend radius and refuses to write a recipe that would make
the tubes touch or fold.

A second recipe, the half (the cutaway), is the same cooler with the radiator cut at its middle plane (x = 15,
keeping the half with the ports): the airflow study simulates that half with a symmetry plane on the
cut, which is also what lets the pictures look into the channels.

Python standard library only. Run from the repository root:
    python3 tools/cooling/aio_recipe.py   # writes crates/engine/samples/cooling-aio{,-half}.json
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "crates" / "engine" / "samples" / "cooling-aio.json"
OUT_CUTAWAY = ROOT / "crates" / "engine" / "samples" / "cooling-aio-half.json"
CUT_X = 15.0   # the radiator's middle plane: the frame and the fin pocket are symmetric about it

TUBE = 6.5     # tube outer diameter, mm
BARB = 3.4     # barb and port radius (the tube end sits 2 mm inside it)
HALF = 4.6     # half the tube centre spacing: 9.2 mm centre to centre, a 2.7 mm gap
TZ = 13.0      # tube plane height: barbs and ports sit on the radiator's solid bottom border
XC, XE = 29.0, -77.0   # pair centre at the pump / at the radiator
Y0, Y1 = 8.0, 74.0     # ends of the curved part (after/before the 10 mm straight runs)
HA, HB = 52.0, 48.0    # Bezier handle lengths (vertical tangents at both ends)


def bezier():
    p = [(XC, Y0), (XC, Y0 + HA), (XE, Y1 - HB), (XE, Y1)]

    def at(t):
        u = 1 - t
        return tuple(u**3 * p[0][i] + 3 * u * u * t * p[1][i] + 3 * u * t * t * p[2][i] + t**3 * p[3][i] for i in (0, 1))

    def d1(t):
        u = 1 - t
        return tuple(3 * u * u * (p[1][i] - p[0][i]) + 6 * u * t * (p[2][i] - p[1][i]) + 3 * t * t * (p[3][i] - p[2][i]) for i in (0, 1))

    def d2(t):
        u = 1 - t
        return tuple(6 * u * (p[2][i] - 2 * p[1][i] + p[0][i]) + 6 * t * (p[3][i] - 2 * p[2][i] + p[1][i]) for i in (0, 1))

    return at, d1, d2


def tube_paths(samples=400, every=25):
    """The two offset centrelines as spline fit points, after checking clearance and bend radius."""
    at, d1, d2 = bezier()
    ts = [i / samples for i in range(samples + 1)]

    def offset(t, side):
        x, y = at(t)
        dx, dy = d1(t)
        n = math.hypot(dx, dy)
        return (x + side * HALF * dy / n, y - side * HALF * dx / n)

    def radius(t):
        (dx, dy), (ex, ey) = d1(t), d2(t)
        cross = dx * ey - dy * ex
        return math.inf if cross == 0 else abs((dx * dx + dy * dy) ** 1.5 / cross)

    right = [offset(t, +1) for t in ts]
    left = [offset(t, -1) for t in ts]
    min_radius = min(radius(t) for t in ts)
    min_gap = min(math.hypot(a[0] - b[0], a[1] - b[1]) for a in right[::2] for b in left[::2])
    if min_radius <= HALF + TUBE / 2 + 5 or min_gap <= TUBE + 1.5:
        sys.exit(f"aio_recipe: tubes would touch or fold (bend radius {min_radius:.1f} mm, centre gap {min_gap:.2f} mm)")
    print(f"aio_recipe: tube centre spacing {min_gap:.2f} mm (gap {min_gap - TUBE:.2f} mm), tightest bend {min_radius:.1f} mm", file=sys.stderr)
    fit = lambda pts: [[round(x, 2), round(y, 2)] for x, y in pts[::every]]
    return fit(right), fit(left)


def recipe(cutaway=False):
    c = []
    # Cold plate and pump
    c.append({"command": "solid.box", "params": {"corner": [0, -64, 0], "length": 56, "width": 56, "height": 5, "body_name": "ColdPlate"}})
    c.append({"command": "solid.fillet", "params": {"edges": [[0, -64, 2.5], [56, -64, 2.5], [0, -8, 2.5], [56, -8, 2.5]], "radius": 3}})
    c.append({"command": "solid.cylinder", "params": {"base": [28, -36, 5], "axis": "Z", "radius": 27, "height": 20, "body_name": "Pump"}})
    c.append({"command": "solid.cylinder", "params": {"base": [28, -36, 25], "axis": "Z", "radius": 23, "height": 6, "operation": "join", "targets": ["Pump"]}})
    xr0, xl0, xr1, xl1 = XC + HALF, XC - HALF, XE + HALF, XE - HALF
    for bx in (xr0, xl0):  # barbs from inside the pump out to a common tip at y = -2
        c.append({"command": "solid.cylinder", "params": {"base": [round(bx, 2), -14, TZ], "axis": [0, 1, 0], "radius": BARB, "height": 12, "operation": "join", "targets": ["Pump"]}})
    # Radiator: filleted core, fin pocket recessed behind a border, air channels through the thickness
    c.append({"command": "solid.box", "params": {"corner": [-100, 92, 5], "length": 230, "width": 30, "height": 118, "body_name": "Radiator"}})
    c.append({"command": "solid.fillet", "params": {"edges": [[-100, 92, 64], [130, 92, 64], [-100, 122, 64], [130, 122, 64],
                                                              [15, 92, 123], [15, 122, 123], [15, 92, 5], [15, 122, 5]], "radius": 5}})
    c.append({"command": "solid.box", "params": {"corner": [-86, 90, 16], "length": 202, "width": 4, "height": 96, "operation": "cut", "targets": ["Radiator"]}})
    c.append({"command": "solid.box", "params": {"corner": [-86, 118, 16], "length": 202, "width": 4, "height": 96, "operation": "cut", "targets": ["Radiator"]}})
    z = 19.0
    while z <= 109:
        c.append({"command": "solid.box", "params": {"corner": [-84, 89, round(z, 2)], "length": 198, "width": 34, "height": 2.4, "operation": "cut", "targets": ["Radiator"]}})
        z += 4.6
    for px in (xr1, xl1):  # ports on the bottom border, pointing at the pump
        c.append({"command": "solid.cylinder", "params": {"base": [round(px, 2), 93, TZ], "axis": [0, -1, 0], "radius": BARB, "height": 9, "operation": "join", "targets": ["Radiator"]}})
    # Tubes: straight 10 mm (ends 2 mm inside the barb/port), spline, straight 10 mm
    right, left = tube_paths()
    for x0, x1, pts, name in ((xr0, xr1, right, "Tube1"), (xl0, xl1, left, "Tube2")):
        sketch = f"{name}Path"
        c += [
            {"command": "sketch.create", "params": {"plane": {"origin": [0, 0, TZ], "x_dir": [1, 0, 0], "y_dir": [0, 1, 0]}, "name": sketch}},
            {"command": "sketch.line", "params": {"points": [[round(x0, 2), -4], [round(x0, 2), Y0]], "ids": ["l1"]}},
            {"command": "sketch.spline.fit_point", "params": {"points": pts}},
            {"command": "sketch.line", "params": {"points": [[round(x1, 2), Y1], [round(x1, 2), 86]], "ids": ["l2"]}},
            {"command": "sketch.finish", "params": {}},
            {"command": "solid.pipe", "params": {"path_sketch": sketch, "path": ["l1", "s1", "l2"], "diameter": TUBE, "body_name": name}},
        ]
    if cutaway:  # keep the half with the ports; the cut face is the symmetry plane of the study
        c.append({"command": "solid.box", "params": {"corner": [CUT_X, 80, 0], "length": 130, "width": 60, "height": 130, "operation": "cut", "targets": ["Radiator"]}})
    return c


def recipe_json(commands):
    """The recipe, one command per line (the form of crates/engine/samples/crank-arm.json)."""
    lines = ",\n".join("  " + json.dumps(cmd) for cmd in commands)
    return '{"commands": [\n' + lines + "\n]}\n"


if __name__ == "__main__":
    for path, cut in ((OUT, False), (OUT_CUTAWAY, True)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(recipe_json(recipe(cut)))
        print(f"aio_recipe: wrote {path.relative_to(ROOT)}", file=sys.stderr)
