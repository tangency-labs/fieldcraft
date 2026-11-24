# NAFEMS FV2 pin-ended double cross, in-plane vibration (planar beams in a 2D model)

Our model of the published problem: eight arms of length 5 m at 0°, 45°, … 315° from a common
joint, each pinned at its far end (x and y held, rotation free); square section 0.125 m; E = 200 GPa,
ν = 0.3, ρ = 8000 kg/m³; SI numbers, no `model.units` command. Each arm is four `beam2` elements.

A plane stress `quad4` (thickness 0.01 m) spans the far ends of the arms at 0°, 90°, 180° and 270°,
which makes the model 2D. All of its nodes are pinned, so it carries no load and adds no stiffness
to the arms; it is there because a model is 2D by its 2D sections (a frame alone is a 3D model,
`asks/b1069-1.md`). The planar beams share the plate's nodes, which is the case the planar-beam
solver is built for.

Targets (NAFEMS TNSB Rev. 3, 1990, test FV2): 11.336 Hz (mode 1), 17.709 Hz (modes 2–8),
45.345 Hz (mode 9), 57.390 Hz (modes 10–16). The source is F. Abbassian, D. J. Dawswell and
N. C. Knowles, *Selected Benchmarks for Natural Frequency Analysis*, NAFEMS 1987, test 2.

## Result

Sixteen frequencies, 1 % relative tolerance, from `solve.frequency` on the `Modes` step.

| Mode | Ours (Hz) | Target (Hz) | Error |
|---|---|---|---|
| 1 | 11.328 | 11.336 | −0.07 % |
| 2, 3 | 17.650 | 17.709 | −0.33 % |
| 4–8 | 17.679 | 17.709 | −0.17 % |
| 9 | 45.361 | 45.345 | +0.04 % |
| 10, 11 | 57.094 | 57.390 | −0.51 % |
| 12–16 | 57.410 | 57.390 | +0.03 % |

## Why the targets differ by a few tenths of a percent

The targets are Euler–Bernoulli values for the pinned–pinned and clamped–pinned arms:
λ²/(2πL²)·√(EI/(ρA)) with λ = π, 3.9266, 2π and 7.0686. The Timoshenko element adds shear
deformation and rotary inertia, which lowers a pinned–pinned arm by 0.10 % on mode 1 and by 0.41 %
on mode 9 in closed form. Raising the arms' cross-section area by 10⁴ at the same second moment of
area and the same mass per length (which removes shear flexibility, axial flexibility and rotary
inertia together) brings all 16 modes to within 0.002 % of the targets at 16 elements per arm, so
the differences are those of the Timoshenko model and of the mesh, not of the model's joint or
boundary conditions. At 16 elements per arm the unscaled arms are −0.10 % (mode 1), −0.41 % (mode 9)
and −1.20 % (modes 10 and 11); at 4 elements the cubic mesh is up to +0.62 % high on modes 10 to 16,
which is why 4 elements show only −0.51 % on modes 10 and 11 (the two effects partly cancel). The
area scaling cannot tell shear from axial flexibility apart.

## Generator

The recipe is produced by this script (`python3 -I gen_fv2.py 4 > recipe.json`):

```python
"""NAFEMS FV2, the pin-ended double cross, as FieldCraft commands (SI numbers).

python3 -I gen_fv2.py [elements per arm, default 4] > recipe.json
"""
import json
import math
import sys

n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
ARM, SIDE = 5.0, 0.125
cmds = [{"command": "mesh.node.add", "params": {"label": 1, "xyz": [0.0, 0.0, 0.0]}}]
beams, ends, element = [], [], 0
for k in range(8):
    angle = k * math.pi / 4.0
    previous = 1
    for i in range(1, n + 1):
        label = 1 + k * n + i
        s = ARM * i / n
        cmds.append({"command": "mesh.node.add", "params": {"label": label, "xyz": [s * math.cos(angle), s * math.sin(angle), 0.0]}})
        element += 1
        cmds.append({"command": "mesh.element.add", "params": {"label": element, "kind": "beam2", "nodes": [previous, label]}})
        beams.append(element)
        previous = label
    ends.append(previous)
plate = element + 1
cmds.append({"command": "mesh.element.add", "params": {"label": plate, "kind": "quad4", "nodes": [ends[0], ends[2], ends[4], ends[6]]}})
cmds += [
    {"command": "model.set.create", "params": {"name": "arms", "elements": beams}},
    {"command": "model.set.create", "params": {"name": "plate", "elements": [plate]}},
    {"command": "model.set.create", "params": {"name": "pins", "nodes": ends}},
    {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 200e9, "poisson": 0.3}, "density": 8000.0}},
    {"command": "model.section.create", "params": {"name": "arms", "set": "arms", "material": "Steel", "kind": {"type": "beam", "profile": {"type": "rectangle", "width": SIDE, "height": SIDE}, "orientation": [0, 0, 1]}}},
    {"command": "model.section.create", "params": {"name": "plate", "set": "plate", "material": "Steel", "kind": {"type": "plane_stress", "thickness": 0.01}}},
    {"command": "model.boundary_condition.create", "params": {"name": "pins", "type": "displacement", "target": {"set": "pins"}, "values": [0, 0, None, None, None, None]}},
    {"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 16}, "boundary_conditions": ["pins"], "loads": []}},
    {"command": "solve.frequency", "params": {"step": "Modes"}, "id": "modes"},
]
print("[\n " + ",\n ".join(json.dumps(c) for c in cmds) + "\n]")
```

## Running it

```sh
cargo xtask oracle --case nafems-fv2
```
