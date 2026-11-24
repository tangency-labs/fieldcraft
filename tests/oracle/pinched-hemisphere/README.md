# Hemispherical shell with an 18° hole under alternating point loads

Our model of the published problem: a hemispherical shell of radius 10 and thickness 0.04, with the
polar 18° cap removed, E = 6.825 × 10⁷, ν = 0.3. Four radial point forces of 2.0 act on the equator at
90° intervals, alternately outward and inward. The reference is the radial displacement under the
load, 0.0940 with lower bound 0.0924 from MacNeal and Harder's reference [7] (*A proposed standard set
of problems to test finite element accuracy*, FEAD 1 (1985) 3–20, Fig. 9 (p. 9), Table 5(b) (p. 11)).

The quadrant model: longitude λ from 0° to 90°, latitude φ from the equator (φ = 0) to the hole edge at
φ = 72° (18° from the pole), with R cos φ cos λ = x, R cos φ sin λ = y and R sin φ = z on the
mid-surface R = 10. Both load points lie on its symmetry planes, so each carries half its force: +1.0
along +x at (10, 0, 0) (outward) and −1.0 along y at (0, 10, 0) (inward). The measured nodes are the two
load points, so the case reports the outward displacement u_x and the inward displacement u_y.

The hemisphere is free to move rigidly along z, and the loads have no z component. The displacement
BC on the hole-edge node (0, n) in the plane y = 0 holds u_z = 0 and removes that rigid motion, so
the reaction along z is zero by equilibrium.

## The mesh

`mesh.block`'s 2D form is planar and takes its corners at z = 0, so it cannot mesh a patch on a
doubly curved surface: the recipe writes the 81 nodes and 64 quad4 elements one by one with
`mesh.node.add` and `mesh.element.add`, the route [docs/verification.md](../../../docs/verification.md)
allows ("nodes and elements written one by one"), as `pinched-cylinder` and `scordelis-lo-roof-shell`
do. Node 1 + i + 9 j lies at (10 cos φ_j cos λ_i, 10 cos φ_j sin λ_i, 10 sin φ_j) with
λ_i = (π/2) i / 8 and φ_j = 72° j / 8; element 1 + i + 8 j is
[node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)], counter-clockwise seen from outside
the shell, so its normal points outward.

The recipe was generated with this script (`python3 gen_hemisphere.py 8`, Python 3 stdlib):

```python
import json, math, sys

R = 10.0
PHI_MAX = math.radians(72.0)   # hole edge: 18 degrees from the pole
LAM_MAX = math.pi / 2          # quadrant, longitude 0 to 90 degrees

def hemisphere(n):
    label = lambda i, j: 1 + i + (n + 1) * j
    commands = []
    for j in range(n + 1):
        for i in range(n + 1):
            lam, phi = LAM_MAX * i / n, PHI_MAX * j / n
            xyz = [R * math.cos(phi) * math.cos(lam), R * math.cos(phi) * math.sin(lam), R * math.sin(phi)]
            commands.append({"command": "mesh.node.add", "params": {"label": label(i, j), "xyz": xyz}})
    for j in range(n):
        for i in range(n):
            nodes = [label(i, j), label(i + 1, j), label(i + 1, j + 1), label(i, j + 1)]
            commands.append({"command": "mesh.element.add",
                             "params": {"label": 1 + i + n * j, "kind": "quad4", "nodes": nodes}})
    return commands + [
        {"command": "model.set.create", "params": {"name": "shell", "elements": list(range(1, n * n + 1))}},
        {"command": "model.set.create", "params": {"name": "sym_y", "nodes": [label(0, j) for j in range(n + 1)]}},
        {"command": "model.set.create", "params": {"name": "sym_x", "nodes": [label(n, j) for j in range(n + 1)]}},
        {"command": "model.set.create", "params": {"name": "load_x", "nodes": [label(0, 0)]}},
        {"command": "model.set.create", "params": {"name": "load_y", "nodes": [label(n, 0)]}},
        {"command": "model.set.create", "params": {"name": "hold_z", "nodes": [label(0, n)]}},
        {"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 6.825e7, "poisson": 0.3}}},
        {"command": "model.section.create", "params": {"name": "shell", "set": "shell", "material": "Mat", "kind": {"type": "shell", "thickness": 0.04}}},
        {"command": "model.boundary_condition.create", "params": {"name": "sym_y", "type": "symmetry", "target": {"set": "sym_y"}, "plane": "y"}},
        {"command": "model.boundary_condition.create", "params": {"name": "sym_x", "type": "symmetry", "target": {"set": "sym_x"}, "plane": "x"}},
        {"command": "model.boundary_condition.create", "params": {"name": "hold_z", "type": "displacement", "target": {"set": "hold_z"}, "values": [None, None, 0, None, None, None]}},
        {"command": "model.load.create", "params": {"name": "push", "type": "force", "target": {"set": "load_x"}, "vector": [1.0, 0, 0]}},
        {"command": "model.load.create", "params": {"name": "pull", "type": "force", "target": {"set": "load_y"}, "vector": [0, -1.0, 0]}},
        {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["sym_y", "sym_x", "hold_z"], "loads": ["push", "pull"]}},
        {"command": "solve.static", "params": {}, "id": "solve"},
    ]

n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
print("[\n " + ",\n ".join(json.dumps(c) for c in hemisphere(n)) + "\n]")
```

## Refinement

Every row is the listing above run with `n` set to the row's divisions through `cargo xtask oracle --case`
(the case's recipe regenerated for the run, then restored); the case mesh is 8 × 8. The reference is 0.0940.

| cells | u | error |
|---|---|---|
| 4 × 4 | 0.086855 | −7.6% |
| 8 × 8 (this case) | 0.092253 | −1.86% |
| 16 × 16 | 0.092828 | −1.25% |
| 32 × 32 | 0.093271 | −0.78% |
| 64 × 64 | 0.093491 | −0.54% |

The MITC4 shell approaches the reference from below and keeps closing the gap, so the mesh
error falls monotonically. The 4 × 4 mesh is too coarse for the 2% grade; the case mesh is the
prescribed 8 × 8 of MacNeal and Harder.

## Tolerance

2%: MacNeal and Harder grade an element A within 2% of the theoretical value on this prescribed coarse
mesh (p. 14), and their reference value is 0.0940 with the lower bound 0.0924 from their ref. [7]. The
8 × 8 mesh is 1.86% low; its own refinement moves the value monotonically towards the reference, so the
tolerance is the source's grading on the mesh it prescribes (docs/verification.md); the 8 × 8 result uses
1.86 of its 2 points. The second quantity, v, is the inward load's displacement along −y: the same
magnitude with the opposite sign, a check of the model's symmetry rather than a second reference.

`fixture.json` is the accepted value at 8 × 8 with the runner's default band, once blessed: the mesh
is written by the recipe, not generated from geometry, so it is bit-identical on every platform.

## Running it

```sh
cargo xtask oracle --case pinched-hemisphere
cargo xtask oracle --bless pinched-hemisphere   # only after the case passes
```
