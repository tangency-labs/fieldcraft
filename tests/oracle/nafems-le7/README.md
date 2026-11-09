# NAFEMS LE7 axisymmetric cylinder/sphere under pressure (axisymmetric quad8)

Our model of the published problem: a thin-walled vessel, a cylinder closed by a hemispherical
cap, of mid-surface radius 1 m and wall thickness 0.025 m, E = 210 GPa, ν = 0.3 (m, MN, MPa),
under an internal pressure of 1 MPa. In FieldCraft's axisymmetric models x is the radius and y the
axis. The cylinder runs from F (y = 0) to the junction C (y = 1.5), the cap from C to the pole A
(y = 2.5). The axial displacement is zero at F (the vessel's plane of symmetry: the whole bottom
edge) and the radial one at A (the pole edge, on the axis). The target is the axial stress σ_yy
on the outer surface at D (R = 1.0125, y = 1.4034), 0.0966 m below the junction, inside the
bending boundary layer the cap's different stiffness causes.

## The mesh

Four `mesh.block` quad8 patches through the wall thickness: three cylinder blocks with E
(y = 1.1136) and D as block corners, so D is a node, and the cap, a quarter ring whose two curved
edges are arcs about (0, 1.5) and whose fourth edge lies on the axis. 2 elements through the
wall, about 50 per metre along it (56, 14 and 5 in the cylinder blocks, 79 round the cap); the
blocks share their edges, so `mesh.block` joins them. The pressure acts on the inner edges of all
four blocks.

The recipe was generated with this script (`python3 le7.py 2 50 > recipe.json`):

```python
import json, math, sys

def le7(nt, per_m):
    """NAFEMS LE7 cylinder with a hemispherical cap, axisymmetric (m, MN, MPa): x is the radius,
    y the axis. Mid-surface radius 1, wall 0.025 thick; the cylinder runs from F (y = 0, the plane
    of symmetry) to C (y = 1.5), the cap from C to the pole A (y = 2.5). nt elements through the
    wall, about per_m elements per metre along it."""
    ri, ro = 1 - 0.0125, 1 + 0.0125
    cuts = [0.0, 1.1136, 1.4034, 1.5]  # F, E, D, C: E and D are nodes
    along = lambda length: max(1, round(length * per_m))
    commands = []
    for k, (y0, y1) in enumerate(zip(cuts, cuts[1:])):
        commands.append({"command": "mesh.block", "params": {"body": f"cylinder_{k + 1}", "kind": "quad8",
            "corners": [[ri, y0, 0], [ro, y0, 0], [ro, y1, 0], [ri, y1, 0]], "divisions": [nt, along(y1 - y0)],
            "edges": [{"corners": [1, 4], "name": "inside"}, {"corners": [1, 2], "name": "bottom"}]}})
    commands.append({"command": "mesh.block", "params": {"body": "cap", "kind": "quad8",
        "corners": [[ri, 1.5, 0], [ro, 1.5, 0], [0, 1.5 + ro, 0], [0, 1.5 + ri, 0]], "divisions": [nt, along(math.pi / 2)],
        "edges": [{"corners": [2, 3], "arc": {"center": [0, 1.5, 0]}}, {"corners": [4, 1], "arc": {"center": [0, 1.5, 0]}, "name": "inside"},
                  {"corners": [3, 4], "name": "pole"}]}})
    bodies = ["cylinder_1", "cylinder_2", "cylinder_3", "cap"]
    return commands + [
        {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": b, "face": "area"} for b in bodies]}},
        {"command": "model.set.create", "params": {"name": "inside", "geometry": [{"body": b, "edge": "inside"} for b in bodies]}},
        {"command": "model.set.create", "params": {"name": "F", "geometry": [{"body": "cylinder_1", "edge": "bottom"}]}},
        {"command": "model.set.create", "params": {"name": "A", "geometry": [{"body": "cap", "edge": "pole"}]}},
        {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}}},
        {"command": "model.section.create", "params": {"name": "vessel", "set": "all", "material": "Steel", "kind": {"type": "axisymmetric"}}},
        {"command": "model.boundary_condition.create", "params": {"name": "F", "type": "symmetry", "target": {"set": "F"}, "plane": "y"}},
        {"command": "model.boundary_condition.create", "params": {"name": "A", "type": "symmetry", "target": {"set": "A"}, "plane": "x"}},
        {"command": "model.load.create", "params": {"name": "pressure", "type": "pressure", "set": "inside", "magnitude": 1.0}},
        {"command": "model.step.create", "params": {"name": "Pressure", "procedure": {"type": "static"}, "boundary_conditions": ["F", "A"], "loads": ["pressure"]}},
        {"command": "solve.static", "params": {}, "id": "solve"},
    ]

nt, per_m = (int(a) for a in sys.argv[1:3]) if len(sys.argv) > 2 else (2, 50)
print("[\n " + ",\n ".join(json.dumps(c) for c in le7(nt, per_m)) + "\n]")
```

## Refinement

Every row is the listing above with the row's divisions, solved and probed at D as the case does;
the case mesh is 2 × 50:

| through the wall × per metre | σ_yy at D | against 25.86 | against 25.27 |
|---|---|---|---|
| 2 × 25 | 25.331 | −2.05% | +0.24% |
| 2 × 50 (this case) | 25.283 | −2.23% | +0.05% |
| 4 × 100 | 25.271 | −2.28% | +0.00% |
| 8 × 200 | 25.269 | −2.29% | −0.00% |

Far from the junction the model is exact: at F the axial stress on the outer surface is 19.503
(p r_i²/(r_o² − r_i²)) and the hoop stress 39.006 (Lamé).

## Two references, one known failure

The published target, 25.86, comes from axisymmetric thin-shell theory. The ESRD guide solves the
same problem with axisymmetric solid elements and reports 25.27 on two meshes, 2.3% below the
target, and our solid model converges to the same 25.27: the difference is the model, not the
mesh. So the case keeps the NAFEMS target with the
tolerance a shell model would earn, 1%, marked `known_failure` until FieldCraft has an
axisymmetric shell element, and checks the solid model against the solid reference, also at 1%.

`fixture.json` holds the passing quantity's 2 × 50 value: the mesh is built by `mesh.block` from
the recipe's corners, so it is the same on every platform.

## Running it

```sh
cargo xtask oracle --case nafems-le7
cargo xtask oracle --bless nafems-le7   # only after the case passes
```


Historical boundary: the frozen recipe omits `model.units` and is therefore read
in the default mm–tonne–s–N–MPa system. Its geometry is 1e-3 of the stated metre
geometry, while E and pressure convert to 210 GPa and 1 MPa. Uniform scaling
preserves this pressure-only linear stress comparison, not displacements,
forces, mass or frequencies. It remains an axisymmetric **solid** recipe and
its shell-reference known failure is retained. The separate axisymmetric-shell
pole/LE7 command test uses explicit SI units and compares outer shell stresses;
it does not overwrite this model or fixture.
