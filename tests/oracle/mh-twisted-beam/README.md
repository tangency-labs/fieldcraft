# MacNeal–Harder twisted beam

Our model of the published problem: a straight strip along x, length 12, width 1.1, thickness (shell depth)
0.32, twisted 90° about the x axis from the root (x = 0) to the tip (x = 12), E = 29 × 10⁶, ν = 0.22. The
root is clamped; the tip carries a unit total force, once in the width direction (in-plane) and once in the
normal direction (out-of-plane). The measured quantities are the tip-centre displacements: u_z for the
in-plane load and u_y for the out-of-plane load. The theoretical values are 0.005424 and 0.001754
(MacNeal and Harder, *A proposed standard set of problems to test finite element accuracy*, FEAD 1 (1985)
3–20, Fig. 6 (p. 8), Table 3 (p. 10)).

The mesh is the prescribed 12 × 2 of quad4 (MITC4) shell elements. `mesh.block`'s 2D form is planar and
takes z = 0, so the recipe writes the 39 nodes and 24 quad4 elements one by one with `mesh.node.add` and
`mesh.element.add`, the route [docs/verification.md](../../../docs/verification.md) allows ("nodes and
elements written one by one"), as the sibling shell cases do.

## The mesh

With nx × ny cells, node 1 + i + (nx + 1) j, i = 0..nx, j = 0..ny, lies at

x = 12 i / nx, s = −0.55 + 1.1 j / ny, θ = (π/2) · x / 12, xyz = [x, s cos θ, s sin θ].

At the root the width runs along y and the normal along z; at the tip the width runs along z and the normal
along y. Element 1 + i + nx j is [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)]. The sets are
`shell` (all elements), `root` (i = 0), `tip_ends` (nodes (nx, 0) and (nx, ny)) and `tip_inner`
(nodes (nx, j), 0 < j < ny). The measured node is (nx, ny / 2), the tip centre at (12, 0, 0).

## The load

The tip force is a unit total, spread as the consistent nodal loads of a uniform shear across the tip edge:
each end node takes 1 / (2 ny) and each inner node 1 / ny. A node set with the default `total`
distribution splits its load equally over its nodes, so the set totals are `tip_ends` 1 / ny and
`tip_inner` (ny − 1) / ny (both 0.5 for ny = 2). The uniform shear is the load a unit tip traction
produces; giving the end nodes the same weight as the inner nodes would load the edges of the tip more than
its middle. The in-plane step loads the sets along +z (the tip's width), the out-of-plane step along +y (the
tip's normal). Both steps hold the root clamped.

The recipe was generated with this script (`python3 gen_twisted.py 12 2 > recipe.json`):

```python
import json, math, sys

L = 12.0
WIDTH = 1.1
THICKNESS = 0.32
YOUNG = 29.0e6
POISSON = 0.22


def twisted_beam(nx, ny):
    node = lambda i, j: 1 + i + (nx + 1) * j
    commands = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = L * i / nx
            s = -WIDTH / 2 + WIDTH * j / ny
            theta = math.pi / 2 * x / L
            xyz = [x, s * math.cos(theta), s * math.sin(theta)]
            commands.append({"command": "mesh.node.add",
                             "params": {"label": node(i, j), "xyz": xyz}})
    for j in range(ny):
        for i in range(nx):
            nodes = [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)]
            commands.append({"command": "mesh.element.add",
                             "params": {"label": 1 + i + nx * j, "kind": "quad4", "nodes": nodes}})
    tip_ends = [node(nx, 0), node(nx, ny)]
    tip_inner = [node(nx, j) for j in range(1, ny)]
    # A total of 1 per load case, split over the tip: each end node gets 1/(2 ny), each inner node 1/ny.
    end_total = 1.0 / ny
    inner_total = (ny - 1) / ny
    return commands + [
        {"command": "model.set.create", "params": {"name": "shell", "elements": list(range(1, nx * ny + 1))}},
        {"command": "model.set.create", "params": {"name": "root", "nodes": [node(0, j) for j in range(ny + 1)]}},
        {"command": "model.set.create", "params": {"name": "tip_ends", "nodes": tip_ends}},
        {"command": "model.set.create", "params": {"name": "tip_inner", "nodes": tip_inner}},
        {"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": YOUNG, "poisson": POISSON}}},
        {"command": "model.section.create", "params": {"name": "shell", "set": "shell", "material": "Mat", "kind": {"type": "shell", "thickness": THICKNESS}}},
        {"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}},
        {"command": "model.load.create", "params": {"name": "in_plane_ends", "type": "force", "target": {"set": "tip_ends"}, "vector": [0, 0, end_total]}},
        {"command": "model.load.create", "params": {"name": "in_plane_inner", "type": "force", "target": {"set": "tip_inner"}, "vector": [0, 0, inner_total]}},
        {"command": "model.load.create", "params": {"name": "out_of_plane_ends", "type": "force", "target": {"set": "tip_ends"}, "vector": [0, end_total, 0]}},
        {"command": "model.load.create", "params": {"name": "out_of_plane_inner", "type": "force", "target": {"set": "tip_inner"}, "vector": [0, inner_total, 0]}},
        {"command": "model.step.create", "params": {"name": "In_plane", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["in_plane_ends", "in_plane_inner"]}},
        {"command": "model.step.create", "params": {"name": "Out_of_plane", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["out_of_plane_ends", "out_of_plane_inner"]}},
        {"command": "solve.static", "params": {"step": "In_plane"}, "id": "in_plane"},
        {"command": "solve.static", "params": {"step": "Out_of_plane"}, "id": "out_of_plane"},
    ]


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    print("[\n " + ",\n ".join(json.dumps(c) for c in twisted_beam(nx, ny)) + "\n]")
```

## Refinement

Each row is the generator run with the listed (nx, ny) and run with `cargo xtask oracle --case
mh-twisted-beam --no-write`; the committed recipe is the 12 × 2 row. Errors are against the reference values.

| cells (nx × ny) | u_z, in-plane | error | u_y, out-of-plane | error |
|---|---|---|---|---|
| 12 × 2 (this case) | 0.0055443 | +2.22% | 0.0016836 | −4.02% |
| 24 × 4 | 0.0055779 | +2.84% | 0.0017620 | +0.45% |
| 48 × 8 | 0.0055912 | +3.08% | 0.0017885 | +1.97% |
| 96 × 16 | 0.0055959 | +3.17% | 0.0017958 | +2.38% |

The element does not converge to the reference at γ = G/100: in-plane it settles at +3.2 % and out-of-plane
at +2.4 % at 96 × 16, and refinement does not remove the offset. Two causes, measured in issue #1077:

- **Drilling penalty on a twisted surface.** On a surface with twist rate τ the penalty adds
  (t τ)² / (6 (1 − ν) γ/G) to the plate-bending compliance at any mesh. With t τ = 0.042, ν = 0.22 and
  γ = G/100 that is 3.75 %. Measured at 96 × 16 in-plane, +3.17 % at γ = G/100 against −0.14 % at γ = G
  (the G/10 row gives +0.17 %). Raising γ locks coarse doubly curved meshes: the pinched hemisphere at 8 × 8
  goes from −1.86 % at γ = G/100 to −33.1 % at γ = G, so no single γ serves both.
- **Bilinear membrane in in-plane bending.** The plain bilinear quadrilateral carries a parasitic shear
  in its plane, with energy relative to bending (G/E)(a/H)² = 0.34 for element length a = 1 and strip
  width H = 1.1. The membrane's share of the beam's compliance is 2.0 % in-plane and 25.8 % out-of-plane,
  which gives about −0.5 % and −6.5 %. At γ = G the 12 × 2 out-of-plane result is −7.79 % and refinement
  along the span removes it (−2.78 % at 24 × 2, −1.03 % at 48 × 2).

The two 12 × 2 quantities stay recorded as known failures with the measured ratios; the case does not claim
a pass.

## Tolerance

2%, the grading MacNeal and Harder give an element within 2% of the theoretical value on a prescribed coarse
mesh (p. 14). At 12 × 2 the in-plane quantity is 1.022 of the theory and the out-of-plane one is 0.960; the
case keeps the tolerance and the reference values unchanged and records both with `known_failure`.

## Running it

```sh
cargo xtask oracle --case mh-twisted-beam
```
