# Filleted L-bracket on the octree all-hex mesh

The corpus part `filleted-bracket` (an 8 mm thick L, 60 × 50 × 40 mm, an R6 fillet inside the
corner and a Ø8 hole through each leg), built with SolidCraft from `design.json` and meshed by
`mesh.hex.octree` (hex8 with incompatible modes, steel E = 210 GPa, ν = 0.3). The bottom face of
the horizontal leg (`F2:side:l1`) is clamped and the top face of the upright leg (`F2:side:l5`)
carries 2 kN along +x.

Checked on the 5 mm mesh: equilibrium; the certified validity of every element; the mesh volume
against the CAD volume; the bending stress of the upright leg at y = 19 on both faces,
σ_yy = M c / I = 2000 · 31 · 4 / 1706.7 = 145.3 MPa (M is exact by statics, the plane-section
stress is the reference, the corner and a hole are near). The von Mises stress at the middle of
the fillet, which has no closed form, is checked for convergence against the 4 mm mesh:

| mesh size | hexahedra | fillet von Mises (MPa) | σ_yy outer / inner at y = 19 (MPa) |
|---|---|---|---|
| 7 mm | 4291 | 137.5 | 149.3 / −146.7 |
| 5 mm | 13277 | 144.5 | 149.3 / −147.1 |
| 4 mm (reference) | 21508 | 142.3 | 147.4 / −145.5 |

A 3 mm mesh needs more than the direct solver's factor cap.
