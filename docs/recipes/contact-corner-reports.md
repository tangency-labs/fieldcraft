# Contact reports at a surface corner

Run the original [command recipe](contact-corner-reports.json) in a fresh model:

```sh
cargo run -p fieldcraft-cli -- run docs/recipes/contact-corner-reports.json
```

The same array can be passed as `commands` to the desktop control channel's `engine.script`.
The recipe sets mm, N and MPa display units. Two touching 10 mm cubes each have one hex8
element, E = 1000 MPa and Poisson ratio zero. Every lower-cube node is fixed, the upper cube's
lateral translations are held, and its top is displaced down by 0.1 mm in one small-sliding
nonlinear static step. This is a static equilibrium, not a physical transient. Deformation is
displayed at 50 times its actual magnitude.

The main surface is the lower cube's top. The secondary surface deliberately includes **all six
faces** of the upper cube. Only its bottom overlaps the main surface. Its four bottom nodes
should each report pressure about 10 MPa and the upper nodes zero; the integrated normal force
is about 1000 N. The small remaining penetration follows the unchanged augmented-contact
convergence criterion. The actual pressure equals the actual normal force divided by the
100 mm² interface area, independently of that stopping gap.

The reporting direction is the unit sum of overlapping face-pair normals weighted by their
nodal mortar support (the integral of the secondary shape function over that overlap).
Uncontacted side and top faces contribute no support and cannot tilt this direction. Small
sliding uses the fixed reference pair's normal; finite sliding uses the current overlap's
projection normal. At a node whose support vectors exactly cancel, there is no reporting
direction: its reported normal component is zero, and any force remains in the shear component
when frictional shear is reported. No arbitrary fallback direction is chosen.

This correction changes the decomposition reported as pressure, normal force and frictional
shear. It leaves the assembled contact force, tangent, multiplier, gap, friction return map and
convergence controls unchanged. `crates/solver/tests/contact.rs` checks this independent flat
interface identity and compares the original closing C-ring's nodal pressure and normal force
between plane strain and a one-layer 3D model. Both rings need the support-based reporting
direction: the 2D surface also turns a corner at its contacting end, while the 3D surface has
additional uncontacted Z faces. Neither model's old report is a reference solution.
