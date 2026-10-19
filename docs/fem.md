# FEM foundation: degrees of freedom and shape functions

What every element in `crates/fem` and every solver builds on. Code: `crates/fem/src/dof.rs`,
`crates/fem/src/basis.rs`, `crates/fem/src/jet.rs`.

## Degrees of freedom

Every unknown is a pair (node, `DofKind`):

| `DofKind` | Meaning | Model / deck number |
|---|---|---|
| `Ux`, `Uy`, `Uz` | translations along X, Y, Z | 1, 2, 3 |
| `Rx`, `Ry`, `Rz` | rotations about X, Y, Z | 4, 5, 6 |
| `P` | pressure (mixed u–p elements) | 8 |
| `T` | temperature | 11 |

The enum order (Ux, Uy, Uz, Rx, Ry, Rz, T, P) is the numbering order within a node.
`DofKind::from_number` / `number` convert to the 1-based numbers boundary conditions use;
`from_direction` / `direction` to the structural directions 0–5.

**Elements declare their unknowns** with `ElementDofs { node: DofSet, internal: Internal }`:
the kinds every node of the element carries (`DofSet::TRANSLATIONS` for solids, `PLANE` for
2D elements, `STRUCTURAL` for shells and beams, `TEMPERATURE` for heat transfer, translations
plus `PRESSURE` for continuous-pressure mixed elements) and the internal unknowns the element
condenses out itself (`incompatible_modes`, discontinuous `pressure`), which never reach the
global system. An element's matrices are ordered node by node, its `node` kinds in `DofKind`
order within a node (a shell's 24 × 24 stiffness: node 1 Ux Uy Uz Rx Ry Rz, node 2 …). Solid
elements: `Formulation::dofs()` / `Solid::dofs()`.

**Numbering** (`DofLayout`): collect what every node needs with `DofLayout::builder()`
(`element(labels, dofs)`, `node(label, kinds)` for constraints and boundary conditions that
need a kind no element gives, `point(label, kinds)` for reference points), then `build()`. A
node carries the union of what was asked for. Unknowns are numbered mesh nodes in label order,
each node's kinds in `DofKind` order, then the points in the order added; a slot's unknowns are
contiguous, so a node pair is one dense block in the assembled matrix. The same model always
gives the same numbering.

Lookups: `node(label)` / `slot(label)` (a slot is a position: mesh nodes `0..nodes()`, points
`point_slot(k)`), `dof(slot, kind)`, `range(slot)`, `kinds(slot)`, `locate(dof) -> (slot,
kind)`, `element_dofs(slots, kinds)`, `describe(dof)`, `stride()` (uniform layouts).

Trusses declare `fieldcraft_fem::truss::DOFS` (translations, on `Basis::Line2`/`Line3`).
Point masses, springs and dashpots are not mesh elements: the implicit solver asks the layout
for what each needs at its nodes (`fieldcraft_solver::discrete::needs`: translations for a mass
or a translational spring, rotations for an inertia or a torsional spring; a reference point
they use is a point with six) and finds its unknowns with `discrete::dof_of`.

Solvers: `fieldcraft_solver::system::Mesh::layout` numbers the implicit system (mesh nodes,
then six per coupled reference point); assembly (`sparse::add_block(k, layout, slots, kinds,
data)`, `add_scalar`, `add_diagonal`) works for any mix of kinds per node, so a solid sharing
nodes with a shell adds its translation block into the shell nodes' six-unknown blocks.
Boundary conditions and loads address (node, kind): a boundary condition holds the kinds a node
carries, and a moment needs rotations. The explicit solver has solid nodes only and checks that
its layout is the uniform node-major one (`NODE_DOFS = 3`).

## Shape functions

`Basis` is every family's interpolation:

| Family | Reference element | Nodes |
|---|---|---|
| `Line2`, `Line3` | ξ ∈ [−1, 1] | ends −1, 1; line3 adds the middle |
| `Tri3`, `Tri6` | ξ, η ≥ 0, ξ + η ≤ 1 | corners (0,0), (1,0), (0,1); mid-sides of 1-2, 2-3, 3-1 |
| `Quad4`, `Quad8`, `Quad9` | [−1, 1]² | corners counter-clockwise from (−1,−1); mid-sides of 1-2, 2-3, 3-4, 4-1; centre |
| `Tet4`, `Tet10` | unit tetrahedron | corners (0,0,0), (1,0,0), (0,1,0), (0,0,1); mid-sides of 1-2, 2-3, 3-1, 1-4, 2-4, 3-4 |
| `Wedge6`, `Wedge15` | triangle × ζ ∈ [−1, 1] | bottom 1-3, top 4-6; mid-sides of the bottom edges, top edges, vertical edges |
| `Pyramid5` | base [−1, 1]² at ζ = 0, apex (0, 0, 1) | base corners as quad4, then the apex |
| `Hex8`, `Hex20`, `Hex27` | [−1, 1]³ | corners 1-4 at ζ = −1, 5-8 above; hex20 mid-sides (bottom, top, vertical edges); hex27 face centres in `ElementKind::Hex8.faces()` order, then the centre |

Node orders are the model's (`ElementKind::edges()` / `faces()` are the single source):
corners, mid-side nodes in edge order, face nodes, body nodes. `Basis::of(ElementKind)` maps a
model element; `nodes()` gives natural node coordinates; `dim()`, `count()`, `corners()`,
`degree()` (the complete polynomial degree reproduced: 1 for linear families including the
rational pyramid, 2 for the quadratic ones).

Evaluation at a natural point `xi: [f64; 3]` (unused coordinates are zero):

- `eval(xi) -> Shape`: values `n` and first derivatives `dn[i] = ∂Nᵢ/∂(ξ, η, ζ)`. tet4, tet10,
  wedge6, hex8, hex20 use the hand-written functions of `shape.rs` (the solid kernels' path);
  the solid kernel's pyramid5 uses this module's.
- `eval2(xi) -> Shape2`: adds `d2n[i][a][b] = ∂²Nᵢ/∂ξₐ∂ξ_b`, from second-order forward-mode
  automatic differentiation of the same formulas (`jet::Jet`): exact to round-off.
- `gradients(&shape, xyz) -> Gradients`: physical `dx[i] = ∇ₓNᵢ = J⁻¹ ∇_ξNᵢ` and `det J`, with
  J_ab = ∂x_b/∂ξₐ; an inverted or collapsed element is `Error::Degenerate`.
- `hessians(&shape2, xyz) -> Hessians`: physical `d2x[i]` =
  J⁻¹ (H_ξ − Σ_b (∂N/∂x_b) X_b″) J⁻ᵀ with X_b″ = Σᵢ x_{i,b} H_ξ(Nᵢ), the map's second
  derivative (T. J. R. Hughes, *The Finite Element Method*, 2000, §3.9, §4.2), and
  `laplacian(i)`. Without the X″ term a field linear in x would have a nonzero Hessian on a
  curved element.

A family of dimension d maps onto the first d physical coordinates (plane elements in x–y,
lines along x). Elements embedded in a higher dimension (shells, beams, trusses in 3D) pass
node coordinates in their own local frame.

## Verification (`basis::tests`, `dof::tests`, `jet::tests`)

For every family: Kronecker delta at the nodes; partition of unity (Σ N = 1, Σ ∇N = 0,
Σ H = 0); exact reproduction of every monomial up to `degree()` with its gradient and Hessian
in natural coordinates; first and second derivatives against central differences;
on an affine map, physical gradients, Hessians and Laplacians exact for a complete polynomial
of the family's degree; on a curved map, linear fields exact with zero Hessian, and the
physical Hessian of a quadratic interpolant consistent with finite differences of its physical
gradient. The solver and explicit migration to `DofLayout` was checked bit for bit: every
result of the solver, explicit, engine and io test suites hashes identically before and after.

## Solver memory limits

`solve.static`, `solve.nonlinear`, `solve.heat`, `solve.frequency` and `solve.buckle`
accept `max_factor_entries`: **400,000,000 by default**, or an explicit whole number from
1 to **2,000,000,000**. For example, `solve.static` with
`{"step": "Load", "max_factor_entries": 500000000}` deliberately raises the cap for that solve;
it does not change the default of later solves. Other work and iteration caps can still only
be lowered.

Each factor value is 8 bytes: the default allows 3.2 GB of values, and the ceiling allows
16 GB. These are factor values alone, not total RAM: indices, assembled matrices, numeric
workspace and other solver storage need additional memory. Raising the cap neither reserves
memory nor guarantees that it is available. Structural and heat cap errors report the required
entry count and its factor-value byte estimate before allocating the factor values.

The same entry budget continues to guard dense constraint reduction and Lanczos storage;
32-bit builds remain limited by addressable allocation sizes. The solve parameter's ceiling
is defined once as `FACTOR_ENTRY_CEILING` in `crates/engine/src/cmd/solve.rs`.

## Heat transfer and thermal strain

A heat conduction element carries one temperature per node (`heat::HEAT_DOFS`, kind `T`);
`heat::Conduction` gives the conductivity ∫ ∇Nᵀ k ∇N dV (full rule), the capacity ∫ ρc N Nᵀ dV
and a source ∫ N Q dV (mass rule), and the heat flux −k∇T at the integration points;
`heat::face_weights` and `face_film` integrate surface flux and convection on tri3/tri6/quad4/quad8
faces. The solver (`fieldcraft_solver::heat`) numbers the temperatures with
`DofLayout::uniform(…, DofSet::TEMPERATURE)`, assembles with the same `sparse::add_block`, and
eliminates prescribed temperatures and ties (see Ties) with the structural solver's
`constraints::Reduction`.

Thermal strain ε_th = α (T − T_ref) along the material axes enters a solid as the load
`Solid::thermal_load(d, &ThermalStrain { alpha, delta })` = ∫ Bᵀ D ε_th dV and its stress as
`Solid::stresses_with(d, u, Some(&thermal))` = D (ε − ε_th); `delta` is the nodal T − T_ref.
2D elements have the same pair (`Plane::thermal_load`, `Plane::stresses_with`): plane stress
leaves εzz free (σzz = 0), plane strain holds it (σzz = −D₃·ε_th), axisymmetric elements strain
the hoop direction too. Trusses (rod2, rod3) take the axial part Σ αₖ tₖ² ΔT along their tangent (`Truss::thermal_load`,
`Truss::stresses_with`). Shells and beams will provide the same two functions for their strains;
springs, dashpots and point masses have no thermal strain.

Explicit dynamics takes a step's temperature loads as a thermal strain rate per element,
α ΔT / end time (the change reached linearly over the step, like a prescribed displacement, so
there is no thermal shock), subtracted from each increment's strain increment before the
constitutive update, so it feeds J2 and Johnson–Cook plasticity alike; isotropic expansion only
(the rate is the same in any frame). Explicit trusses take no thermal strain yet (a warning
says so).

Model side: `Material::thermal` (`conductivity`, `specific_heat`, `expansion`, each isotropic
or per axis, and `reference_temperature`), `BoundaryCondition::Temperature`, loads
`SurfaceFlux`, `Convection`, `HeatSource` (heat transfer) and `Temperature` (a temperature field
in a static step), procedures `HeatTransfer { transient }` and `CoupledTemperatureDisplacement`.
Commands: `solve.heat`, and `solve.static` for coupled steps. The deck writers refuse heat
transfer by name for now.
## Hybrid (mixed u–p) solids

A section with `"formulation": {"hybrid": true}` gives `Formulation::Hybrid { pressures }`
(`crate::hybrid`): the isotropic elasticity splits as D = D_dev + K m mᵀ, the mean stress is an
independent field discontinuous between elements (hex20: 1, x, y, z; hex8 and tet10: constant),
and it is condensed out element by element, p̂ = K M⁻¹ Gᵀ u, leaving the stiffness
K_dev + K G M⁻¹ Gᵀ. `ElementDofs.internal.pressure` counts the condensed unknowns; nothing
global changes, and Cholesky still applies. Stresses are D_dev ε + p m with the element's
pressure field. The choice of pressure space follows the inf-sup analysis (docs in
`crates/fem/src/hybrid.rs`): Q2–P1disc and P2–P0 are stable; Q1–P0 does not lock but may show
checkerboard pressures; linear tetrahedra and wedges lock even so and are refused. A thermal
strain splits the same way: ∫ Bᵀ D_dev ε_th on the nodes and, through the pressure,
K G M⁻¹ ∫ N_p tr ε_th, so a heated nearly incompressible solid takes its bulk stress −3KαΔT
when confined. Plane strain and axisymmetric quad4, quad8 and tri6 sections take one constant
pressure (the mean dilatation, the hoop strain in the divergence): K_dev + κ g gᵀ / V with
g = ∫ Bᵀ m dV (`crates/fem/src/plane.rs`); the linear triangle locks and plane stress has no
volume constraint, so both are refused. A continuous
nodal pressure (`DofKind::P`, Taylor–Hood) needs a symmetric indefinite factorization and is
planned.

## Symmetric indefinite factorization (LDLᵀ)

`fieldcraft_solver::ldlt` factors P A Pᵀ = L D Lᵀ (L unit lower triangular, D with 1 × 1 and
2 × 2 pivots) on the Cholesky factor's symbolic factorization: the same nested dissection
ordering, elimination tree, supernodes, fill, factor-size cap and parallel subtree tasks, so it
is deterministic on any number of threads. Pivots are chosen by the bounded Bunch–Kaufman (rook)
search of Ashcraft, Grimes & Lewis (1998) among the columns of the supernode being factored
(Schenk & Gärtner's supernode Bunch–Kaufman, 2006): an interchange inside a supernode changes no
row structure, so nothing is delayed and the static structure holds. Its limits, stated in the
module: L's entries below the supernode are not bounded by the search (`solve_refined` adds
iterative refinement), and a column whose supernode offers no pivot above round-off is a located
`Singular` error (a matrix singular at that column, or a zero diagonal coupled only to later
supernodes, which no FieldCraft system has: hybrid pressures are condensed per element).

The factor reports its **inertia** (positive and negative eigenvalue counts, Sylvester's law):
a 2 × 2 pivot of the search has a negative determinant (one of each). `Symmetric::new(a, cap,
kind)` picks Cholesky (`Definite`, the default everywhere a matrix is positive definite, so
linear results are bit for bit unchanged), LDLᵀ (`Indefinite`), or Cholesky with LDLᵀ when the
matrix turns out not to be positive definite (`Auto`, one symbolic factorization for both).
On a hex8 stiffness pattern of 140 000 equations (1.5 · 10⁸ entries) the LDLᵀ factor costs 1.1
times the Cholesky factor, and 1.2 times on a shift with 10 600 negative eigenvalues and 7 700
2 × 2 pivots (`ldlt::tests::speed_against_cholesky`).

## Linear buckling

A step `{"type": "buckle", "modes": 4}` (`Procedure::Buckle`, solved by `solve.buckle`,
`fieldcraft_solver::buckle`) solves its loads and prescribed motions as a linear static preload,
forms the initial-stress stiffness K_σ = ∫ (∇N_a · σ ∇N_b) I dV of every solid and 2D element
from the preload's stresses at the stiffness points (`Solid::geometric_stiffness`,
`Plane::geometric_stiffness`, with the axisymmetric hoop term), and finds the load factors λ
with (K + λ K_σ) φ = 0 by block Lanczos in the K inner product (Grimes, Lewis & Simon 1994):
without a shift on K⁻¹(−K_σ) (θ = 1/λ, Cholesky of K), the lowest positive λ; with `"shift": σ`
on the buckling transformation (K + σ K_σ)⁻¹ K (θ = λ/(λ − σ), LDLᵀ), the λ above σ, and the
LDLᵀ inertia counts the load factors between 0 and σ, reported. Negative load factors (the
reversed load buckles) are reported, not returned. Modes are scaled to a largest component of 1
and display and probe like frequency modes (their value is λ; `FrequencyResult::buckling`).
Boundary conditions hold their directions at zero in the eigenproblem; pressures are dead loads
(no load stiffness); trusses, shells and beams are refused (no K_σ yet). Decks: Abaqus
`*BUCKLE, EIGENSOLVER=LANCZOS` (refused after a general step, whose loaded state Abaqus would
buckle from, and with a shift); Nastran SOL 105 and LS-DYNA refuse it by name. Verified by Euler
columns of hex20 (pinned and clamped–free, 0.18 % from Engesser's shear-corrected load) and the
shifted form against the plain one.

**Imperfections.** A nonlinear static step's `"imperfection": {"step": "Buckle", "mode": 1,
"scale": 0.01}` solves that buckling step first and adds its mode, scaled to a largest component
of `scale`, to the node coordinates; the step is solved and reported on that geometry, with a
warning naming the mode. Abaqus and LS-DYNA decks refuse it by name (their imperfection reads a
results file). Verified: a clamped–free column with its first mode as imperfection under half its
buckling load deflects by δ₀ α/(1 − α) (2.8 %).

## Nonlinear static analysis (finite strain)

A step `{"type": "nonlinear_static", …}` (`Procedure::NonlinearStatic`, solved by
`solve.nonlinear`, `fieldcraft_solver::nonlinear`) follows equilibrium in load increments with
Newton–Raphson iterations. The step's load factor λ runs from 0 to 1: dead loads are λ times
their values, a pressure follows the deformed faces with λ times its magnitude, and prescribed
displacements are λ times theirs (a linear ramp).

**Material point** (`fieldcraft_fem::constitutive`). Given the deformation gradient F and the
history committed at the end of the last converged increment, a material returns the
**Kirchhoff stress** τ = J σ and the **spatial tangent** c of its Lie derivative (L_v τ = c : d),
in Voigt order with engineering-shear columns, so c = D at F = I. This is the push-forward of
the second Piola–Kirchhoff stress and its tangent: τ = F S Fᵀ, c_ijkl = F_iI F_jJ F_kK F_lL
C_IJKL. History variables (`MaterialModel::history_len` per integration point) are read from the
committed slice and written to a separate one, which the solver commits when the increment
converges, so a cut-back increment leaves no trace. The materials:

| Material in a nonlinear step | Model |
|---|---|
| Elastic (isotropic or orthotropic), no option | St Venant–Kirchhoff: S = D E, E = ½(FᵀF − I). Right for large rotations with small strains; it softens without bound in large compression. |
| Isotropic elastic with `"hyperelastic": {"type": "neo_hookean"}` | Compressible neo-Hookean with the isochoric–volumetric split, W = μ/2 (Ī₁ − 3) + κ/2 (J − 1)², μ and κ from E and ν (C₁₀ = μ/2, D₁ = 2/κ). Its small-strain limit is the linear material, which linear steps use. |
| Isotropic elastic with `{"type": "mooney_rivlin", "c01": …}` or `{"type": "yeoh", "c20": …, "c30": …}` | Decoupled polynomials W̄ = C₁₀ (Ī₁ − 3) + C₀₁ (Ī₂ − 3) + C₂₀ (Ī₁ − 3)² + C₃₀ (Ī₁ − 3)³ with the same volume κ/2 (J − 1)²: Mooney–Rivlin with C₁₀ = μ/2 − C₀₁ (0 ≤ C₀₁ ≤ μ/2), Yeoh with C₁₀ = μ/2, so 2(C₁₀ + C₀₁) is the elastic μ and the small-strain material stays the elastic one. Stress and tangent from the fictitious τ̄ and c̄ (Holzapfel 2000, eqs. 6.168, 6.193). Abaqus `*HYPERELASTIC, MOONEY-RIVLIN` / `YEOH`; LS-DYNA refuses them by name. |

| Isotropic elastic with `"plastic"` (bilinear isotropic or kinematic, table, Johnson–Cook) | J2 plasticity (`fieldcraft_fem::elastoplastic`): at finite strain Simo's multiplicative model (CMAME 99, 1992) with a Hencky elastic energy, τ = K ln J I + 2G dev ½ ln bₑ, the history Cₚ⁻¹ and the equivalent plastic strain, integrated by the exponential return map on the principal logarithmic strains; at small strain the additive radial return with the back stress. Both call `J2::radial_return`, the update of explicit dynamics, and its consistent tangent `J2::tangent`. |

**Plasticity's tangent.** With Cₚ⁻¹ frozen over an iteration, L_v τ = 𝔻ᵉᵖ : ½ (∂ ln b/∂b) : (d b + b d) − (d τ + τ d), 𝔻ᵉᵖ the algorithmic tangent of the return
(K I⊗I + 2Gθ ℙ − 2Gθ̄ n⊗n, Simo & Hughes 1998 Box 3.2) and the derivative of the logarithm the
Daleckii–Krein divided differences θ_AB = artanh(r)/r, r = (λ²_A − λ²_B)/(λ²_A + λ²_B), in the
principal frame of the trial bₑ: smooth through equal stretches, so uniaxial and undeformed
states need no special case. The tangent is symmetric (isotropic hardening) and checked against
central differences of the stress and of every element's forces from a plastic history. Hybrid
elements take the deviatoric return and the volume K/2 (ln J)². Refused by name in nonlinear
steps: kinematic hardening at finite strain (a back stress is not coaxial with bₑ, which the
principal-frame return needs; it works with small strain), the Johnson–Cook strain-rate term (a
static step has no time; C = 0 is the rate-independent law), and incompatible-mode elements of a
plastic material. Damage, erosion and adiabatic heating act in explicit dynamics only (a warning
says so). A hardening slope below −3G leaves the point without a tangent: the increment cuts
back and the step stops at the minimum increment, naming the element.

"Isotropic elastic" keeps meaning the same thing everywhere: the neo-Hookean model is chosen
only by the explicit option, and explicit dynamics refuses a hyperelastic material by name.

**Elements** (`fieldcraft_fem::finite`), total Lagrangian over the reference element with its
own quadrature: f_a = ∫ Bₓ,aᵀ τ dV₀ and K_ab = ∫ (Bₓ,aᵀ c Bₓ,b + ∇ₓN_a·τ∇ₓN_b I) dV₀ with the
current gradients ∇ₓN = F⁻ᵀ∇₀N (material plus geometric stiffness), identical term by term to
∫ B₀ᵀ S dV₀ and its linearisation. Every element's tangent is checked against central
differences of its forces.

| Element | Nonlinear steps |
|---|---|
| tet4, tet10, wedge6, hex8, hex20 (full), hex20 (reduced) | total Lagrangian |
| hybrid quad4, quad8, tri6 (plane strain, axisymmetric) | mean dilatation: θ = ∫ J dV₀ / V₀ with the hoop stretch, p = U′(θ), stiffness U″ g gᵀ / V₀ |
| hybrid hex8, hex20, tet10 | Simo–Taylor–Pister three-field form on the linear hybrid elements' pressure fields (θ the projection of J, p the projection of U′(θ)); needs the neo-Hookean option (St Venant–Kirchhoff has no volumetric split) |
| hex8 and quad4 incompatible modes | refused: their finite-strain form (enhanced assumed strain, Simo–Armero 1992) needs the modes as unknowns carried between iterations and hourglasses in large compression; use full integration, hybrid or quadratic elements |
| plane strain, axisymmetric (tri3, tri6, quad4, quad8) | total Lagrangian with F₃₃ = 1 or r/R (hoop column N/r, geometric term τ_θθ N_a N_b/r²) |
| plane stress | refused (σ_zz = 0 needs a local iteration for the thickness) |
| trusses, shells, beams | refused (large rotations of structural elements come later in M15) |

Without nonlinear geometry (`"nonlinear_geometry": false`) an elastic element is the linear one
(K u), and an element of a plastic material takes the same loops with small strain
(`Kinematics::Small`: ε = sym ∇u, the reference gradients, f = ∫ Bᵀσ dV, K = ∫ Bᵀ𝔻B dV with the
consistent tangent, no geometric stiffness; hybrid elements project 1 + tr ε and K tr ε), so an
elastic model is still the linear solution reached by Newton.

**Follower pressure** (`fieldcraft_fem::face::follower_pressure`, `follower_edge_pressure`):
fᵢ = −p ∫ Nᵢ (x,ξ × x,η) dξ dη on the current face, with its load stiffness
K = p ∫ N_a (N_b,η [x,ξ]× − N_b,ξ [x,η]×). The load stiffness is not symmetric. The solver keeps
its symmetric part and leaves the skew part out, so the tangent stays symmetric for the
supernodal Cholesky factorization: the residual is exact (the converged answer does not
change), and over a closed surface or one whose boundary is held the skew parts of neighbouring
faces cancel (Schweizerhof & Ramm 1984), so nothing is lost there; elsewhere Newton may converge
linearly instead of quadratically, and the result says how large the left-out part is. An LU
factorization would double the factor's memory and lose the symmetric solver for a gain only in
the rate of convergence of open loaded surfaces; it is revisited with the symmetric indefinite
(LDLᵀ) factorization that arc-length and Taylor–Hood need. Without nonlinear geometry a pressure
is a dead load.

**Load path.** A step's `load_path`, `[[fraction, factor], …]` from `[0, 0]` to fraction 1, makes
the load factor piecewise linear in the step fraction (without one it equals the fraction), so
loads and prescribed displacements can come down or reverse within a step: cyclic loading.
Increments end at every point of the path, and the extrapolated predictor scales with the change
of the load factor, so it turns with the path. LS-DYNA writes the path as the step's curve;
Abaqus refuses it by name for now (it would be an amplitude on every load and boundary
condition).

**Predictors.** An increment starts from the linear extrapolation of the last one; with
path-dependent materials an increment without one (the first, after a turn of the path to zero
change) starts from the tangent predictor: one solve with the converged state's tangent for the
change of the loads and of the prescribed displacements, so a prescribed motion is spread into the
interior before the material sees it (started from the boundary's jump alone, the elements next
to it would flow plastically far beyond the solution). Its work is then the energy criterion's
reference. A first tangent that is not positive definite is a cutback when the undeformed
model's tangent is (the strain of the increment, not a model that can move freely).

**Newton–Raphson** with the consistent tangent refactored every iteration, an optional line
search (off by default: on large elastic rotations the full step converges in fewer
iterations), and convergence when the force residual (∞-norm over the largest force seen), the
displacement correction (over the increment's largest displacement) and the energy (over the
first iteration's) are all below their tolerances, or the residual is round-off. Rotations and
their moments are measured apart from translations and forces (see *Convergence norms* below).
Each increment
starts from the linear extrapolation of the last one. Automatic incrementation grows an
increment by half after two fast ones (≤ 4 iterations), shrinks it after a slow one (> 8), and
cuts it to a quarter when an element inverts, the residual is not finite or grows twice, the
tangent is singular, or the iterations run out; the step stops with the load factor
and the reason at the minimum increment, the increment cap, or `Limits::max_newton_iterations`.

**Convergence norms.** The residual and the correction are measured per kind of freedom, each
against its own scale, so the test is the same in any consistent unit system and at any size of
the part (a moment scales as a force times a length, a rotation not at all):

| Check | Measured | Against |
|---|---|---|
| force | largest residual force | the largest force seen (applied, internal, reaction), at least 10⁻⁸ of the mean diagonal stiffness times L |
| moment | largest residual moment | the largest moment seen, at least the force scale times L |
| translation | largest translation of the correction | the increment's largest translation (or prescribed one), at least 10⁻¹² L |
| rotation | largest rotation of the correction, in radians | the increment's largest rotation (or prescribed one), at least the translation reference over L |

L is the diagonal of the model's bounding box. The force and moment checks share the `force`
tolerance, the translation and rotation checks the `displacement` tolerance, and each iteration
records the larger measure of each pair. The energy criterion's work, s δq · r, is in N·m for
both kinds and stays one sum. Today only coupled reference points carry rotations in a nonlinear
step (shells and beams are refused); a model without them converges exactly as with one norm.
In the round-off floor of the force scale a rotation's diagonal stiffness counts over L² (the
stiffness of the translation it gives at the model's size).

**Indefinite tangents.** A tangent that is not positive definite (a strained nearly
incompressible mixed element, an unstable material, a state past a limit or bifurcation point)
is factored by LDLᵀ and Newton goes on; each iteration records the tangent's negative
eigenvalues, and the summary warns when a converged state keeps some (an unstable equilibrium
under load control). This lifts the old limit on hybrid elements near ν = ½ under prescribed
displacement (at ν = 0.4999 every increment used to be cut back to the cap). Only the first
tangent of a model that can move freely is an error. The tangent predictor of path-dependent
materials keeps Cholesky and is skipped when the tangent is indefinite.

**Arc length.** `"arc_length": {"constraint": "cylindrical"|"spherical", "max_load_factor": 1,
"displacement_limit": {"node", "direction", "displacement"}}` makes the load factor an unknown
(Riks 1979, Crisfield 1981): each increment keeps ‖Δq‖² + w Δλ² = Δl² (w = 0 cylindrical, w the
first tangent displacement's ‖q_f‖² for the spherical arc; the norm weighs a rotation by L², the
model's size squared, so the arc is a length² in any unit system and at any size of the part,
and a step driven by rotations is still measured by them), predicts along the tangent with the
sign of its angle to the last increment (Feng, Perić & Owen 1996), and corrects with
δq = δq_r + δλ δq_g, δλ the root of the constraint's quadratic whose increment turns least
(Crisfield's angle criterion). The first increment's load factor change is `initial_increment`,
which sets the arc; arcs grow, shrink and cut back as increments do, within the minimum and
maximum increments' arcs. Prescribed motions scale with λ (∂r/∂λ includes −K T_p u_p) and
follower pressures with their current faces. The tangent at every converged state is factored
for the next predictor (LDLᵀ when indefinite); its negative eigenvalues are recorded per
increment (`negative_pivots`), and where their count changes the step reports a **critical
point**: a limit point where the load factor turns, a bifurcation otherwise
(`History::critical_points`, in the `solve.nonlinear` summary). The step ends at the maximum load
factor, at the displacement limit, or at its increment cap (a warning). Refused: a load path,
temperatures. Decks: Abaqus `*STATIC, RIKS` (cylindrical only); LS-DYNA refuses arc-length steps
by name. Verified by the snap-through of a shallow sinusoidal arch of hex20 solids
(`crates/solver/tests/arc_length.rs`, `tests/oracle/arc-shallow-arch`): the path within 1.9 % of
the one-mode closed form, both limit points, and the inertia 0 → 1 → 0 exactly between them.

**Results**: the final displacements and reactions, the **Cauchy stress** and the
**logarithmic strain** (ln V, engineering shears; the small strain without nonlinear geometry) at
the integration points and averaged at the nodes, the **equivalent plastic strain** at the
integration points of plastic elements (`solve.probe` with `elements`, the largest in the step
summary), every increment's load factor, iterations and residuals, and a displacement frame per
converged increment (`results.frame`, the slider in the results panel).

**Couplings** stay the linear constraints of the linear solver (small rotations of the
reference points), with a warning; springs act along their fixed directions. Centrifugal loads
are refused by name in nonlinear steps.

**Thermal strain.** A nonlinear step's temperatures are its temperature loads, or the steady
conduction of its temperature boundary conditions and heat loads solved first (one-way, as the
linear coupled temperature–displacement step); they ramp with the load factor from each
material's reference temperature. With nonlinear geometry F = F_m F_θ, F_θ = diag(1 + αᵢ ΔT) at
each integration point (V. A. Lubarda, Appl. Mech. Rev. 57:95–108, 2004): the material sees F_m,
the loops run on the expanded volume J_θ dV₀ (mixed elements project det F_m), and the reported
strain is the total logarithmic one; the linear stretch 1 + α ΔT keeps α the engineering
coefficient, so a free body grows by exactly α ΔT. Without nonlinear geometry the split is
additive (ε_m = ε − α ΔT), and elastic elements on their linear stiffness take the equivalent
nodal forces as a dead load. Johnson–Cook plasticity softens at the point's temperature. The
expansion is constant in temperature (the model has no temperature-dependent coefficient yet).

**Decks**: Abaqus `*STEP, NLGEOM=YES|NO, INC=` with `*STATIC` (initial, period 1, minimum,
maximum) and `*CONTROLS, PARAMETERS=LINE SEARCH`; the neo-Hookean option as `*HYPERELASTIC, NEO
HOOKE`; plasticity as for explicit dynamics (`*PLASTIC`); a load path, and a plastic deck of
several general steps (Abaqus carries the plastic state from one to the next, every model step
starts undeformed), are refused. LS-DYNA: the ramped implicit run with full Newton (`ilimit` 1, `maxref`) and
`*CONTROL_IMPLICIT_AUTO`; hyperelastic materials in nonlinear steps and steps without nonlinear
geometry are refused; a load path is the step's curve. Nastran refuses nonlinear steps (SOL 400) by name. Iteration caps and
tolerances are left to each solver's own controls, whose norms differ from ours.

**Verification** (`crates/solver/tests/nonlinear.rs`, `crates/fem/src/{constitutive,finite,face}.rs`,
and `tests/oracle/nl-*`): consistent tangents of both materials and of every element (central
differences); uniaxial stretch to λ = 2 of a neo-Hookean block in every element kind, plane
strain and axisymmetric (closed form); St Venant–Kirchhoff uniaxial stretch; the nonlinear patch
test on a distorted mesh (all solids, hybrid ones included); a 90° rigid rotation (zero stress);
the elastica (Bisshopp–Drucker) at P L²/(E I) = 1, 2, 5, 10 with hex20 (0.2–0.5 %); the inflation of a thick neo-Hookean tube by a follower pressure (Rivlin's closed form);
Newton's quadratic convergence; the linear limit; a finite-strain Cook's membrane is pending
its transcribed reference; and the failures (an element inverted by
crushing stops at the minimum increment, naming it; the caps; what is refused). Plasticity
(`crates/solver/tests/plasticity.rs`, `crates/fem/src/elastoplastic.rs`, `tests/oracle/`
`nl-j2-uniaxial`, `j2-bauschinger`, `plastic-thick-cylinder`, `plastic-hinge-cantilever`): the
consistent tangents against central differences from a plastic history (material point and
every element, finite and small strain, full and hybrid); uniaxial tension to a logarithmic
strain of 0.2 in every element kind, hybrid ones included (Kirchhoff stress exact to 10⁻⁷, and
0.13 % from explicit dynamics' hypoelastic Cauchy stress, the elastic volume change between the
two laws); the Bauschinger effect on a load path back to zero strain (exact); Hill's thick
cylinder, the bore displacement at a plastic radius of 1.5a (0.009 %) and the collapse pressure
(2·10⁻⁸); a cantilever's collapse 3.1 % above beam theory's hinge, converged under refinement;
Newton's order with the consistent tangent 1.9 (small strain) and 2.3 (finite strain) in a
plastic cantilever; and the refusals and the softening material that stops at the minimum
increment. Mooney–Rivlin and Yeoh (`nl-mooney-rivlin-modes`, `nl-yeoh-uniaxial`,
`nl-mooney-rivlin-tube`): tangents by central differences (material point and hybrid
elements), Treloar's uniaxial, equibiaxial and pure-shear modes against the compressible closed
form (exact), and a tube inflated by a follower pressure while stretched axially to λ_z = 1.5
(Rivlin's incompressible solution, 0.6 % at ν = 0.4995). Thermal strain (`nl-thermal-expansion`,
`j2-heated-bar`, `crates/solver/tests/nonlinear_thermal.rs`): free expansion and a fully held
block at a thermal stretch of 1.2 (exact, full and hybrid), a bar held at its ends heated past
yield (exact), the linear limit against the linear and coupled solvers, and tangents with a
temperature field.

crushing stops at the minimum increment, naming it; the caps; what is refused).

## Contact and friction (explicit dynamics)

**Model** (`fieldcraft_model::contact`): a contact names its surfaces as sets of element faces —
a **pair** (`main`, `secondary`) or **general** contact on one surface (self contact) or on every
exterior face of the solids — with Coulomb **friction** μ = μ_k + (μ_s − μ_k) e^{−d γ̇} (static,
kinetic, decay per slip rate), a **penalty scale** s ≤ 1 and what happens to nodes that start
inside the other surface (**adjust**: moved onto it strain-free and reported; **refuse**: the
analysis stops naming the node, the face and the depth). Planar rigid walls take the same friction.
Commands `interaction.contact.pair|general|edit|delete|find`; Abaqus writes `*SURFACE
INTERACTION` with `*FRICTION` (exponential decay), `*CONTACT PAIR` and general `*CONTACT`, LS-DYNA
`*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE` and `_SINGLE_SURFACE` (`fs`, `fd`, `dc`, `sfsa`,
`sfsb`); Nastran and the linear implicit analyses refuse contact by name; nonlinear static steps
take contact pairs (next section).

**Geometry** (`fieldcraft_fem::contact`, shared with the implicit solver to come): the closest
point of a triangle (exact region tests) or of a bilinear quadrilateral (projected Newton with the
exact Hessian) to a point, with the signed gap along the face's normal and the lateral offset;
a bucket grid of face boxes whose cell is twice the mean box, grown until the grid and its entries
fit their caps, queried by point or by box with the candidates in face order.

**Kinematics** (`fieldcraft_explicit::contact`): segment to segment. Each secondary face and
each main face facing it (normals within 84° of opposite) are projected onto the secondary face's
plane, split into triangles and clipped against each other; on each overlap triangle a
three-point rule takes the gap g = (x_s − x_m)·n_m between the faces' closest points. Where
g < 0 (no deeper than half the main face's element), a penalty pressure p = ε(−g) acts on the
secondary face's nodes with weights N_i^s along n_m and on the main face's nodes with −N_j^m. For
a uniform pressure across flat faces both meshes receive exactly the consistent nodal forces,
so the contact patch test passes on non-matching meshes. Every contact runs both ways (a pair
with each surface as the secondary, single-surface contact with every face against every other
face that shares no node with it) at half the penalty each.

**Penalty**: ε = s γ min_i(m_i/A_i)/Δt² over the contact's nodes (m_i lumped mass, A_i the node's
area facing one way, Δt the stable increment, γ = 0.3). By Gershgorin's and Weyl's inequalities
the contact adds at most 2 s γ/Δt² to the square of any frequency, which the default increment
scale 0.9 leaves room for: contact never shortens the increment.

**Friction**: at the secondary nodes, by a regularised return map against the main face that
pressed each node most: the elastic predictor t* = P t_n − k_t Δu_T (k_t = ε A_i/2), stick when
|t*| ≤ μ f_n, otherwise slip back to the cone, dissipating μ f_n (|t*| − μ f_n)/k_t. On rigid
walls, a node's tangential velocity changes by at most μ times its normal velocity change.

**Energy**: the contact forces' work (penalty energy plus friction dissipated) is integrated with
the trapezoidal rule and enters the balance; friction dissipation (contact and walls) is
reported apart. **Search** every increment, in parallel by faces, forces added in face order: the
results are the same on any number of threads.

**Verification** (`crates/explicit/tests/contact.rs`, `tests/oracle/explicit-*`): two bars
colliding (contact for 2L/c, interface force ρcAv₀/2, velocities exchanged, kinetic energy kept to
0.01%, the contact energy back to 10⁻⁴ of it), a block sliding to a stop under constant and
decaying friction on a rigid wall and on a deformable plate (distance and time to 0.1–0.5%), the
contact patch test on non-matching meshes (σ_zz uniform to 0.1%), two spheres colliding against
Hertz (peak force +2.3%, duration −0.4% on a graded quarter-sphere mesh), a free body untouched
by general contact on all its faces, initial penetrations adjusted or refused, the face cap, and
bit-for-bit results on one and four threads.

## Contact in nonlinear static steps (mortar, augmented Lagrangian)

**Discretisation** (`fieldcraft_solver::contact`): mortar contact, one pass, the pair's
`secondary` surface the non-mortar side. Each secondary node A carries a pressure p_A and the
weighted gap g̃_A = Σ ∫ Φ_A (x_m − x_s)·n dA over the overlaps of its faces with the main faces
facing them, Φ_A the secondary face's own shape functions; ĝ_A = g̃_A / A_A with A_A = ∫ Φ_A over
the node's faces. The pressure field is interpolated like the displacements, so a uniform
pressure crosses a non-matching interface as the consistent nodal forces of both meshes: the
contact patch test. General (single-surface) contact uses one face list, a set's or every
exterior face of the solids, as both sides: each face against every other one facing it that
shares no node with its element (its neighbours) and that it does not reach through more than
half an element (the far side of a thin wall), paired both ways, so the folds of one body press
on each other symmetrically. Linear faces (tri3, quad4: tet4, wedge6, hex8), and in 2D models
(plane strain, plane stress, axisymmetric) the linear edges of tri3 and quad4 elements (sets'
edges; general contact defaults to every exterior edge); quadratic faces and edges and
arc-length steps are refused by name.

**2D models**: each edge is oriented with its normal (t_y, −t_x) out of its element; the main
edge is projected along the secondary edge's normal onto its line, which maps main points to the
secondary parameter linearly, so the overlap is an interval integrated with three-point
Gauss–Legendre (exact for the degree-3 Φ_a N_b r of the axisymmetric patch test). The
out-of-plane measure (the thickness, or 2πr at the secondary point) is in the weights, so A_A,
g̃_A and the forces are totals over the slab or the ring, as the 2D elements' are (B. Yang,
T. A. Laursen, X. Meng, IJNME 62:1183–1225, 2005). Points keep z = 0 and the same 24 variables
as faces, and the tangent blocks keep their in-plane rows. The element's thickness under an edge
is its area over the edge's length.

**Segmentation** (`fieldcraft_fem::mortar`): both faces are projected along the secondary
face's unit normal onto its plane, the main polygon is clipped by the secondary one
(Sutherland–Hodgman; a projected quadrilateral that is not convex goes as its two triangles), and
the overlap's fan of triangles is integrated with Radon's seven-point degree-5 rule. At each point
both faces' coordinates come from inverting their projected maps (linear, or bilinear by the root
of its quadratic), so the main point is on the main surface along the normal, warped or not. A
vertex within 10⁻⁸ of a clip edge's length of its line is on it and never cut: the intersection
of two nearly coincident edges is ill-conditioned and its derivative swamped by round-off. The
kernel is written over a generic scalar: `f64` for values, forward-mode `Dual` for gradients and
`Hyper` (gradient and packed Hessian in the 24 corner coordinates) for the consistent tangent,
the exact linearisation of the segmentation, normal and projection (`fieldcraft_fem::dual`).

**Enforcement**: augmented Lagrangian. With the multipliers λ_A fixed, p_A = ⟨λ_A − ε_A ĝ_A⟩,
the derivative of Π = Σ A_A/(2ε_A) ⟨λ_A − ε_A ĝ_A⟩²: internal forces −p_A ∂g̃_A/∂x and the tangent
Σ_closed (ε_A/A_A) ∂g̃_A ⊗ ∂g̃_A − p_A ∂²g̃_A (the second term the contact's geometric stiffness).
Once Newton converges an increment, λ_A ← p_A (Uzawa) and the increment is solved again from its
solution, until every closed node's gap is within 10⁻⁵ of the thickness of the element under it
and no open node keeps a multiplier (at most 25 augmentations, then a warning). The penalty per
length is ε_A = s · 10 · 3Ḡ/h (Ḡ the mean shear modulus, h = V/A of the element under the face,
s the contact's `penalty_scale`): about ten times the local element stiffness, and 3Ḡ rather than
the constrained modulus keeps it so for nearly incompressible rubber.

**Active set and Newton**: the max() decides at every evaluation which nodes are closed
(semismooth Newton); an increment whose active set alternates between two states three times is
cut back (chattering). The tangent's pattern takes the blocks of the closed nodes' face pairs on
top of the mesh's. With contact, the force scale of the convergence test is that of the
converged states, the loads and the current state rather than the largest force of every trial
(a trial can press a body deep into another before the active set settles). The augmentations'
re-solutions are measured against the increment's start and its first correction; the next
increment's size follows the increment's own iterations.

**Small and finite sliding**: without nonlinear geometry the face pairs (those within a main
element's thickness) and their mortar integrals D_ai = ∫ Φ_a N_i, M_aj = ∫ Φ_a N_j are fixed in the
reference configuration, so g̃_a = n · (Σ M_aj x_j − Σ D_ai x_i) is linear: exact and quadratic-free.
With nonlinear geometry the pairs are searched (bucket grid of the main faces' boxes grown by
their depth limit) and integrated on the current configuration at every evaluation; a pair deeper
than half the main element's thickness has passed through and is left out.

**Friction**: Coulomb's law at the closed secondary nodes by a return map (Wriggers 2006, §5.2;
Simo & Laursen 1992). The node's weighted tangential slip since the last converged increment is
w_A = B_A Δx with B_A = [∫ Φ_A N_i P, −∫ Φ_A N_j P] (P = I − n n; the rate form of the mortar slip,
McDevitt & Laursen 2000), ĝ_T = w_A / A_A; the trial traction t* = P t_n − ε_A ĝ_T sticks when
|t*| ≤ μ p_A and otherwise returns to the cone, t = μ p_A t*/|t*|, with the static coefficient for
a node that stuck in the last increment and the kinetic one for one that slipped (a static step
has no slip rate for the exponential decay). The traction does the virtual work t · B_A δx, so a
slipping interface carries exactly μ times its normal force. The tangent is symmetrised: it
keeps (ε_A/A_A) Bᵀ B in stick and (ε_A μ p_A/(A_A |t*|)) Bᵀ (I − t̂ t̂) B in slip and leaves out
the coupling to the pressure (μ t̂ ⊗ ∂p_A/∂x) and the change of B with the geometry, so the
factorization stays the symmetric one, like the follower pressure's; the residual is exact, so
only the rate of convergence suffers while nodes slip. Tractions and slip states are committed
when an increment converges. If an increment's active set chatters (a node touching with no
pressure, at the edge of a contact that lifts under friction), the increment is cut back and the
nodes open by less than 10⁻³ of their element's thickness keep their penalty in the tangent for
the rest of the step.

**Initial overclosure**: nodes whose gap starts negative are reported; `adjust` takes the overlap
out of their weighted gap (an offset: strain-free, nodes are not moved), `refuse` stops naming
the node and the depth.

**Output**: per secondary node the pressure (its normal contact force over the area it carries,
∫ Φ_A over the overlaps, which is smoother than the multipliers at the contact's edge), the
multiplier, the gap ĝ, the status (open, closed; with friction stick or slip), and with friction
the shear traction and the slip over the last increment; per contact the force on the secondary
surface, the normal force Σ p_A ∫Φ_A, the deepest penetration and the augmentations.
`solve.nonlinear` summarises each contact, `solve.probe` gives a secondary node's `contact`
entries, and the results display draws `contact_pressure`, `contact_gap`, `contact_status`,
`contact_shear` and `contact_slip`.

**Verification** (`crates/solver/tests/contact.rs`, `tests/oracle/nl-contact-patch`,
`tests/oracle/nl-hertz-cylinder`, and their 2D forms `nl-contact-patch-plane-strain`,
`nl-contact-patch-axisymmetric`, `nl-hertz-cylinder-plane-strain`): the contact patch test on
non-matching hex8, tet4 and (plane strain and axisymmetric) quad4 meshes
(uniform σ_zz, pressure equal to it, no lateral stress: to 10⁻⁹ with small sliding, 10⁻⁴ with
finite sliding, where the gaps the augmentation leaves enter through their geometric terms);
Hertz's line contact, plane strain (peak pressure within 0.3%, the pressure at two nodes within
3% and 6%, the vertical stress under the contact within 2%, the half-width bracketed); Hertz's
sphere on a flat, a quarter model (peak pressure within 4%, contact radius bracketed; heavy);
Newton's quadratic rate near the solution with finite sliding; Coulomb friction between pressed
blocks (stick well below the cone with the forces of glued blocks, slip on it carrying μ N to 0.1%
with small sliding and 3% with finite sliding, where the block turns a little, and μ_k N once
slipping; `tests/oracle/nl-friction-sliding-blocks`); initial overclosure adjusted and
refused; bit-for-bit results on one and four threads; the refusals. In 2D: the patch tests pass
to 10⁻⁹ (small sliding) and 10⁻⁴ (finite), Hertz's line contact on quad4 edges matches the
hex8 model to 10⁻⁵, and a hybrid Mooney–Rivlin C-ring closing on itself in plane strain moves as
the one-layer solid ring with z held, node for node to 10⁻⁵.

## Ties (non-matching meshes)

**Model** (`fieldcraft_model::tie`): a tie names a `main` and a `secondary` surface (sets of
element faces) and an optional position tolerance; commands `interaction.tie|tie.edit|tie.delete`.
Abaqus writes `*TIE, ADJUST=NO` (secondary surface first, `POSITION TOLERANCE` when given),
LS-DYNA `*CONTACT_TIED_SURFACE_TO_SURFACE` between segment sets; Nastran refuses ties by name.

**Constraints** (`fieldcraft_solver::tie`): dual mortar. Each secondary face's dual shape
functions ψ_a = Σ_b A_ab N_b, A = D M⁻¹ (M_ab = ∫ N_a N_b, D = diag ∫ N_a), are biorthogonal to its
own, so the secondary mortar matrix is diagonal and each secondary node becomes the weighted
average of the main nodes under it, u_A = Σ_J (M_AJ/Σ_K M_AK) u_J with M_AJ = ∫ ψ_A N_J over the
overlaps of the reference configuration (`fieldcraft_fem::mortar`). The rows go to the constraint
elimination like a coupling's, the secondary node the intended dependent, in linear static,
frequency, buckling and nonlinear steps (small-displacement constraints). Dividing by the sum
keeps rigid motions exact at a surface's edge where a face is only partly covered. Main faces
count when the secondary face overlaps them along its normal with a mean gap within the position
tolerance (default half the main element's thickness); uncovered secondary nodes stay free with a
warning. Explicit dynamics refuses ties by name.

**Heat** (`fieldcraft_solver::heat`): the same weights tie temperatures, T_A = Σ_J w_AJ T_J (one row
per secondary node, `tie::rows_of` with the temperature kind), in steady and transient heat
transfer and in the heat solution of coupled temperature–displacement steps. The heat solver
eliminates them with the prescribed temperatures through the same transformation
(`constraints::Reduction`), so the interface has no resistance, the heat leaving one side enters
the other, and the reactions at prescribed temperatures (T_pᵀ (K T + C Ṫ − f)) balance the heat
supplied.

**Verification** (`crates/solver/tests/contact.rs`, `tests/oracle/tie-patch`,
`tests/oracle/tie-heat-series`): a linear field across non-matching hex8/tet4 meshes reproduced to
round-off (to 10⁻⁶ on distorted quadrilateral faces, where the overlap integration is accurate to
its degree) in linear and nonlinear static steps; uniaxial tension across a tie, σ_zz = Eε in
every element to 10⁻⁹. Heat: steady conduction between held ends through tied hex8/tet4 blocks
(matching, graded and mixed meshes), the linear profile, the flux and the heat flow k A ΔT / L to
10⁻⁹; with conductivities k₁ ≠ k₂ the series-resistance interface temperature
(k₁ T₀ + k₂ T₁)/(k₁ + k₂) and heat flow ΔT/(L₁/k₁ + L₂/k₂); a transient tied bar equal to the
single-body bar; a coupled step with a tie at a uniform ΔT expands freely (u = α ΔT x, no stress).
Ties take linear faces only (tet4, wedge6, hex8), for heat as for displacement.

## Hyperelastic materials in explicit dynamics

Hyperelastic materials (neo-Hookean, Mooney–Rivlin, Yeoh) use the shared finite-strain material
point (`fieldcraft_fem::constitutive`) in explicit dynamics too, in total form, which is objective
by construction (no stress rate): the one-point hexahedron takes F = I + ∂u/∂X, its mean over the
undeformed element, and the Cauchy stress of the nodal forces τ/J; the tetrahedron takes the
isochoric stress of its own F and the **average nodal pressure** of a decoupled material (each
node's U′(J_a) from its tetrahedra's volume ratio, the element's the mean of its nodes': J. Bonet,
A. J. Burton, Commun. Numer. Methods Eng. 14:437–449, 1998, at finite strain), which does not lock
when the material is nearly incompressible. (Averaging J into an F-bar instead was unstable on
the Kuhn mesh of the tests.) The stable increment follows the tangent: each element's frequency
bound is scaled by the ratio of its tangent's largest normal stiffness, plus the stress, to the
small-strain one. Logarithmic strain is the strain output. Hyperelastic elements need nonlinear
geometry and take no thermal strain yet (both refused by name). Verified by uniaxial strain to
λ = 1.5 of all three laws on hex8 and tet4 against the material point (`crates/explicit/tests/hyperelastic.rs`)
and by `tests/oracle/explicit-neo-hookean-implicit`, a quasi-static explicit pull against the
implicit nonlinear solver.

## Mass scaling (explicit dynamics)

Opt-in per explicit step (`mass_scaling: {target_increment, set}`): every solid element (of the
set, or of the model) whose stable increment times the step's scale is below the target has its
density raised by (target / increment)², which brings its increment to the target (its lumped
masses, frequency bound and hourglass and bulk-viscosity scales follow; adiabatic heating keeps
the physical density). The run reports the elements scaled, the added mass in all and per
section, its fraction of the model's mass, a warning, and the kinetic energy of the added mass
(`energy.added_kinetic`, part of the kinetic energy, so the balance holds). Verified on a uniform
bar scaled to twice its increment: added mass three times the bar's, the increment doubled, the
added mass's share of the kinetic energy three quarters (`crates/explicit/tests/one_d.rs`).
