# Surface foundations

`model.foundation.create` and `.edit` take a complete definition; `.delete` takes `name`.

```json
{"command":"model.foundation.create","params":{"name":"Bed","set":"bottom","normal":"100 N/mm^3","tangential":"50 N/mm^3","normal_damping":"0.2 N*s/mm^3","tangential_damping":0,"loss_factor":0.02}}
```

Stiffness is force per volume (pressure per displacement), stored in N/m³.
Viscous coefficients are stored in N·s/m³; loss factor is dimensionless. Bare
numbers use the display system; explicit quantities pass through the unit boundary.
The reply identifies the definition by name. Edit replaces all fields; omitted
optional coefficients reset to zero. Coefficients must be finite and nonnegative.
The set must exist and `model.check` requires that it resolve to element faces.
Definitions act in every structural step and persist with the model; a referenced set cannot be deleted. They have no per-step activation switch.

At each reference-face integration point, n is the unit normal. The local operator
is k_n nnᵀ + k_t(I−nnᵀ), integrated consistently with the face shape functions.
The tangential coefficient is equal in both in-plane directions. Reversing the
normal does not change the matrix. Foundations are linear and bilateral; they do
not represent lift-off, friction or a finite-thickness substrate. Shells and
composite shells use reference-surface translations, with no added rotational
spring or thickness-offset lever arm.

Linear statics and undamped natural frequencies include foundation stiffness.
Direct and modal harmonic response include C_f from the two viscous coefficients
and H_f = loss_factor K_f. For exp(iωt), Z=K−ω²M+i(ωC+H).
Step Rayleigh damping applies to total stiffness, including foundations; material
Rayleigh damping still applies only to the material elements. Nonlinear and
explicit dynamics refuse foundations. Planar/axisymmetric foundations are refused.

The design parameter is `{"type":"foundation","name":"Bed","component":"normal"}`
or `"component":"tangential"`. Static, simple-mode frequency and direct harmonic
sensitivities use exact dK/dk, including at zero stiffness; loss-factor stiffness
derivatives are retained. Joint identification accepts these parameters together
with material parameters. For static observations, each trial shares one baseline
factorization per step across parameters and observation points; tangent vectors
retain the same assembled derivatives as scalar sensitivities. This is not a
uniform-material scaling approximation. Repeated eigenvalues retain the existing sensitivity
refusal; modal harmonic sensitivities remain outside this branch's scope.

Abaqus export uses `*FOUNDATION` for normal-only undamped support. General linear
support writes consistent sparse K/C/H using `*MATRIX INPUT` and `*MATRIX ASSEMBLE`,
with shared physical nodes, global translation DOFs and a lower-triangle format.
It does not also write the same stiffness as `*FOUNDATION`. Matrix export supports
linear static, frequency and direct harmonic procedures; modal-superposition harmonic
matrix export is refused until the writer explicitly enables the required SIM architecture.
Nastran and LS-DYNA writers refuse foundations by name. Text export tests do not
constitute evidence that a commercial solver has accepted or solved the deck.

References: Hetényi, *Beams on Elastic Foundation* (1946); Zienkiewicz, Taylor & Zhu,
*The Finite Element Method* (2013), virtual work; Abaqus 2025 public Keywords Guide,
[FOUNDATION](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-foundation.htm),
[MATRIX INPUT](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-matrixinput.htm),
and [Using Matrices](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEANLRefMap/simaanl-c-definingmatrices.htm).
