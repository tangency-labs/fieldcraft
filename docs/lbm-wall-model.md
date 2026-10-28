# Wall-modelled large-eddy simulation

A fluid step with Smagorinsky's model can model the wall instead of resolving it. Each fluid cell
next to a no-slip wall is rebuilt from a law of the wall, so a grid whose first cell lies in the log
layer still carries the wall's shear stress. Without a wall model such a grid has too little
resolution for the bounce-back wall to carry it.

The method is the lattice Boltzmann wall model of Malaspinas & Sagaut (2014) with Musker's
boundary-layer profile. It runs with the multiple-relaxation-time collision, on the CPU and the
GPU.

## Synopsis

```json
{"name": "Flow", "procedure": {
  "type": "fluid", "domain": "Channel", "end_time": 2.0,
  "collision": "mrt",
  "turbulence": {"type": "smagorinsky", "constant": 0.1, "wall_model": "musker"}
}}
```

`model.step.create` takes these fields. `solve.fluid` runs the step with `"backend": "cpu"` or
`"backend": "gpu"`.

## Parameters

| Field | Value | Meaning |
|---|---|---|
| `collision` | `"mrt"` | Multiple relaxation times (d'Humières et al. 2002). The wall model is validated with this collision only. |
| `turbulence.type` | `"smagorinsky"` | Required: the wall model is part of Smagorinsky's model. |
| `turbulence.constant` | 0.08–0.1 | Smagorinsky's C. The validation used 0.08 at Re_τ 1000 and 0.1 at Re_τ 2000. |
| `turbulence.wall_model` | `"musker"` | Rebuild the cells next to no-slip walls from Musker's profile. Omitted: the bounce-back alone. |
| `turbulence.van_driest` | omit | Van Driest damping is for wall-resolved grids. Don't combine it with a wall model. |

With the wall model, the fluid's own relaxation time may go down to τ = 0.5001. A grid that needs
a lower τ is refused with the value it needs.

## When to use it

- **Use it** when the first fluid cell's centre lies in the log layer, y⁺ ≈ 30–100. Typical
  cases are internal flows at Re_τ 10³ and above: ducts, heat-sink channels and fin passages.
- **Don't use it** for a wall-resolved grid, with the first cell at y⁺ ≲ 2. Use Van Driest damping
  or no wall treatment there.
- **Don't use it in the buffer layer** (y⁺ 5–30). The method is not designed or validated for a
  first cell there.

y⁺ = y u_τ / ν, where y is the distance of the cell centre from the wall, half a cell for a plane
wall on the grid.

## Method

Before each collision, for every wall cell x₁ (a fluid cell with a no-slip bounce-back link) and the
cell x₂ one link further along the wall normal:

1. **Friction velocity.** The velocity at x₂, relative to the wall, gives the direction e parallel
   to the wall and the speed U₂. u_τ is the root of U₂ = ũ(y₂; u_τ), where ũ is Musker's profile
   and y₂ = 3 y₁. The profile covers the viscous sublayer, the buffer layer and the log law in one
   expression. Each root is solved by Newton's method, guarded by bisection.
2. **Velocity.** x₁'s velocity is set to ũ(y₁; u_τ) e.
3. **Density.** x₁'s density is the one its known populations give with that velocity, by the
   Zou–He relation.
4. **Stress.** x₁'s non-equilibrium stress is the one its known populations carry. The unknown
   populations, those streaming in from the wall, are taken from their opposites.
5. **Reconstruction.** All of x₁'s populations are rebuilt as the equilibrium at that density and
   velocity, plus the regularised non-equilibrium part of that stress.
6. **Relaxation.** x₁ relaxes with ν + ν_t, where ν_t is a damped mixing-length viscosity from the
   profile's shear at y₁: (κ y₁ (1 − e^(−y₁⁺/A⁺)))² |∂ũ/∂y|, with κ = 0.384 and A⁺ = 26. Smagorinsky's
   viscosity applies everywhere else.

**Collision.** The D3Q19 multiple-relaxation-time operator maps the populations to 19 orthogonal
moments, each relaxed at its own rate:

| Moments | Rate |
|---|---|
| energy | 1.19 |
| energy square, quartic modes | 1.4 |
| energy flux | 1.2 |
| cubic modes | 1.98 |
| the five stresses | 1 / (3 (ν + ν_t) + ½), per cell |

Guo's forcing enters through the moments. The rates are those Malaspinas & Sagaut use, from
d'Humières et al.'s stability analysis.

**Walls.** The wall distance and normal come from each wall cell's bounce-back links: exact for
plane walls, and within a fraction of a cell for curved ones. A moving wall is taken relative to its
own velocity. A wall cell whose x₂ is not a fluid cell inside the domain keeps the bounce-back.

**Reported wall force.** The wall force that `solve.fluid` reports is the momentum the flow
transfers to the wall, the drag the simulation actually carries. It is not the stress the wall
model infers at x₂. The two differ on coarse grids, as the validation below shows.

The force has two parts:
- the bounce-back's momentum exchange on the wall cell's links, after the rebuilt cell's
  collision;
- the momentum the reconstruction (step 5) takes from the cell: the momentum of its streamed
  populations less that of its rebuilt ones.

Together they close the fluid's momentum budget every step, so at steady state the walls carry
the body force exactly. A wall cell on an edge or corner, with links into more than one wall,
puts all of the reconstruction's share into the wall of its nearest link.

The density the reconstruction sets (step 3) does not conserve mass exactly. In a coarse
channel the total mass drifts by 1e-7 to 1e-5 over the first 40 steps.

**Implementation.** `crates/lbm/src/wall.rs` and `crates/lbm/src/mrt.rs` (CPU, f64), and
`crates/lbm-gpu/src/lbm.wgsl` (GPU, f32). Both are written from the papers.

## Validation

**Case.** Periodic turbulent channel flow between two plane walls, half-height δ. The flow rate is
held constant by a body force updated every step, targeting the bulk velocity of Dean's
correlation, as in Malaspinas & Sagaut (2014).

| Item | Setting |
|---|---|
| Reynolds number | Re_τ 1000 and 2000 |
| Grid | N = 10, 20 and 30 cells per δ |
| Smagorinsky C | 0.08 at Re_τ 1000, 0.1 at Re_τ 2000 |
| Precision | f32 on the GPU (NVIDIA GB10) |
| Start-up discarded | 20 δ/u_τ |
| Statistics averaged | ≥ 50 δ/u_τ |

**References.**
- Lee & Moser (2015) DNS, for the mean profile and the skin friction C_f.
- Dean's correlation, for C_f at the run's bulk Reynolds number.

**Acceptance**, fixed before the runs: C_f and the log-layer velocity within 5 % of Lee & Moser.

**Two measures of C_f.**
- **Force balance:** the drag the flow carries, which equals the body force at steady state. Users
  get this one.
- **Wall-model u_τ:** the friction velocity the wall model infers at x₂. Malaspinas & Sagaut's
  results are reproduced with this measure.

| Re_τ | N | first cell y⁺ | C_f, force balance: vs Lee & Moser / vs Dean | C_f, wall-model u_τ: vs Dean | profile error E | log layer vs Lee & Moser |
|---|---|---|---|---|---|---|
| 1000 | 10 | 46 | −11.4 % / −14.5 % | +1.3 % | 2.65 | 7.6 % |
| 1000 | 20 | 23 | −8.9 % / −12.1 % | −0.9 % | 1.21 | 6.7 % |
| 1000 | 30 | 16 | −8.7 % / −11.9 % | +1.4 % | 0.99 | 8.3 % |
| 2000 | 10 | 100 | refused: needs τ = 0.50007 | — | — | — |
| 2000 | 20 | 52 | +10.2 % / +8.9 % | −3.1 % | 1.01 | 9.0 % |
| 2000 | 30 | 33 | **−0.7 % / −1.8 %** | −1.8 % | 0.29 | **3.3 %** |

E is the profile error measure defined by Malaspinas & Sagaut (2014). Their values with Musker's
profile are 2.11, 1.07 and 0.83 at Re_τ 950, and 1.15 and 0.88 at Re_τ 2000, for N = 10, 20 and 30.
At Re_τ 1000 the 20- and 30-cell grids put the first cell below y⁺ 30, outside the range the method
is meant for. They're listed for completeness.

A run of Re_τ 1000 at N = 10 on the CPU in f64 gives a force-balance C_f of −14.8 % against Dean,
against the GPU's −14.5 %.

**Results.**
- **The implementation reproduces Malaspinas & Sagaut (2014).** C_f by the wall model's u_τ is
  within the paper's 0.3–4.2 % of Dean on every grid, and E is close to the paper's values.
- **The drag the flow carries meets the acceptance at Re_τ 2000, N = 30.** There the second cell
  sits at y⁺ ≈ 100, well inside the log layer, and the inferred and carried u_τ agree.
- **On coarser grids the drag is off by 9–15 %,** with a sign that depends on Re_τ and N. The
  velocity the LES gives at x₂ is under-resolved, so the stress the model infers differs from the
  one the flow carries, the log-layer mismatch described by Kawai & Larsson (2012).

**Reproduce.**

```sh
FIELDCRAFT_PAPER_CHANNEL=1 cargo test --release -p fieldcraft-lbm-gpu --test verification \
  -- --ignored study_malaspinas_and_sagaut_s_channel_against_lee_and_moser --nocapture
```

`FIELDCRAFT_PAPER_CHANNEL=2000:30` runs one grid. Add `FIELDCRAFT_PAPER_CHANNEL_CPU=1` to run it on
the CPU in f64. The runs take 2 minutes to an hour each on a GB10.

## Limitations

- **Drag accuracy depends on the grid.** Expect skin friction and pressure drop within 5 % only
  near the validated set-up, with about 30 cells per δ and the second cell well inside the log
  layer. Plan for ±15 % on coarser grids.
- **Equilibrium boundary layers.** The profile assumes an attached, equilibrium boundary layer.
  Separation, reattachment and strong pressure gradients are not validated.
- **Smooth walls only.** No roughness.
- **No wall heat flux model.** The thermal lattice keeps its bounce-back walls, so wall heat transfer
  on a wall-modelled grid is not validated.
- **Validated with `"mrt"` only.** Other collisions accept the wall model but are not validated
  with it.

## References

- O. Malaspinas, P. Sagaut. *Wall model for large-eddy simulation based on the lattice Boltzmann
  method.* J. Comput. Phys. 275:25–40, 2014.
  doi:[10.1016/j.jcp.2014.06.020](https://doi.org/10.1016/j.jcp.2014.06.020)
- D. d'Humières, I. Ginzburg, M. Krafczyk, P. Lallemand, L.-S. Luo. *Multiple-relaxation-time
  lattice Boltzmann models in three dimensions.* Phil. Trans. R. Soc. Lond. A 360:437–451, 2002.
  doi:[10.1098/rsta.2001.0955](https://doi.org/10.1098/rsta.2001.0955)
- A. J. Musker. *Explicit expression for the smooth wall velocity distribution in a turbulent
  boundary layer.* AIAA J. 17(6):655–657, 1979.
  doi:[10.2514/3.61193](https://doi.org/10.2514/3.61193)
- M. Lee, R. D. Moser. *Direct numerical simulation of turbulent channel flow up to Re_τ ≈ 5200.*
  J. Fluid Mech. 774:395–415, 2015.
  doi:[10.1017/jfm.2015.268](https://doi.org/10.1017/jfm.2015.268)
- R. B. Dean. *Reynolds number dependence of skin friction and other bulk flow variables in
  two-dimensional rectangular duct flow.* J. Fluids Eng. 100(2):215–223, 1978.
  doi:[10.1115/1.3448633](https://doi.org/10.1115/1.3448633)
- S. Kawai, J. Larsson. *Wall-modeling in large eddy simulation: length scales, grid resolution,
  and accuracy.* Phys. Fluids 24:015105, 2012.
  doi:[10.1063/1.3678331](https://doi.org/10.1063/1.3678331)
- U. Piomelli, E. Balaras. *Wall-layer models for large-eddy simulations.* Annu. Rev. Fluid Mech.
  34:349–374, 2002.
  doi:[10.1146/annurev.fluid.34.082901.144919](https://doi.org/10.1146/annurev.fluid.34.082901.144919)

See also [the collision operators](manual/fluid/collision.md), [verification and validation](verification.md)
and [the solver literature](solvers/literature.md), entries C5.1–C5.4.
