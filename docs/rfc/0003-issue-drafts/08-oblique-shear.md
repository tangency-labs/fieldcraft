# Measure D3Q19 HRR shear-wave damping on a moving stream

Diagnostic for [#1154](https://github.com/tangency-labs/fieldcraft/issues/1154). A maximum growth
rate helps identify instability. A stable collision can also damp an under-resolved physical shear
wave excessively, depending on its direction. Wissocq et al. studied two-dimensional models. The
proposed measurement is for FieldCraft's D3Q19 operator.
[Wissocq, Coreixas & Boussuge](https://arxiv.org/abs/2006.07353)

## Suggested change

- Extend the bounded CPU shear-wave test to a periodic 3D box with uniform background velocity.
  Choose a lattice-periodic wavevector and initialize a small perturbation perpendicular to it.
  Project the transverse velocity onto **both** sine and cosine and use the magnitude, so advection
  through the box does not appear as damping. Reuse the existing axial, zero-stream case as a
  control.
  [Current axial projection and step cap](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/hrr.rs#L295-L337)
- For example, in a cubic periodic box, integer mode vectors (5, 0, 0) and (3, 4, 0) have equal
  length and admit the same transverse polarization along z. Use a box large enough for both modes
  to be resolved, and specify the angle between the background stream and each wavevector. Compare
  their measured rates at matched physical wavenumber magnitude. Include the full finite-difference
  strain stencil if predicting HRR amplification. The local-collision TRT matrix alone omits this
  neighbour contribution. Start with fixed σ and the speeds and relaxation times relevant to the
  planned levels.
  [HRR stencil](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/src/hrr.rs#L26-L43) ·
  [TRT analysis for comparison](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/trt_stability.rs#L59-L98)

## Acceptance

- [ ] Recover the existing axial control and the long-wave continuum decay rate ν|k|² within a
  justified tolerance. Show the projection is insensitive to a change of initial phase.
- [ ] Publish measured decay or growth, |k|, direction, stream speed, τ and σ per case, including
  early stops. Keep physical shear damping distinct from any other growing mode. The results inform
  the setup-warning decision and a possible collision replacement.
