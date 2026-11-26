# Isolate triangulation error in interpolated-wall Stokes drag

Coordinate with [#798](https://github.com/tangency-labs/fieldcraft/issues/798). The existing
interpolated-wall test explicitly attaches an exact sphere to a triangle mesh. The LES cylinder
driver supplies a triangulated obstacle, so a passing exact-sphere check leaves a geometry-path
difference.
[Exact sphere setup](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/interpolated.rs#L179-L208) ·
[LES cylinder construction](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L200-L223)

## Suggested change

- Extend `sphere_drag` with a centre offset and a choice of exact-surface metadata. Run the same
  generated triangles with and without that metadata. At each comparison keep the sphere centre,
  nominal radius, box size, viscosity, forcing and wall rule identical. Record the mesh tessellation
  separately from lattice resolution. Refine the tessellation at fixed lattice spacing to isolate
  faceting.
- Reuse the R = 3, 6, 9 controls, the steady force-balance check and Hasimoto comparison. Use
  matching phases from the existing #798 work where available, and one finer radius to check the
  trend. Keep the established TRT case as the control. A low-τ HRR run is an extension once its
  creeping-flow and steady-state requirements fit the work budget.
  [Drag measurement](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/interpolated.rs#L209-L232) ·
  [Radius and viscosity controls](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/interpolated.rs#L245-L270) ·
  [Halfway study in PR #990](https://github.com/tangency-labs/fieldcraft/pull/990)

## Acceptance

- [ ] Report drag error by geometry path, lattice radius, tessellation and sub-cell phase at fixed
  nominal solid fraction, including Re, maximum Mach and force-balance residual. Mark any
  unconverged rows.
- [ ] Preserve the existing exact-sphere regression. Show the effect of varying tessellation at
  fixed lattice spacing separately from that of refining the lattice at fixed geometry.
