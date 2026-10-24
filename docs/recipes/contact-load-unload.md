# Contact loading and unloading

Run the original [command recipe](contact-load-unload.json) in a fresh model:

```sh
cargo run -p fieldcraft-cli -- run docs/recipes/contact-load-unload.json
```

The same array can be passed as `commands` to the desktop control channel's `engine.script`.
Coordinates are in mm, forces in N and pressure in MPa, set by the recipe's `model.units` command.
The model has two touching 10 mm cubes, each one hex8 element. The entire lower cube is fixed;
the upper cube's lateral translations are held. Its top has a prescribed vertical displacement
of -0.1 mm. Frictionless mortar contact joins the touching faces; both cubes have isotropic
elasticity with E = 1000 MPa and Poisson ratio zero.

The two nonlinear static steps start independently:

- **Press** ramps the prescribed displacement from zero to its full value. The four secondary
  nodes close, with pressure about 10 MPa. The small remaining penetration follows the unchanged
  default augmented-contact convergence criterion.
- **Cycle** follows load factors 0 → 0.5 → -0.05 at step fractions 0, 0.5 and 1. The final upper
  displacement and open contact gap are 0.005 mm, with zero pressure. `solve.history` returns the
  accepted frames' load factors, mean top displacement and reaction sum. At the fully unloaded
  state, a relative reaction-balance warning can report tiny absolute roundoff against zero load.

These are static equilibria, not a physical transient. The recipe finishes with the Cycle contact
gap contour and a deformation scale of 1000 so the small opening is visible. Use scale 1 for the
actual geometry; change the view without changing the solution through `results.deformed`.

This is a working contact workflow, not a naturally chattering case. The deterministic solver
regression `a_committed_contact_increment_does_not_carry_its_chatter_band_into_the_next` separately
seeds the existing temporary band state, exercises a failed retry and contact augmentations,
commits the accepted increment, then checks the next predictor/Newton solve and a near-open
finite-difference tangent. No production state-injection command is provided.
